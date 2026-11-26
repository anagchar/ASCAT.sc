## Allele-specific copy number from haplotagged reads (HP tag), per bin.
## The inputs are the per-bin
## hap1/hap2 read counts that getCoverageTrack.10XBAM(hp_tag="HP") adds
## to each cell's track.

## Unfolded BAF of one segment: pooled hap1 / (hap1 + hap2), p in [0, 1].
## Only meaningful when the HP labels are consistent over the segment (chromosome-scale phasing).
fitBinom.hap.noswitch <- function(h1, h2, values=seq(0, 1, .001))
{
    ok <- !is.na(h1) & !is.na(h2)
    n1 <- sum(h1[ok])
    n <- n1 + sum(h2[ok])
    if (n == 0) return(rep(NA_real_, 3))
    llh <- dbinom(n1, size=n, prob=values, log=TRUE)
    llh <- llh - max(llh)
    post <- exp(llh) / sum(exp(llh))
    cs <- cumsum(post)
    c(values[which(cs > .05)[1]], values[which.max(llh)], values[which(cs > .95)[1]])
}

## purity/ploidy grid search on segment BAF + logR (sizes as weights); nA/nB from BAF and total CN
.hp_searchGrid <- function(baf, logr, sizes, purs=seq(.90,.99,.01), ploidies=seq(2,6,.01))
{
    getNANB <- function(baf, logr, purity, ploidy)
    {
        K <- (2^logr*ploidy-2*(1-purity))/purity
        Na <- (1-purity-baf*(2-2*purity)-(purity*baf-purity)*K)/purity
        Nb <- K-Na
        list(Na=Na, Nb=Nb)
    }
    errors <- function(purity, ploidy)
    {
        nanb <- getNANB(baf, logr, purity, ploidy)
        sum(((nanb$Na-round(nanb$Na))^2+(nanb$Nb-round(nanb$Nb))^2)*(sizes/sum(sizes)))
    }
    errs <- matrix(sapply(ploidies, function(psi) sapply(purs, function(rho) errors(rho, psi))),
                   nrow=length(purs), ncol=length(ploidies), dimnames=list(purs, ploidies))
    mins <- arrayInd(which.min(errs), dim(errs))
    purity <- purs[mins[1]]
    ploidy <- ploidies[mins[2]]
    list(purity=purity, ploidy=ploidy, profile=getNANB(baf, logr, purity, ploidy), errs=errs)
}

.hp_getNANBfromTot <- function(baf, tot, purity, retNa=TRUE)
{
    Na <- round(baf*(purity*tot+(1-purity)*2)-1+purity)
    Nb <- round(tot-Na)
    if(retNa) return(Na)
    Nb
}

.hp_chr <- function(x)
{
    x <- gsub("chr", "", as.character(x))
    x[x == "23"] <- "X"
    x[x == "24"] <- "Y"
    x
}

## ---------------------------------------------------------------------------------------------
## Allele-specific calls: state likelihoods straight from the per-bin hap1/hap2 counts.
##   per bin j of a segment, state nMaj+nMin with expected major-haplotype fraction
##       f = (purity*(nMaj*(1-eps) + nMin*eps) + (1-purity)) / (purity*total + 2*(1-purity))
##   L_j = 0.5*BB(h1_j | n_j, f, rho) + 0.5*BB(h1_j | n_j, 1-f, rho)
##     - the 0.5/0.5 mixture per bin: HP labels may switch between bins (phase switch errors)
##     - eps: fraction of reads tagged with the wrong haplotype (LOH/1 copy gives f = 1-eps, not 1)
##     - rho: beta-binomial overdispersion (single-cell amplification: reads in a bin are not independent)
##   segment log-likelihood = sum_j log L_j
## call_rule = "logr_first" (default): logR picks the total, the haplotypes pick the split
##   1. eps and rho per cell on a grid, maximising sum over segments of the best state's llh (any total:
##      restricting to logR-compatible states lets logR/BAF-discordant segments inflate eps)
##   2. BAF_fit = the segment's observed major-haplotype fraction under the same model (beta-binomial,
##      switch-tolerant, this cell's rho), i.e. on the same scale as the expected f of a state
##   3. logR-compatible states: |total - ntot| <= logr_tol (none -> "logR_discordant")
##   4. call = logR-compatible state with the highest llh, accepted if |f_state - BAF_fit| <= baf_tol or its
##      llh is within delta of the best state of any total; otherwise "BAF_discordant"
##   The check is on the BAF scale on purpose: a llh tolerance shrinks with depth, so at high depth a
##   0.02 error in eps/purity would reject the right state (e.g. 4+1 at logR CN 5.07, BAF 0.735, eps 0.15).
## call_rule = "baf_first": BAF-compatible states = llh >= best llh - delta; among them the total closest
##   to ntot (ties: higher llh); |ntot - total| > logr_tol -> "logR_discordant".
## No call -> 1+1 (1+0 on chrX/Y of a male cell).
## hap_llr (fitted segments): evidence for allelic imbalance from the haplotype counts alone, log-likelihood
##   of the best folded fraction (BAF_fit) minus that of 0.5 (balanced), same per-bin model and rho. It does not
##   depend on eps (BAF_fit vs the state's expected fraction is checked by baf_tol). ~0 when the counts cannot
##   tell (few reads) or the segment is balanced; >= log(20) is used in the plots to call an imbalanced state
##   (1+0, 2+0, 2+1, ...) "verified", i.e. logR and haplotypes agree, as in hap_stats.R.
## ---------------------------------------------------------------------------------------------
.hp_bbll <- function(k, n, f, rho)
{
    if(rho < 1e-6) return(dbinom(k, n, f, log=TRUE))
    a <- f*(1-rho)/rho
    b <- (1-f)*(1-rho)/rho
    lchoose(n, k) + lbeta(k+a, n-k+b) - lbeta(a, b)
}

## Per-bin log-likelihood matrix, bins x f: log(0.5*BB(h1 | n, f, rho) + 0.5*BB(h1 | n, 1-f, rho)).
## Summing its rows over the bins of a segment gives the segment log-likelihood for each f.
.hp_binll <- function(h1, n, f, rho)
{
    ## sparse counts: many bins share the same (h1, n); compute each distinct pair once
    key <- paste(h1, n)
    u <- !duplicated(key)
    k1 <- h1[u]; n1 <- n[u]
    nu <- length(k1)
    k <- rep(k1, times=length(f))
    nn <- rep(n1, times=length(f))
    ff <- rep(f, each=nu)
    a <- .hp_bbll(k, nn, ff, rho)
    b <- .hp_bbll(k, nn, 1-ff, rho)
    m <- pmax(a, b)
    matrix(m + log(.5*exp(a-m) + .5*exp(b-m)), nrow=nu, ncol=length(f))[match(key, key[u]), , drop=FALSE]
}

## Segment sums of a bins x f matrix: segments (one per element of `rows`) x f
.hp_segsum <- function(M, rows)
    do.call("rbind", lapply(rows, function(r) colSums(M[r, , drop=FALSE])))

## Folded fraction of each segment on a grid (flat prior): c(q05, best, q95) per row of a segments x fgrid llh
.hp_fitci <- function(SL, fgrid)
    t(apply(SL, 1, function(llh)
    {
        llh <- llh - max(llh)
        cs <- cumsum(exp(llh) / sum(exp(llh)))
        c(fgrid[which(cs > .05)[1]], fgrid[which.max(llh)], fgrid[which(cs > .95)[1]])
    }))

## segbins: list (one element per segment) of data.frames with hap1 and depth of its bins (depth > 0);
## use: segments to fit. eps, rho given: no grid search, the calls use these values.
## Returns list(calls = data.frame, one row per segment;
##              baf = matrix q05 / BAF / q95 of BAF_fit, NA for unused segments).
callAS_HPmodel <- function(segbins, use, ntot, chr, purity=1, ismale=FALSE,
                           logr_tol=0.5, max_cn=6, delta=2,
                           baf_tol=0.05, call_rule=c("logr_first", "baf_first"),
                           eps_grid=c(.01,.02,.03,.05,.075,.1,.125,.15,.2,.25),
                           rho_grid=c(0,.01,.02,.05,.1,.2),
                           eps=NULL, rho=NULL)
{
    if(!is.null(eps) && !is.null(rho) && !is.na(eps) && !is.na(rho))
    {
        eps_grid <- eps; rho_grid <- rho
    }
    call_rule <- match.arg(call_rule)
    states <- do.call("rbind", lapply(1:max_cn, function(tot)
        cbind(nMaj=ceiling(tot/2):tot, nMin=tot-(ceiling(tot/2):tot))))
    tots <- states[,"nMaj"] + states[,"nMin"]
    fexp <- function(eps)
        pmin(pmax((purity*(states[,"nMaj"]*(1-eps) + states[,"nMin"]*eps) + 1 - purity) /
                  (purity*tots + 2*(1-purity)), 1e-3), 1-1e-3)
    idx <- which(use)
    ## logR-compatible states per segment (logical matrix segments x states)
    compat <- outer(ntot, tots, function(n, t) !is.na(n) & abs(t - n) <= logr_tol)

    ## all bins of the fitted segments stacked; rows[[r]] = the bins of segment idx[r]
    h1 <- unlist(lapply(segbins[idx], function(b) b$hap1))
    dp <- unlist(lapply(segbins[idx], function(b) b$depth))
    nb <- vapply(segbins[idx], nrow, integer(1))
    rows <- split(seq_along(h1), rep(seq_along(idx), nb))
    ## expected f of every state for every eps (rounded as the distinct values the llh is computed at)
    fr <- lapply(eps_grid, function(e) round(fexp(e), 8))
    uf <- unique(unlist(fr))

    eps <- NA_real_; rho <- NA_real_; LL <- NULL
    if(length(idx) > 0)
    {
        ## segments x distinct f, once per rho; a grid point (eps, rho) picks its states' columns
        SF <- lapply(rho_grid, function(r) .hp_segsum(.hp_binll(h1, dp, uf, r), rows))
        grid <- expand.grid(eps=seq_along(eps_grid), rho=seq_along(rho_grid))
        obj <- sapply(seq_len(nrow(grid)), function(g)
            sum(apply(SF[[grid$rho[g]]][, match(fr[[grid$eps[g]]], uf), drop=FALSE], 1, max)))
        best <- which.max(obj)                   # ties -> smallest eps / rho (grid order)
        eps <- eps_grid[grid$eps[best]]; rho <- rho_grid[grid$rho[best]]
        LL <- SF[[grid$rho[best]]][, match(fr[[grid$eps[best]]], uf), drop=FALSE]
    }
    fbest <- fexp(if(is.na(eps)) 0 else eps)
    ## observed major-haplotype fraction per segment under the fitted model (same scale as fbest)
    fgrid <- seq(.5, .999, .001)
    baf <- matrix(NA_real_, length(ntot), 3, dimnames=list(NULL, c("q05", "BAF", "q95")))
    hap_llr <- rep(NA_real_, length(ntot))
    if(length(idx) > 0)
    {
        SL <- .hp_segsum(.hp_binll(h1, dp, fgrid, rho), rows)
        baf[idx, ] <- .hp_fitci(SL, fgrid)
        hap_llr[idx] <- apply(SL, 1, max) - SL[, 1]
    }
    baf_fit <- baf[, "BAF"]

    out <- lapply(seq_along(ntot), function(i)
    {
        sexchr <- ismale && .hp_chr(chr[i]) %in% c("X", "Y")
        normal <- if(sexchr) c(1, 0) else c(1, 1)
        nocall <- function(status, fb=NA) list(normal[1], normal[2], NA, fb, status)
        r <- match(i, idx)
        if(is.na(r) || is.na(ntot[i])) return(nocall("no_BAF"))
        ll <- LL[r, ]
        fb <- fbest[which.max(ll)]
        if(call_rule == "logr_first")
        {
            cand <- which(compat[i, ])
            if(!length(cand)) return(nocall("logR_discordant", fb))
            sel <- cand[which.max(ll[cand])]
            baf_ok <- abs(fbest[sel] - baf_fit[i]) <= baf_tol
            if(!baf_ok && ll[sel] < max(ll) - delta) return(nocall("BAF_discordant", fb))
            return(list(states[sel,"nMaj"], states[sel,"nMin"], fbest[sel], fb, "called"))
        }
        ok <- which(ll >= max(ll) - delta)
        sel <- ok[order(abs(tots[ok] - ntot[i]), -ll[ok])[1]]
        if(abs(tots[sel] - ntot[i]) > logr_tol) return(nocall("logR_discordant", fb))
        list(states[sel,"nMaj"], states[sel,"nMin"], fbest[sel], fb, "called")
    })
    g <- function(k) sapply(out, function(o) o[[k]])
    calls <- data.frame(nMaj_call=as.numeric(g(1)),
                        nMin_call=as.numeric(g(2)),
                        total_call=as.numeric(g(1)) + as.numeric(g(2)),
                        BAF_expected_call=as.numeric(g(3)),
                        BAF_model_best=as.numeric(g(4)),
                        BAF_fit=baf_fit,
                        call_status=as.character(g(5)),
                        tag_error=eps,
                        rho=rho,
                        hap_llr=hap_llr,
                        stringsAsFactors=FALSE)
    list(calls=calls, baf=baf)
}

## BAF change-points inside one logR segment, on the per-bin folded beta-binomial.
## h1, n: hap1 and depth of the segment's bins with reads, in genomic order; rho: the cell's overdispersion.
## Score of a stretch of bins = max_f sum_j log L_j(f), f on fgrid (folded fraction, 0.5-1).
##   method = "binary":   split at the bin boundary that maximises LL(left) + LL(right) - LL(whole)
##   method = "circular": as circular binary segmentation (DNAcopy): the best middle piece a..b against the
##                        rest of the segment (both flanks share one f), LL(middle) + LL(rest) - LL(whole);
##                        a single cut is the case of an empty flank. Finds a short change inside a long
##                        segment, which a single cut dilutes.
## Accepted if the gain >= llr and every new piece has >= min_bins bins and >= min_reads reads; the pieces
## are then searched again. Cumulative sums of the bins x f matrix give each stretch's score by subtraction.
## Returns the indices of the last bin of every piece but the last (integer(0): no split).
.hp_cosegment <- function(h1, n, rho, llr=10, min_bins=5, min_reads=20, fgrid=seq(.5, .999, .001),
                          method=c("circular", "binary"))
{
    method <- match.arg(method)
    if(length(h1) < 2*min_bins) return(integer(0))
    M <- .hp_binll(h1, n, fgrid, rho)
    rowmax <- function(X) X[cbind(seq_len(nrow(X)), max.col(X, "first"))]
    cuts <- integer(0)
    search <- function(lo, hi)
    {
        nb <- hi - lo + 1
        if(nb < 2*min_bins) return(invisible())
        C <- rbind(0, apply(M[lo:hi, , drop=FALSE], 2, cumsum))   # C[k+1, ] = sum of bins 1..k
        D <- c(0, cumsum(n[lo:hi]))
        tot <- C[nb+1, ]
        ok <- function(k) k >= min_bins                            # piece size in bins
        best <- c(gain=-Inf, a=NA, b=NA)
        ## middle piece = bins a..b (1-based); flanks 1..a-1 and b+1..nb: empty or a full piece each
        starts <- if(method == "binary") 1 else seq_len(nb)
        for(a in starts)
        {
            if(a > 1 && (!ok(a - 1) || D[a] < min_reads)) next
            b <- if(method == "binary" || a > 1) (a + min_bins - 1):nb else (a + min_bins - 1):(nb - min_bins)
            b <- b[b <= nb]
            if(method == "binary") b <- b[b <= nb - min_bins]
            right <- nb - b
            b <- b[(right == 0 | (right >= min_bins & D[nb+1] - D[b+1] >= min_reads)) &
                   D[b+1] - D[a] >= min_reads & !(a == 1 & right == 0)]
            if(!length(b)) next
            mid <- C[b+1, , drop=FALSE] - matrix(C[a, ], length(b), length(tot), byrow=TRUE)
            rest <- matrix(tot, length(b), length(tot), byrow=TRUE) - mid
            g <- rowmax(mid) + rowmax(rest) - max(tot)
            if(max(g) > best["gain"]) best <- c(gain=max(g), a=a, b=b[which.max(g)])
        }
        if(!is.finite(best["gain"]) || best["gain"] < llr) return(invisible())
        new <- c(if(best["a"] > 1) lo + best["a"] - 2, if(best["b"] < nb) lo + best["b"] - 1)
        cuts <<- c(cuts, new)
        bnd <- c(lo - 1, sort(new), hi)
        for(k in seq_len(length(bnd) - 1)) search(bnd[k] + 1, bnd[k + 1])
    }
    search(1, length(h1))
    sort(unname(cuts))
}

## One cell: allele-specific calls from the hap1/hap2 counts of the bins inside each total-CN segment
## (callAS_HPmodel), then the purity/ploidy search, nA/nB and ntot_* per segment, so that
## nprof.fixed (columns BAF, ntot_fixed) feeds getAS_CNA_smoothed() unchanged.
## BAF / q05 / q95 = the segment's folded major-haplotype fraction under the model (BAF_fit) and its 5-95%
## interval. BAF_noswitch / q05_noswitch / q95_noswitch = pooled hap1 / (hap1 + hap2) (unfolded).
## coseg = TRUE: eps and rho are fitted on the logR segments, then each logR segment is split at BAF
## change-points (.hp_cosegment) and the calls are made on the pieces with the same eps and rho. The pieces
## keep the logR values (and so the total copy number) of their logR segment; an imbalanced piece at CN 2
## is then called 2+0 (CN-LOH). seg_source: "logR" (segment as in the total-CN profile) or "BAF" (piece of a
## logR segment split at BAF change-points).
## Extra columns at the end: nBins_hap, hap1, hap2, num.mark, calls, seg_source. Per-bin counts in $bin_baf.
getAS_CNA_HP_sample <- function(track,
                                profile,
                                purity,
                                ploidy,
                                purs,
                                ploidies,
                                min_reads=20,
                                min_bins=2,
                                ismale=FALSE,
                                baf_tol=0.05,
                                logr_tol=0.5,
                                max_cn=6,
                                delta=2,
                                call_rule="logr_first",
                                eps_grid=NULL,
                                rho_grid=NULL,
                                coseg=TRUE,
                                coseg_llr=10,
                                coseg_min_bins=5,
                                coseg_method="circular",
                                cell_name=NULL)
{
    suppressPackageStartupMessages(require(GenomicRanges))
    if(!is.data.frame(profile)) stop("no total copy-number profile")
    if(!all(c("hap1", "hap2") %in% colnames(track[[1]])))
        stop("no hap1/hap2 columns in nlCTS.tumour - run run_sc_sequencing with hp_tag")

    bins <- do.call("rbind", lapply(track, function(df)
        data.frame(chr=.hp_chr(df$space), start=df$start, end=df$end,
                   hap1=df$hap1, hap2=df$hap2, stringsAsFactors=FALSE)))
    bins$depth <- bins$hap1 + bins$hap2
    bins$BAF_noswitch <- ifelse(bins$depth > 0, bins$hap1 / bins$depth, NA)
    if(!is.null(cell_name)) bins$cell <- cell_name

    nprof <- data.frame(chr=as.character(profile[,"chromosome"]),
                        startpos=as.numeric(profile[,"start"]),
                        endpos=as.numeric(profile[,"end"]),
                        total_copy_number=as.numeric(profile[,"total_copy_number"]),
                        total_copy_number_logr=as.numeric(profile[,"total_copy_number_logr"]),
                        logr=as.numeric(profile[,"logr"]),
                        logr.sd=as.numeric(profile[,"logr.sd"]),
                        fitted=as.numeric(profile[,"total_copy_number"]),
                        q05=as.numeric(NA), BAF=as.numeric(NA), q95=as.numeric(NA),
                        q05_noswitch=as.numeric(NA), BAF_noswitch=as.numeric(NA), q95_noswitch=as.numeric(NA),
                        stringsAsFactors=FALSE)
    nprof[nprof[,1]=="chr23",1] <- "chrX"
    nprof[nprof[,1]=="23",1] <- "X"
    nprof[nprof[,1]=="chr24",1] <- "chrY"
    nprof[nprof[,1]=="24",1] <- "Y"

    ## each bin goes to the segment that contains its midpoint
    grseg <- GRanges(.hp_chr(nprof$chr), IRanges(nprof$startpos, nprof$endpos))
    mid <- round((bins$start + bins$end) / 2)
    ovs <- findOverlaps(grseg, GRanges(bins$chr, IRanges(mid, mid)))
    qH <- queryHits(ovs)
    sH <- subjectHits(ovs)

    ## bins of each segment (all, and with haplotagged reads)
    segall <- vector("list", nrow(nprof))
    for(i in unique(qH)) segall[[i]] <- bins[sort(sH[qH == i]), ]
    nmark <- as.numeric(profile[,"num.mark"])
    seg_source <- rep("logR", nrow(nprof))

    segstats <- function()
    {
        segbins <<- lapply(segall, function(b) if(is.null(b)) NULL else b[b$depth > 0, ])
        nBins_hap <<- vapply(segbins, function(b) if(is.null(b)) 0 else nrow(b), numeric(1))
        hap1_seg <<- vapply(segbins, function(b) if(is.null(b)) 0 else sum(b$hap1), numeric(1))
        hap2_seg <<- vapply(segbins, function(b) if(is.null(b)) 0 else sum(b$hap2), numeric(1))
        use <<- nBins_hap >= min_bins & hap1_seg + hap2_seg >= min_reads
        nprof[, c("q05_noswitch","BAF_noswitch","q95_noswitch")] <<- NA_real_
        for(i in which(use))
            nprof[i, c("q05_noswitch","BAF_noswitch","q95_noswitch")] <<-
                fitBinom.hap.noswitch(segbins[[i]]$hap1, segbins[[i]]$hap2)
    }
    segbins <- nBins_hap <- hap1_seg <- hap2_seg <- use <- NULL
    segstats()

    if(sum(use) == 0) stop("no segment passed min_bins/min_reads - no haplotagged reads?")
    ## logR copy number at the cell's purity/ploidy (used by the calls as validation)
    nprof[,"ntot_fixed"] <- transform_bulk2tumour(nprof[,"logr"], purity, ploidy)

    ## state likelihoods from the per-bin hap1/hap2 counts; logR (ntot_fixed) as validation
    callfit <- function(...)
        callAS_HPmodel(segbins=segbins, use=use, ntot=nprof[,"ntot_fixed"], chr=nprof[,"chr"],
                       purity=purity, ismale=ismale, logr_tol=logr_tol, max_cn=max_cn, delta=delta,
                       baf_tol=baf_tol, call_rule=call_rule,
                       eps_grid=if(is.null(eps_grid)) eval(formals(callAS_HPmodel)$eps_grid) else eps_grid,
                       rho_grid=if(is.null(rho_grid)) eval(formals(callAS_HPmodel)$rho_grid) else rho_grid, ...)
    fit <- callfit()

    ## BAF change-points inside each logR segment, with the cell's rho; calls again on the pieces with the
    ## same eps / rho
    if(coseg && !is.na(fit$calls$rho[1]))
    {
        eps <- fit$calls$tag_error[1]
        rho <- fit$calls$rho[1]
        cuts <- lapply(seq_len(nrow(nprof)), function(i)
        {
            b <- segbins[[i]]
            if(is.null(b) || (ismale && .hp_chr(nprof$chr[i]) %in% c("X", "Y"))) return(integer(0))
            .hp_cosegment(b$hap1, b$depth, rho, llr=coseg_llr, min_bins=coseg_min_bins, min_reads=min_reads,
                          method=coseg_method)
        })
        if(any(lengths(cuts) > 0))
        {
            parent <- integer(0); starts <- ends <- nm <- numeric(0); pall <- list()
            for(i in seq_len(nrow(nprof)))
            {
                k <- cuts[[i]]
                if(!length(k))
                {
                    parent <- c(parent, i); starts <- c(starts, nprof$startpos[i]); ends <- c(ends, nprof$endpos[i])
                    nm <- c(nm, nmark[i]); pall <- c(pall, list(segall[[i]]))
                    next
                }
                b <- segbins[[i]]
                bnd <- b$start[k + 1]                  # first bin of each right-hand piece
                st <- c(nprof$startpos[i], bnd)
                en <- c(bnd - 1, nprof$endpos[i])
                a <- segall[[i]]
                pieces <- lapply(seq_along(st), function(p) a[a$start >= st[p] & a$start <= en[p], ])
                n_p <- vapply(pieces, nrow, numeric(1))
                n_p[length(n_p)] <- nmark[i] - sum(n_p[-length(n_p)])   # keeps the sum of num.mark
                parent <- c(parent, rep(i, length(st))); starts <- c(starts, st); ends <- c(ends, en)
                nm <- c(nm, n_p); pall <- c(pall, pieces)
            }
            split_parent <- parent %in% which(lengths(cuts) > 0)
            nprof <- nprof[parent, ]
            rownames(nprof) <- NULL
            nprof$startpos <- starts
            nprof$endpos <- ends
            nmark <- nm
            segall <- pall
            seg_source <- ifelse(split_parent, "BAF", "logR")
            segstats()
            fit <- callfit(eps=eps, rho=rho)
        }
    }
    calls <- fit$calls
    calls$call_rule <- call_rule
    calls$baf_mode <- "model"
    nprof[, c("q05","BAF","q95")] <- fit$baf

    nona <- !is.na(nprof[,"BAF"])
    exceptNA <- function(vec) { newvec <- rep(NA, length(nona)); newvec[nona] <- vec; newvec }
    sizes <- nprof[nona,"endpos"] - nprof[nona,"startpos"]
    sG <- .hp_searchGrid(nprof[nona,"BAF"], nprof[nona,"logr"], sizes, purs=purs, ploidies=ploidies)
    sG.fixed <- .hp_searchGrid(nprof[nona,"BAF"], nprof[nona,"logr"], sizes, purs=purity, ploidies=ploidy)
    nprof[,"ntot_free"] <- transform_bulk2tumour(nprof[,"logr"], sG$purity, sG$ploidy)
    nprof <- nprof[, c(setdiff(colnames(nprof), c("ntot_free", "ntot_fixed")), "ntot_free", "ntot_fixed")]

    build <- function(s)
        cbind(nprof,
              nA=.hp_getNANBfromTot(baf=nprof[,"BAF"], tot=nprof[,"fitted"], purity=purity),
              nB=.hp_getNANBfromTot(baf=nprof[,"BAF"], tot=nprof[,"fitted"], purity=purity, retNa=FALSE),
              nA_sc_raw=exceptNA(s$profile$Nb),
              nB_sc_raw=exceptNA(s$profile$Na),
              nA_sc=exceptNA(round(s$profile$Nb)),
              nB_sc=exceptNA(round(s$profile$Na)),
              nBins_hap=nBins_hap,
              hap1=hap1_seg,
              hap2=hap2_seg,
              num.mark=nmark,
              calls,
              seg_source=seg_source)
    list(nprof.free=build(sG),
         purity.free=sG$purity,
         ploidy.free=sG$ploidy,
         nprof.fixed=build(sG.fixed),
         purity.fixed=sG.fixed$purity,
         ploidy.fixed=sG.fixed$ploidy,
         bin_baf=bins)
}

## All cells: fills res$allProfiles_AS (one data.frame per cell, columns read by getAS_CNA_smoothed()
## and the plotting code), writes one table per cell to
## outdir/as_cna_profile_hp/ and all profiles to outdir/all_as_cna_profiles_hp_<projectname>.pdf.
## Cells that fail get an error string instead of a list (getAS_CNA_smoothed() skips them).
getAS_CNA_HP <- function(res,
                         purs,
                         ploidies,
                         projectname="project",
                         outdir="./",
                         min_reads=20,
                         min_bins=2,
                         baf_tol=0.05,
                         logr_tol=0.5,
                         max_cn=6,
                         delta=2,
                         call_rule="logr_first",
                         eps_grid=NULL,
                         rho_grid=NULL,
                         coseg=TRUE,
                         coseg_llr=10,
                         coseg_min_bins=5,
                         coseg_method="circular",
                         mc.cores=1)
{
    cell_names <- names(res$allTracks)
    sex <- rep(if(is.null(res$sex)) "female" else res$sex, length.out=length(cell_names))
    if(is.null(cell_names)) stop("res$allTracks must be named")
    if(!is.list(purs)) purs <- lapply(cell_names, function(x) purs)
    if(!is.list(ploidies)) ploidies <- lapply(cell_names, function(x) ploidies)
    refitted <- any(grepl("refitted", names(res)))

    print(paste0("## derive Allele-specific Profiles from HP-tag haplotype counts (call_rule=", call_rule,
                 ", coseg=", coseg, if(coseg) paste0(", method=", coseg_method, ", llr=", signif(coseg_llr, 3),
                                          ", min_bins=", coseg_min_bins), ")"))
    res$allProfiles_AS <- parallel::mclapply(setNames(seq_along(cell_names), cell_names), function(x)
    {
        cell <- cell_names[x]
        sol <- if(refitted) res$allSolutions.refitted.auto[[cell]] else res$allSolutions[[cell]]
        tryCatch(getAS_CNA_HP_sample(track=res$allTracks[[cell]]$nlCTS.tumour,
                                     profile=res$allProfiles[[cell]],
                                     purity=sol$purity,
                                     ploidy=sol$ploidy,
                                     purs=purs[[x]],
                                     ploidies=ploidies[[x]],
                                     min_reads=min_reads,
                                     min_bins=min_bins,
                                     ismale=sex[x] == "male",
                                     baf_tol=baf_tol,
                                     logr_tol=logr_tol,
                                     max_cn=max_cn,
                                     delta=delta,
                                     call_rule=call_rule,
                                     eps_grid=eps_grid,
                                     rho_grid=rho_grid,
                                     coseg=coseg,
                                     coseg_llr=coseg_llr,
                                     coseg_min_bins=coseg_min_bins,
                                     coseg_method=coseg_method,
                                     cell_name=cell),
                 error=function(e) paste0("getAS_CNA_HP failed for ", cell, ": ", conditionMessage(e)))
    }, mc.cores=mc.cores)

    failed <- !sapply(res$allProfiles_AS, function(x) is.list(x) && "nprof.fixed" %in% names(x))
    if(all(failed))
        stop("getAS_CNA_HP: no cell got an allele-specific profile, e.g. ", res$allProfiles_AS[[1]])
    if(any(failed))
        warning(sum(failed), " / ", length(failed), " cells without HP allele-specific profile, e.g. ",
                res$allProfiles_AS[failed][[1]])

    print("## write to disk and plot Allele-specific Profiles (HP)")
    as_cna_dir <- file.path(outdir, "as_cna_profile_hp")
    dir.create(as_cna_dir, showWarnings=FALSE, recursive=TRUE)
    pdf(file.path(outdir, paste0("all_as_cna_profiles_hp_", projectname, ".pdf")), width=15, height=5)
    for(x in which(!failed))
    {
        cell <- cell_names[x]
        try({
            plot_AS_profile(res$allProfiles_AS[[x]]$nprof.fixed)
            title(cell, cex=.5)
        })
        write.table(res$allProfiles_AS[[x]]$nprof.fixed, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE,
                    file=file.path(as_cna_dir, paste0("as_cna_profile_hp_", cell, ".txt")))
    }
    dev.off()
    res
}
