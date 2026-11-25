## Allele-specific copy number from haplotagged reads (HP tag), per bin.
## Counterpart of the SNP-based getAS_CNA(): the inputs are the per-bin
## hap1/hap2 read counts that getCoverageTrack.10XBAM(hp_tag="HP") adds
## to each cell's track.

## BAF of one segment from per-bin hap1/hap2 counts, robust to phase switch errors.
##
## Each block of `block` consecutive bins contributes
##     P(h1 | n, p) + P(h1 | n, 1 - p)
## i.e. the haplotype labels of a block are treated as exchangeable. A phase
## switch between blocks therefore cannot bias the estimate, but the labels
## (which haplotype is the amplified one) are not used. The result is the
## folded BAF in [0.5, 1].
##   block = 1  : fully switch-robust, least power (sparse single-cell counts)
##   block > 1  : more power per term; a switch inside a block dilutes it
## Returns c(lower 5%, best, upper 95%); NA if the segment has no phased reads.
fitBinom.hap <- function(h1,
                         h2,
                         block=1,
                         maxdepth=1000,
                         values=seq(.5, 1, .001))
{
    ok <- !is.na(h1) & !is.na(h2) & (h1 + h2) > 0
    h1 <- h1[ok]
    h2 <- h2[ok]
    if (length(h1) == 0) return(rep(NA_real_, 3))
    grp <- ceiling(seq_along(h1) / block)
    counts <- as.numeric(tapply(h1, grp, sum))
    depths <- as.numeric(tapply(h1 + h2, grp, sum))
    big <- depths > maxdepth
    counts[big] <- round(counts[big] / depths[big] * maxdepth)
    depths[big] <- maxdepth
    llh <- vapply(values, function(p)
    {
        a <- dbinom(counts, size=depths, prob=p, log=TRUE)
        b <- dbinom(counts, size=depths, prob=1 - p, log=TRUE)
        m <- pmax(a, b)                       # log-sum-exp: no underflow at high depth
        m[is.infinite(m)] <- 0                # both terms impossible (e.g. p = 1, mixed bin): -Inf, not NaN
        sum(m + log(exp(a - m) + exp(b - m)))
    }, numeric(1))
    llh <- llh - max(llh)
    post <- exp(llh) / sum(exp(llh))
    cs <- cumsum(post)
    c(values[which(cs > .05)[1]], values[which.max(llh)], values[which(cs > .95)[1]])
}

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

## identical to .searchGrid / getNANBfromTot inside getAS_CNA()
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

## BAF-driven allele-specific call of each segment, logR as validation.
##   1. candidate states nMaj+nMin (nMaj >= nMin, 1 <= total <= max_cn); expected folded BAF with purity rho:
##        (rho*nMaj + 1 - rho) / (rho*total + 2*(1 - rho))
##   2. BAF-compatible states: expected BAF inside [q05 - baf_tol, q95 + baf_tol] of the segment's BAF fit
##   3. among those, the state whose total is closest to the logR copy number ntot (ties: closest BAF)
##   4. validation: |ntot - total| <= logr_tol, otherwise no call
## No call (no BAF, no BAF-compatible state, or logR discordant) -> the normal state: 1+1 (CN 2),
## 1+0 on chrX/chrY of a male cell.
## call_status: "called" | "no_BAF" | "no_BAF_state" | "logR_discordant"
callAS_BAF_logR <- function(baf, q05, q95, ntot, chr, purity=1, ismale=FALSE,
                            baf_tol=0.05, logr_tol=0.5, max_cn=6)
{
    states <- do.call("rbind", lapply(1:max_cn, function(tot)
        cbind(nMaj=ceiling(tot/2):tot, nMin=tot-(ceiling(tot/2):tot))))
    tots <- states[,"nMaj"] + states[,"nMin"]
    expBAF <- (purity*states[,"nMaj"] + 1 - purity) / (purity*tots + 2*(1 - purity))

    out <- t(sapply(seq_along(baf), function(i)
    {
        sexchr <- ismale && .hp_chr(chr[i]) %in% c("X", "Y")
        normal <- if(sexchr) c(1, 0) else c(1, 1)
        nocall <- function(status) c(normal, NA, status)
        if(is.na(baf[i]) || is.na(ntot[i])) return(nocall("no_BAF"))
        lo <- if(is.na(q05[i])) baf[i] else q05[i]
        hi <- if(is.na(q95[i])) baf[i] else q95[i]
        ok <- which(expBAF >= lo - baf_tol & expBAF <= hi + baf_tol)
        if(length(ok) == 0) return(nocall("no_BAF_state"))
        best <- ok[order(abs(tots[ok] - ntot[i]), abs(expBAF[ok] - baf[i]))[1]]
        if(abs(tots[best] - ntot[i]) > logr_tol) return(nocall("logR_discordant"))
        c(states[best,"nMaj"], states[best,"nMin"], expBAF[best], "called")
    }))
    data.frame(nMaj_call=as.numeric(out[,1]),
               nMin_call=as.numeric(out[,2]),
               total_call=as.numeric(out[,1]) + as.numeric(out[,2]),
               BAF_expected_call=as.numeric(out[,3]),
               call_status=out[,4],
               stringsAsFactors=FALSE)
}

## Segment BAF used for the call, folded to [0.5, 1] with its interval:
##   "switch"   BAF / q05 / q95 of fitBinom.hap (already folded; robust to phase switches, but biased
##              upward at low depth per block: balanced segments come out > 0.5)
##   "noswitch" BAF_noswitch / q05_noswitch / q95_noswitch (pooled hap1/(hap1+hap2), unbiased when the
##              HP labels are consistent over the segment, i.e. chromosome-scale phasing), folded here;
##              an interval spanning 0.5 folds to [0.5, max(1 - q05, q95)]
.hp_callBAF <- function(nprof, baf_mode=c("switch", "noswitch"))
{
    baf_mode <- match.arg(baf_mode)
    if(baf_mode == "switch")
        return(data.frame(baf=nprof$BAF, lo=nprof$q05, hi=nprof$q95))
    b <- nprof$BAF_noswitch; l <- nprof$q05_noswitch; h <- nprof$q95_noswitch
    l[is.na(l)] <- b[is.na(l)]
    h[is.na(h)] <- b[is.na(h)]
    lo <- ifelse(l <= .5 & h >= .5, .5, ifelse(l > .5, l, 1 - h))
    hi <- ifelse(l <= .5 & h >= .5, pmax(1 - l, h), ifelse(l > .5, h, 1 - l))
    data.frame(baf=pmax(b, 1 - b), lo=lo, hi=hi)
}

## ---------------------------------------------------------------------------------------------
## baf_mode = "model": state likelihoods straight from the per-bin hap1/hap2 counts.
##   per bin j of a segment, state nMaj+nMin with expected major-haplotype fraction
##       f = (purity*(nMaj*(1-eps) + nMin*eps) + (1-purity)) / (purity*total + 2*(1-purity))
##   L_j = 0.5*BB(h1_j | n_j, f, rho) + 0.5*BB(h1_j | n_j, 1-f, rho)
##     - the 0.5/0.5 mixture per bin: HP labels may switch between bins (phase switch errors)
##     - eps: fraction of reads tagged with the wrong haplotype (LOH/1 copy gives f = 1-eps, not 1)
##     - rho: beta-binomial overdispersion (single-cell amplification: reads in a bin are not independent)
##   segment log-likelihood = sum_j log L_j
## call_rule = "logr_first" (default): logR picks the total, the haplotypes pick the split
##   1. eps and rho per cell on a grid, maximising sum over segments of the best state's llh (any total, as
##      before: restricting to logR-compatible states lets logR/BAF-discordant segments inflate eps)
##   2. BAF_fit = the segment's observed major-haplotype fraction under the same model (beta-binomial,
##      switch-tolerant, this cell's rho), i.e. on the same scale as the expected f of a state
##   3. logR-compatible states: |total - ntot| <= logr_tol (none -> "logR_discordant")
##   4. call = logR-compatible state with the highest llh, accepted if |f_state - BAF_fit| <= baf_tol or its
##      llh is within delta of the best state of any total; otherwise "BAF_discordant"
##   The check is on the BAF scale on purpose: a llh tolerance shrinks with depth, so at high depth a
##   0.02 error in eps/purity would reject the right state (e.g. 4+1 at logR CN 5.07, BAF 0.735, eps 0.15).
## call_rule = "baf_first" (previous behaviour): BAF-compatible states = llh >= best llh - delta; among them
##   the total closest to ntot (ties: higher llh); |ntot - total| > logr_tol -> "logR_discordant".
## No call -> 1+1 (1+0 on chrX/Y of a male cell).
## ---------------------------------------------------------------------------------------------
.hp_bbll <- function(k, n, f, rho)
{
    if(rho < 1e-6) return(dbinom(k, n, f, log=TRUE))
    a <- f*(1-rho)/rho
    b <- (1-f)*(1-rho)/rho
    lchoose(n, k) + lbeta(k+a, n-k+b) - lbeta(a, b)
}

.hp_segll <- function(h1, n, f, rho)
{
    a <- .hp_bbll(h1, n, f, rho)
    b <- .hp_bbll(h1, n, 1-f, rho)
    m <- pmax(a, b)
    sum(m + log(.5*exp(a-m) + .5*exp(b-m)))
}

callAS_HPmodel <- function(segbins, use, ntot, chr, purity=1, ismale=FALSE,
                           logr_tol=0.5, max_cn=6, delta=2,
                           baf_tol=0.05, call_rule=c("logr_first", "baf_first"),
                           eps_grid=c(.01,.02,.03,.05,.075,.1,.125,.15,.2,.25),
                           rho_grid=c(0,.01,.02,.05,.1,.2))
{
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
    ## segments x states log-likelihoods (computed once per distinct f)
    llmat <- function(eps, rho)
    {
        f <- fexp(eps)
        uf <- unique(round(f, 8))
        m <- matrix(sapply(idx, function(i)
            sapply(uf, function(ff) .hp_segll(segbins[[i]]$hap1, segbins[[i]]$depth, ff, rho))),
            nrow=length(uf))
        t(m[match(round(f, 8), uf), , drop=FALSE])
    }

    eps <- NA_real_; rho <- NA_real_; LL <- NULL
    if(length(idx) > 0)
    {
        grid <- expand.grid(eps=eps_grid, rho=rho_grid)
        obj <- sapply(seq_len(nrow(grid)), function(g)
            sum(apply(llmat(grid$eps[g], grid$rho[g]), 1, max)))
        best <- which.max(obj)                   # ties -> smallest eps / rho (grid order)
        eps <- grid$eps[best]; rho <- grid$rho[best]
        LL <- llmat(eps, rho)
    }
    fbest <- fexp(if(is.na(eps)) 0 else eps)
    ## observed major-haplotype fraction per segment under the fitted model (same scale as fbest)
    fgrid <- seq(.5, .999, .001)
    baf_fit <- rep(NA_real_, length(ntot))
    for(i in idx)
        baf_fit[i] <- fgrid[which.max(sapply(fgrid, function(ff)
            .hp_segll(segbins[[i]]$hap1, segbins[[i]]$depth, ff, rho)))]

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
    data.frame(nMaj_call=as.numeric(g(1)),
               nMin_call=as.numeric(g(2)),
               total_call=as.numeric(g(1)) + as.numeric(g(2)),
               BAF_expected_call=as.numeric(g(3)),
               BAF_model_best=as.numeric(g(4)),
               BAF_fit=baf_fit,
               call_status=as.character(g(5)),
               tag_error=eps,
               rho=rho,
               stringsAsFactors=FALSE)
}

## One cell: segment BAF from the hap1/hap2 counts of the bins inside each total-CN segment, then the
## purity/ploidy search, nA/nB and ntot_* exactly as getProfile() inside getAS_CNA(), so that
## nprof.fixed (columns BAF, ntot_fixed) feeds getAS_CNA_smoothed() unchanged.
## Extra columns at the end: nBins_hap, hap1, hap2. Per-bin counts in $bin_baf.
getAS_CNA_HP_sample <- function(track,
                                profile,
                                purity,
                                ploidy,
                                purs,
                                ploidies,
                                block=1,
                                min_reads=20,
                                min_bins=2,
                                maxdepth=1000,
                                ismale=FALSE,
                                baf_tol=0.05,
                                logr_tol=0.5,
                                max_cn=6,
                                baf_mode="switch",
                                delta=2,
                                call_rule="logr_first",
                                eps_grid=NULL,
                                rho_grid=NULL,
                                cell_name=NULL)
{
    baf_mode <- match.arg(baf_mode, c("switch", "noswitch", "model"))
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

    nBins_hap <- hap1_seg <- hap2_seg <- rep(0, nrow(nprof))
    segbins <- vector("list", nrow(nprof))
    use <- rep(FALSE, nrow(nprof))
    for(i in unique(qH))
    {
        b <- bins[sort(sH[qH == i]), ]
        b <- b[b$depth > 0, ]
        nBins_hap[i] <- nrow(b)
        hap1_seg[i] <- sum(b$hap1)
        hap2_seg[i] <- sum(b$hap2)
        segbins[[i]] <- b
        if(nrow(b) >= min_bins && sum(b$depth) >= min_reads)
        {
            use[i] <- TRUE
            nprof[i, c("q05","BAF","q95")] <- fitBinom.hap(b$hap1, b$hap2, block=block, maxdepth=maxdepth)
            nprof[i, c("q05_noswitch","BAF_noswitch","q95_noswitch")] <- fitBinom.hap.noswitch(b$hap1, b$hap2)
        }
    }

    nona <- !is.na(nprof[,"BAF"])
    if(sum(nona) == 0) stop("no segment passed min_bins/min_reads - no haplotagged reads?")
    exceptNA <- function(vec) { newvec <- rep(NA, length(nona)); newvec[nona] <- vec; newvec }
    sizes <- nprof[nona,"endpos"] - nprof[nona,"startpos"]
    sG <- .hp_searchGrid(nprof[nona,"BAF"], nprof[nona,"logr"], sizes, purs=purs, ploidies=ploidies)
    sG.fixed <- .hp_searchGrid(nprof[nona,"BAF"], nprof[nona,"logr"], sizes, purs=purity, ploidies=ploidy)
    nprof[,"ntot_free"] <- transform_bulk2tumour(nprof[,"logr"], sG$purity, sG$ploidy)
    nprof[,"ntot_fixed"] <- transform_bulk2tumour(nprof[,"logr"], sG.fixed$purity, sG.fixed$ploidy)

    ## BAF-driven call, logR (ntot of the cell's purity/ploidy fit) as validation
    if(baf_mode == "model")
    {
        calls <- callAS_HPmodel(segbins=segbins, use=use, ntot=nprof[,"ntot_fixed"], chr=nprof[,"chr"],
                                purity=purity, ismale=ismale, logr_tol=logr_tol, max_cn=max_cn, delta=delta,
                                baf_tol=baf_tol, call_rule=call_rule,
                                eps_grid=if(is.null(eps_grid)) eval(formals(callAS_HPmodel)$eps_grid) else eps_grid,
                                rho_grid=if(is.null(rho_grid)) eval(formals(callAS_HPmodel)$rho_grid) else rho_grid)
        calls$call_rule <- call_rule
    }
    else
    {
        cb <- .hp_callBAF(nprof, baf_mode)
        calls <- callAS_BAF_logR(baf=cb$baf, q05=cb$lo, q95=cb$hi,
                                 ntot=nprof[,"ntot_fixed"], chr=nprof[,"chr"],
                                 purity=purity, ismale=ismale,
                                 baf_tol=baf_tol, logr_tol=logr_tol, max_cn=max_cn)
        calls$BAF_call <- cb$baf
    }
    calls$baf_mode <- baf_mode

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
              num.mark=as.numeric(profile[,"num.mark"]),
              calls)
    list(nprof.free=build(sG),
         purity.free=sG$purity,
         ploidy.free=sG$ploidy,
         nprof.fixed=build(sG.fixed),
         purity.fixed=sG.fixed$purity,
         ploidy.fixed=sG.fixed$ploidy,
         bin_baf=bins)
}

## All cells: fills res$allProfiles_AS (same structure as getAS_CNA()), writes one table per cell to
## outdir/as_cna_profile_hp/ and all profiles to outdir/all_as_cna_profiles_hp_<projectname>.pdf.
## Cells that fail get an error string instead of a list (getAS_CNA_smoothed() skips them).
getAS_CNA_HP <- function(res,
                         purs,
                         ploidies,
                         projectname="project",
                         outdir="./",
                         block=1,
                         min_reads=20,
                         min_bins=2,
                         maxdepth=1000,
                         baf_tol=0.05,
                         logr_tol=0.5,
                         max_cn=6,
                         baf_mode="switch",
                         delta=2,
                         call_rule="logr_first",
                         eps_grid=NULL,
                         rho_grid=NULL,
                         mc.cores=1)
{
    cell_names <- names(res$allTracks)
    sex <- rep(if(is.null(res$sex)) "female" else res$sex, length.out=length(cell_names))
    if(is.null(cell_names)) stop("res$allTracks must be named")
    if(!is.list(purs)) purs <- lapply(cell_names, function(x) purs)
    if(!is.list(ploidies)) ploidies <- lapply(cell_names, function(x) ploidies)
    refitted <- any(grepl("refitted", names(res)))

    print(paste0("## derive Allele-specific Profiles from HP-tag haplotype counts (block=", block,
                 ", baf_mode=", baf_mode, if(baf_mode == "model") paste0(", call_rule=", call_rule), ")"))
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
                                     block=block,
                                     min_reads=min_reads,
                                     min_bins=min_bins,
                                     maxdepth=maxdepth,
                                     ismale=sex[x] == "male",
                                     baf_tol=baf_tol,
                                     logr_tol=logr_tol,
                                     max_cn=max_cn,
                                     baf_mode=baf_mode,
                                     delta=delta,
                                     call_rule=call_rule,
                                     eps_grid=eps_grid,
                                     rho_grid=rho_grid,
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
