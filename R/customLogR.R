## Run the ASCAT.sc fit on a logR that you normalise yourself.
##
## run_sc_sequencing() always corrects log2(records+1) with a loess on GC (smoothCoverageTrack), and
## there is no argument to switch it off. These helpers start from the raw counts kept in
## res$allTracks[[cell]]$nlCTS.tumour, normalise them however you want (nothing, GC loess, tangent/PCA
## on a reference panel, both, or your own code) and put the result back as
## res$allTracks.processed[[cell]] = list(lCTS, lSegs), then segment and fit purity/ploidy the same
## way run_sc_sequencing() does. No BAM is read again.
##
##   cm    <- getCountMatrix(res)                                  ## bins x cells counts + bin table
##   panel <- buildPanelCounts(cm$counts, pool=diploid_cells)      ## bins x n_boot pseudobulk counts
##   logr  <- normaliseCounts(cm$counts, cm$bins, steps="tangent", panel=panel)
##   res2  <- fitFromLogR(res, logr)
##
## Any bins x cells matrix works as logr (rows in the order of cm$bins, colnames = cell names).

getCountMatrix <- function(res, cells=NULL)
{
    if(is.null(res$allTracks)) stop("res$allTracks not found")
    if(is.null(cells)) cells <- names(res$allTracks)
    if(!all(cells%in%names(res$allTracks))) stop("some cells are not in names(res$allTracks)")
    lCT <- res$allTracks[[cells[1]]]$nlCTS.tumour
    nb <- sapply(lCT, nrow)
    chrs <- if(!is.null(res$chr)) res$chr else if(!is.null(names(lCT))) names(lCT) else paste0(seq_along(lCT))
    bins <- data.frame(chr=rep(chrs, nb),
                       chr_index=rep(seq_along(lCT), nb),
                       start=unlist(lapply(lCT, function(x) x$start)),
                       end=unlist(lapply(lCT, function(x) x$end)),
                       stringsAsFactors=FALSE)
    gc <- unlist(res$lGCT)
    if(length(gc)==nrow(bins))
        bins$gc <- as.numeric(gc)
    else
        warning("res$lGCT does not match the bins of res$allTracks - no gc column")
    counts <- vapply(cells, function(cc)
        as.numeric(unlist(lapply(res$allTracks[[cc]]$nlCTS.tumour, function(x) x$records))),
        numeric(nrow(bins)))
    rownames(counts) <- paste0(bins$chr, ":", bins$start, "-", bins$end)
    list(counts=counts, bins=bins)
}

## logR currently stored in res$allTracks.processed (e.g. the GC-corrected one from run_sc_sequencing)
getProcessedLogR <- function(res, cells=NULL, column="smoothed")
{
    if(is.null(cells)) cells <- names(res$allTracks.processed)
    nbins <- sum(sapply(res$allTracks.processed[[cells[1]]]$lCTS, nrow))
    vapply(cells, function(cc)
        as.numeric(unlist(lapply(res$allTracks.processed[[cc]]$lCTS, function(x) x[[column]]))),
        numeric(nbins))
}

## Reference panel for the tangent step: n_boot pseudobulks, each the sum of cells_per_bulk cells
## sampled from pool. Use cells you trust to be diploid as pool: CNAs shared by the pool end up in
## the panel and are removed from every cell by the projection.
buildPanelCounts <- function(counts,
                             pool=NULL,
                             n_boot=200,
                             cells_per_bulk=100,
                             seed=1)
{
    cells <- if(is.null(pool)) colnames(counts) else intersect(pool, colnames(counts))
    if(length(cells)<2) stop("need at least 2 cells of the pool in colnames(counts)")
    k <- min(cells_per_bulk, length(cells))
    if(k<cells_per_bulk) warning(paste0("cells_per_bulk > pool size, using ", k))
    if(!is.null(seed)) set.seed(seed)
    panel <- vapply(seq_len(n_boot), function(b) rowSums(counts[, sample(cells, k), drop=FALSE]),
                    numeric(nrow(counts)))
    rownames(panel) <- rownames(counts)
    colnames(panel) <- paste0("pseudobulk_", seq_len(n_boot))
    panel
}

## steps are applied in the given order to log2(counts+1) (median-centred per cell):
##   "none"    nothing (raw log coverage)
##   "gc"      per-cell loess on the bin covariates (default "gc", as smoothCoverageTrack); add any other
##             bin-level column of bins to covariates, e.g. a pseudobulk read-length track
##   "tangent" remove the projection on the first k left singular vectors of the panel; k=NULL takes
##             the smallest k reaching var_explained. The panel goes through the steps before "tangent"
##             too, so that cells and reference are normalised the same way.
normaliseCounts <- function(counts,
                            bins,
                            steps="gc",
                            panel=NULL,
                            k=NULL,
                            var_explained=0.9,
                            center=FALSE,
                            covariates="gc",
                            MC.CORES=1)
{
    steps <- match.arg(steps, c("none","gc","tangent"), several.ok=TRUE)
    steps <- steps[steps!="none"]
    if("tangent"%in%steps & is.null(panel))
        stop("steps include 'tangent' but no panel is given (see buildPanelCounts)")
    if(!is.null(panel) && nrow(panel)!=nrow(counts))
        stop("panel and counts must have the same bins (rows)")
    if("gc"%in%steps & !all(covariates%in%colnames(bins)))
        stop(paste0("covariates not found in bins: ", paste(setdiff(covariates, colnames(bins)), collapse=", ")))
    .logc <- function(m)
    {
        m <- log2(m+1)
        sweep(m, 2, apply(m, 2, median), "-")
    }
    logr <- .logc(counts)
    ref <- if(is.null(panel)) NULL else .logc(panel)
    info <- list(steps=if(length(steps)==0) "none" else steps)
    for(i in seq_along(steps))
    {
        if(steps[i]=="gc")
        {
            print(paste0("## loess on ", paste(covariates, collapse=" + ")))
            logr <- .loessCorrect(logr, bins[, covariates, drop=FALSE], MC.CORES=MC.CORES)
            if(!is.null(ref) & "tangent"%in%steps[-seq_len(i)])
                ref <- .loessCorrect(ref, bins[, covariates, drop=FALSE], MC.CORES=MC.CORES)
            info$covariates <- covariates
        }
        if(steps[i]=="tangent")
        {
            tg <- .tangentCorrect(logr, ref, k=k, var_explained=var_explained, center=center)
            print(paste0("## tangent: k=", tg$k, " components (", round(100*tg$var_explained, 1),
                         "% of panel variance)"))
            logr <- tg$logr
            info$k <- tg$k
            info$center <- center
            info$n_panel <- ncol(panel)
        }
    }
    attr(logr, "normalisation") <- info
    ## logR relative to a reference panel: male X is then handled as with a PON (see transform_bulk2tumour)
    attr(logr, "isPON") <- "tangent"%in%steps
    logr
}

.loessCorrect <- function(m, covars, MC.CORES=1)
{
    covars <- as.data.frame(lapply(covars, function(v)
    {
        v <- as.numeric(v)
        v[!is.finite(v)] <- median(v[is.finite(v)])
        v
    }))
    form <- stats::as.formula(paste("y ~", paste(colnames(covars), collapse=" + ")))
    out <- parallel::mclapply(seq_len(ncol(m)), function(j)
    {
        df <- cbind(y=m[,j], covars)
        df$y-stats::fitted(myloess(LL=nrow(df), form, data=df))
    }, mc.cores=MC.CORES)
    bad <- sapply(out, inherits, "try-error")
    if(any(bad)) stop(paste0("loess failed for: ", paste(colnames(m)[bad], collapse=", ")))
    out <- do.call(cbind, out)
    dimnames(out) <- dimnames(m)
    out
}

.tangentCorrect <- function(logr, ref, k=NULL, var_explained=0.9, center=FALSE)
{
    mu <- if(center) rowMeans(ref) else rep(0, nrow(ref))
    sv <- svd(ref-mu, nv=0)
    ve <- cumsum(sv$d^2)/sum(sv$d^2)
    if(is.null(k)) k <- which(ve>=var_explained)[1]
    k <- max(1L, min(as.integer(k), length(sv$d)))
    U <- sv$u[, seq_len(k), drop=FALSE]
    y <- logr-mu
    list(logr=y-U%*%crossprod(U, y), k=k, var_explained=ve[k])
}

## Per-cell summaries of a logR matrix, to compare normalisations before fitting:
## mapd (median absolute difference of neighbouring bins, within chromosomes), sd, and the Spearman
## correlation with GC left after normalisation
getLogRStats <- function(logr, bins)
{
    same <- bins$chr[-1]==bins$chr[-nrow(bins)]
    st <- data.frame(cell=colnames(logr),
                     mapd=apply(logr, 2, function(x) median(abs(diff(x)[same]))),
                     sd=apply(logr, 2, sd),
                     stringsAsFactors=FALSE)
    if(!is.null(bins$gc))
        st$cor_gc <- apply(logr, 2, function(x) suppressWarnings(cor(x, bins$gc, method="spearman")))
    rownames(st) <- NULL
    st
}

## Segment logr (bins x cells, rows in the order of getCountMatrix()$bins) and fit purity/ploidy as
## run_sc_sequencing() does, for all columns of logr or only cells. Returns res restricted to those
## cells, with new allTracks.processed, allSolutions and allProfiles; results derived from the previous
## logR (refits, allele-specific profiles, QC, filters) are dropped. purs/ploidies/maxtumourpsi/
## segmentation_alpha/segmentation default to the values stored in res, and segmentation_penalty too when
## segmentation is the one of res (see run_sc_sequencing). isPON defaults to TRUE after a tangent step
## of normaliseCounts(); that tag is lost when logr is subset (logr[, x]), so use cells= or set isPON.
fitFromLogR <- function(res,
                        logr,
                        cells=NULL,
                        purs=NULL,
                        ploidies=NULL,
                        maxtumourpsi=NULL,
                        segmentation_alpha=NULL,
                        segmentation=NULL,
                        segmentation_penalty=NULL,
                        isPON=NULL,
                        predict_refit=FALSE,
                        MC.CORES=1)
{
    if(is.null(colnames(logr)) || !all(colnames(logr)%in%names(res$allTracks)))
        stop("colnames(logr) must be cell names from names(res$allTracks)")
    if(is.null(cells)) cells <- colnames(logr)
    if(!all(cells%in%colnames(logr))) stop("some cells are not in colnames(logr)")
    if(any(!is.finite(logr[, cells]))) stop("logr has NA/Inf values")
    idx <- match(cells, names(res$allTracks))
    jdx <- match(cells, colnames(logr))
    nb <- sapply(res$allTracks[[cells[1]]]$nlCTS.tumour, nrow)
    if(sum(nb)!=nrow(logr)) stop("nrow(logr) does not match the bins of res$allTracks (see getCountMatrix)")
    ends <- cumsum(nb)
    starts <- ends-nb+1
    allchr <- if(!is.null(res$chr)) res$chr else names(res$allTracks[[cells[1]]]$nlCTS.tumour)
    pick <- function(x, default) if(is.null(x)) default else x
    segmentation_alpha <- pick(segmentation_alpha, pick(res$segmentation_alpha, 0.01))
    ## [[ ]]: res$segmentation would partially match res$segmentation_alpha in results without segmentation
    segmentation <- match.arg(pick(segmentation, pick(res[["segmentation"]], "cbs")), c("cbs","pcf","quickseg"))
    if(is.null(segmentation_penalty) && identical(segmentation, res[["segmentation"]]))
        segmentation_penalty <- res[["segmentation_penalty"]]
    maxtumourpsi <- pick(maxtumourpsi, pick(res$maxtumourpsi, 5))
    if(is.null(isPON))
    {
        isPON <- attr(logr, "isPON")
        if(is.null(isPON))
        {
            isPON <- isTRUE(res$isPON)
            print(paste0("## logr has no normalisation tag - using isPON=", isPON,
                         " from res (set isPON=TRUE if logr is relative to a reference panel)"))
        }
    }
    sex <- pick(res$sex, "female")
    if(length(sex)!=length(res$allTracks)) sex <- rep(sex[1], length(res$allTracks))
    sex <- sex[idx]
    if(is.null(purs)) purs <- if(is.list(res$purs)) res$purs[idx] else pick(res$purs, seq(0.1, 1, 0.01))
    if(is.null(ploidies)) ploidies <- if(is.list(res$ploidies)) res$ploidies[idx] else pick(res$ploidies, seq(1.7, 5, 0.01))
    if(!is.list(purs)) purs <- lapply(cells, function(x) purs)
    if(!is.list(ploidies)) ploidies <- lapply(cells, function(x) ploidies)
    SBDRYs <- NULL
    data("SBDRYs_precomputed", package="ASCAT.sc", envir=environment())
    if(as.character(segmentation_alpha)%in%names(SBDRYs))
        SBDRY <- SBDRYs[[as.character(segmentation_alpha)]]
    else
        SBDRY <- DNAcopy::getbdry(eta=0.05, nperm=10000, max.ones=floor(10000*segmentation_alpha)+1)
    .fitOne <- function(j)
    {
        lCTS <- lapply(seq_along(nb), function(i)
        {
            x <- res$allTracks[[cells[j]]]$nlCTS.tumour[[i]]
            x <- x[, !colnames(x)%in%c("fitted","smoothed"), drop=FALSE]
            sm <- as.numeric(logr[starts[i]:ends[i], jdx[j]])
            cbind(x, fitted=log2(x$records+1)-sm, smoothed=sm)
        })
        names(lCTS) <- allchr
        tseg <- system.time(
        {
            sigma <- if(segmentation%in%c("pcf","quickseg")) getTrackNoise(lapply(lCTS, function(x) x$smoothed)) else NULL
            lSegs <- lapply(seq_along(lCTS), function(i)
                segmentTrack(lCTS[[i]]$smoothed,
                             chr=paste0(i),
                             starts=lCTS[[i]]$start,
                             ends=lCTS[[i]]$end,
                             SBDRY=SBDRY,
                             ALPHA=segmentation_alpha,
                             segmentation=segmentation,
                             penalty=segmentation_penalty,
                             sigma=sigma))
        })[["elapsed"]]
        names(lSegs) <- paste0(seq_along(lCTS))
        track <- list(lCTS=lCTS, lSegs=lSegs)
        ismale <- sex[j]=="male"
        sol <- try(searchGrid(track,
                              purs=purs[[j]],
                              ploidies=ploidies[[j]],
                              maxTumourPhi=maxtumourpsi,
                              ismale=ismale,
                              isPON=isPON), silent=TRUE)
        ## isPON is also passed to fitProfile (run_sc_sequencing does not), so that male X in the
        ## profile matches the grid search
        prof <- try(getProfile(fitProfile(track,
                                          purity=sol$purity,
                                          ploidy=sol$ploidy,
                                          ismale=ismale,
                                          isPON=isPON),
                               CHRS=allchr), silent=TRUE)
        list(track=track, sol=sol, prof=prof, time=tseg)
    }
    print(paste0("## segment (", segmentation, ") and fit ", length(cells), " cells (isPON=", isPON, ")"))
    out <- parallel::mclapply(seq_along(cells), function(j) try(.fitOne(j), silent=TRUE), mc.cores=MC.CORES)
    out <- lapply(out, function(o) if(inherits(o, "try-error")) list(track=o, sol=o, prof=o) else o)
    nfail <- sum(sapply(out, function(o) inherits(o$sol, "try-error")))
    if(nfail>0) print(paste0("## ", nfail, " cells failed (try-error in allSolutions)"))
    stale <- grepl("^allSolutions|^allProfiles|_AS|refitted|^summary$|^QC|^pass$|^fail$|^cells$|^filters|^timeto(processed|fit|segment)$",
                   names(res))
    res <- res[!stale]
    res$allTracks <- res$allTracks[cells]
    res$allTracks.processed <- lapply(out, function(o) o$track)
    res$allSolutions <- lapply(out, function(o) o$sol)
    res$allProfiles <- lapply(out, function(o) o$prof)
    names(res$allTracks.processed) <- names(res$allSolutions) <- names(res$allProfiles) <- cells
    ## segmentation time per cell (s), NA for cells that failed
    res$timetosegment <- setNames(sapply(out, function(o) if(is.null(o$time)) NA else o$time), cells)
    res$sex <- sex
    res$purs <- purs
    res$ploidies <- ploidies
    res$maxtumourpsi <- maxtumourpsi
    res$segmentation_alpha <- segmentation_alpha
    res$segmentation <- segmentation
    res$segmentation_penalty <- segmentation_penalty
    res$isPON <- isPON
    res$multipcf <- FALSE
    res$normalisation <- pick(attr(logr, "normalisation"), list(steps="custom"))
    if(predict_refit)
    {
        print("## predict Refit")
        res <- predictRefit_all(res)
    }
    res
}
