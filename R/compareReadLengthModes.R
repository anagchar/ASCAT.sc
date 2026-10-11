## Compare read-length correction modes on an existing single-cell result object (10X mode), using
## the end-points of the pipeline: for every cell and mode the corrected logR is segmented, purity and
## ploidy are fitted and the copy-number profile is derived exactly as in run_sc_sequencing(), without
## a normal reference. Meant to be run on cells expected to be diploid (normal cells): a better
## correction gives lower MAPD, fewer breakpoints and less genome called non-diploid.
##
## res:         result of run_sc_sequencing() (uses res$allTracks[[cell]]$nlCTS.tumour, res$lGCT, res$lSe,
##              res$chr, res$sex, res$purs/res$ploidies/res$maxtumourpsi and, if by_cluster=TRUE,
##              res$clustering$clusters)
## cells:       cells to evaluate, e.g. the normal barcodes (NULL = n_cells random cells)
## modes:       any of "none" (GC only), "cell" (per-cell mean read length), "pseudobulk"
## stats:       pseudobulk statistics to try ("short_fraction" needs the n_len_* band columns, i.e. BAMs
##              read with this version; otherwise it falls back to "mean_length")
## purs, ploidies: purity/ploidy grids for the fit (NULL = the grids stored in res, as in the pipeline);
##              for normal cells a restricted grid (e.g. purs=1, ploidies=seq(1.8,2.2,0.01)) is faster
## expected_cn: copy number expected over the genome (2; chrX/Y of male cells expect 1)
##
## Returns one row per cell x mode:
##   mapd:          median |diff| of the corrected logR over adjacent bins (noise)
##   n_breakpoints: number of CBS breakpoints over all chromosomes (segments - chromosomes)
##   frac_altered:  fraction of the genome (bins) with a called total copy number != expected_cn
##   purity, ploidy, ambiguous: the fitted solution
##   sd_logr, leak_slope, leak_r2: secondary; leak = linear dependence of the corrected logR on the
##                  pooled read-length covariate (0 = read-length effect removed)
## A per-mode summary (medians over cells, and the proportion of cells with <= 3% of the genome
## altered) is printed and attached as attr(, "summary").
compareReadLengthModes <- function(res,
                                   cells=NULL,
                                   n_cells=100,
                                   modes=c("none","cell","pseudobulk"),
                                   stats=c("short_fraction","mean_length"),
                                   short_threshold=NULL,
                                   by_cluster=FALSE,
                                   min_reads_per_bin=10,
                                   purs=NULL,
                                   ploidies=NULL,
                                   expected_cn=2,
                                   segmentation_alpha=NULL,
                                   mc.cores=1,
                                   seed=1)
{
    suppressPackageStartupMessages(require(parallel))
    suppressPackageStartupMessages(require(DNAcopy))
    if(is.null(res$allTracks)) stop("compareReadLengthModes: res$allTracks not found")
    lGCT <- if(!is.null(res$nlGCT)) res$nlGCT else res$lGCT
    lSe  <- if(!is.null(res$nlSe)) res$nlSe else res$lSe
    allchr <- if(!is.null(res$chr)) res$chr else names(lSe)
    allcells <- names(res$allTracks)
    if(is.null(allcells)) allcells <- paste0("cell", seq_along(res$allTracks))
    lTracks <- lapply(res$allTracks, function(x) x$nlCTS.tumour)
    names(lTracks) <- allcells
    if(sum(sapply(lGCT, length))!=sum(sapply(lTracks[[1]], nrow)))
        stop("compareReadLengthModes: GC track and coverage tracks have different numbers of bins")
    set.seed(seed)
    if(is.null(cells)) cells <- sample(allcells, min(n_cells, length(allcells)))
    cells <- cells[cells%in%allcells]
    if(length(cells)==0) stop("compareReadLengthModes: none of the requested cells are in res$allTracks")
    ## pseudobulk tracks, pooled over ALL cells of res
    groups <- if(by_cluster) res$clustering$clusters else NULL
    tracks <- list()
    if("pseudobulk"%in%modes)
        for(st in stats)
        {
            tr <- try(getReadLengthTrack(lTracks, stat=st,
                                         short_threshold=if(st=="short_fraction") short_threshold else NULL,
                                         groups=groups, min_reads_per_bin=min_reads_per_bin), silent=TRUE)
            if(inherits(tr, "try-error"))
            {
                print(paste0("compareReadLengthModes: skipping stat ", st, ": ", conditionMessage(attr(tr, "condition"))))
                next
            }
            tracks[[attr(tr, "stat")]] <- tr   # the stat actually used (short_fraction may fall back to mean_length)
        }
    ## reference covariate for the leak metric: the global short-fraction track if available, else mean length
    refname <- if("short_fraction"%in%names(tracks)) "short_fraction" else if(length(tracks)>0) names(tracks)[1] else NULL
    if(is.null(refname))
        refcov <- getReadLengthTrack(lTracks, stat="mean_length", min_reads_per_bin=min_reads_per_bin, verbose=FALSE)
    else
    {
        refcov <- tracks[[refname]]
        if(!is.null(attr(refcov, "global"))) refcov <- attr(refcov, "global")
    }
    refcov <- unlist(refcov)
    ## segmentation and fitting settings, as in run_sc_sequencing()
    alpha <- if(!is.null(segmentation_alpha)) segmentation_alpha else if(!is.null(res$segmentation_alpha)) res$segmentation_alpha else 0.01
    data("SBDRYs_precomputed", package="ASCAT.sc", envir=environment())
    if(as.character(alpha)%in%names(SBDRYs))
        SBDRY <- SBDRYs[[as.character(alpha)]]
    else
    {
        nperms <- 10000
        SBDRY <- DNAcopy::getbdry(eta=0.05, nperm=nperms, max.ones=floor(nperms*alpha)+1)
    }
    maxphi <- if(!is.null(res$maxtumourpsi)) res$maxtumourpsi else 5
    ## purity/ploidy grid: user-given, else the one stored in res (a list per cell after the pipeline), else default
    getgrid <- function(user, stored, cellidx, default)
    {
        if(!is.null(user)) return(user)
        if(is.list(stored)) return(stored[[cellidx]])
        if(!is.null(stored)) return(stored)
        default
    }
    sexv <- res$sex
    if(is.null(sexv) || length(sexv)!=length(allcells)) sexv <- rep(if(is.null(sexv)) "female" else sexv[1], length(allcells))
    names(sexv) <- allcells
    runs <- list()
    for(m in modes)
    {
        if(m=="pseudobulk")
            for(st in names(tracks)) runs[[paste0("pseudobulk_", st)]] <- list(readlen="pseudobulk", track=tracks[[st]])
        else
            runs[[m]] <- list(readlen=m, track=NULL)
    }
    if(length(runs)==0) stop("compareReadLengthModes: nothing to run (check modes/stats)")
    isXY <- function(chr) gsub("chr","",as.character(chr))%in%c("X","Y","23","24")
    evalcell <- function(cell)
    {
        cellidx <- match(cell, allcells)
        lCT <- lTracks[[cell]]
        nreads <- sum(unlist(lapply(lCT, function(x) x$records)))
        ismale <- isTRUE(sexv[cell]=="male")
        purs_cell <- getgrid(purs, res$purs, cellidx, seq(0.1, 1, 0.01))
        ploidies_cell <- getgrid(ploidies, res$ploidies, cellidx, seq(1.7, 5, 0.01))
        do.call(rbind, lapply(names(runs), function(rn)
        {
            r <- runs[[rn]]
            invisible(capture.output({
                tr <- getTrackForAll(bamfile=NULL, window=NULL, lCT=lCT, lSe=lSe, lGCT=lGCT, lNormals=NULL,
                                     allchr=allchr, sdNormalise=0, SBDRY=SBDRY, segmentation_alpha=alpha,
                                     lRLT=if(is.null(r$track)) NULL else pickReadLengthTrack(r$track, cell),
                                     readlen=r$readlen)
                sol <- try(searchGrid(tr, purs=purs_cell, ploidies=ploidies_cell, maxTumourPhi=maxphi,
                                      ismale=ismale, isPON=FALSE), silent=TRUE)
                prof <- if(inherits(sol, "try-error")) sol else
                    try(getProfile(fitProfile(tr, purity=sol$purity, ploidy=sol$ploidy, ismale=ismale), CHRS=allchr), silent=TRUE)
            }))
            sm <- unlist(lapply(tr$lCTS, function(x) x$smoothed))
            fit <- lm(sm ~ refcov)
            nbp <- sum(sapply(tr$lSegs, function(s) nrow(s$output)-1))
            fitok <- !inherits(sol, "try-error") && !inherits(prof, "try-error") && is.data.frame(prof)
            frac_altered <- NA_real_
            if(fitok)
            {
                cn <- prof$total_copy_number
                w <- prof$num.mark
                expected <- ifelse(ismale & isXY(prof$chromosome), 1, expected_cn)
                ok <- !is.na(cn) & !is.na(w)
                frac_altered <- sum(w[ok & cn!=expected])/sum(w[ok])
            }
            data.frame(cell=cell,
                       mode=rn,
                       n_reads=nreads,
                       cluster=if(is.null(res$clustering$clusters)) NA else unname(res$clustering$clusters[cell]),
                       purity=if(fitok) sol$purity else NA_real_,
                       ploidy=if(fitok) sol$ploidy else NA_real_,
                       ambiguous=if(fitok) isTRUE(sol$ambiguous) else NA,
                       mapd=median(abs(diff(sm))),
                       n_breakpoints=nbp,
                       frac_altered=frac_altered,
                       sd_logr=sd(sm),
                       leak_slope=unname(coef(fit)[2]),
                       leak_r2=summary(fit)$r.squared,
                       stringsAsFactors=FALSE)
        }))
    }
    print(paste0("compareReadLengthModes: ", length(cells), " cells x ", length(runs), " modes (", paste(names(runs), collapse=", "), ")"))
    out <- do.call(rbind, mclapply(cells, evalcell, mc.cores=mc.cores))
    rownames(out) <- NULL
    med <- function(x) median(x, na.rm=TRUE)
    summ <- do.call(rbind, lapply(names(runs), function(rn)
    {
        o <- out[out$mode==rn,]
        data.frame(mode=rn,
                   mapd=med(o$mapd),
                   n_breakpoints=med(o$n_breakpoints),
                   frac_altered=med(o$frac_altered),
                   prop_cells_le3pct_altered=mean(o$frac_altered<=0.03, na.rm=TRUE),
                   ploidy=med(o$ploidy),
                   purity=med(o$purity),
                   n_failed_fits=sum(is.na(o$frac_altered)),
                   leak_slope=med(o$leak_slope),
                   stringsAsFactors=FALSE)
    }))
    print("Median over cells per mode (lower mapd / n_breakpoints / frac_altered = better on diploid cells; leak_slope ~ 0 = read-length effect removed):")
    print(summ, row.names=FALSE, digits=3)
    attr(out, "summary") <- summ
    attr(out, "tracks") <- tracks
    out
}
