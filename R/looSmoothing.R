## Leave-one-out (LOO) pseudobulk normalization for single-cell coverage tracks.
##
## For each target cell i we build ONE reference pseudobulk from all *other*
## reference cells and normalize i against it with ASCAT.sc's existing lm-based
## projection (smoothCoverageTrack -> smoothNormals). The reference pool is the
## set of diploid cells (identify_diploid_cells) so clonal signal cannot leak
## into the normal panel. Segmentation and profile fitting are left untouched.

## Sum the per-bin coverage (records) of the given cells into a single
## per-chromosome pseudobulk, reusing combineDiploid on their nlCTS.tumour. The
## returned object keeps the bin layout (space/start/end/width) of nlCTS.tumour
## and only the records column is summed.
make_one_pseudobulk <- function(cell_names, res)
{
    lNormals <- lapply(cell_names, function(cc) res$allTracks[[cc]]$nlCTS.tumour)
    combineDiploid(lNormals)
}

## LOO-normalize every cell in names(res$allTracks) against a diploid pseudobulk.
## The total pool sum S is computed ONCE; the reference for a pool cell i is
## obtained by subtraction (ref_i = S - records_i per chromosome), and cells
## outside the pool are normalized against the full pool sum S. Each reference is
## passed as a length-1 panel (lNormals = list(ref_i)) so smoothCoverageTrack
## routes to the smoothNormals lm projection rather than the single-normal
## subtraction branch. Returns a list of normalized tracks (lCTS shape, with
## fitted/smoothed columns) named by cell.
run_loo_smoothing <- function(res,
                              reference_cells=identify_diploid_cells(res),
                              mc.cores=1)
{
    require(parallel)
    targets <- names(res$allTracks)
    ## Total reference pseudobulk over the whole diploid pool, summed ONCE.
    S <- make_one_pseudobulk(reference_cells, res)
    out <- mclapply(targets, function(i)
    {
        nlCTS.tumour <- res$allTracks[[i]]$nlCTS.tumour
        ## Leave cell i out by subtraction; full pool sum for non-pool cells.
        ref_i <- S
        if(i%in%reference_cells)
            for(j in 1:length(ref_i))
                ref_i[[j]]$records <- S[[j]]$records-nlCTS.tumour[[j]]$records
        lCTS <- smoothCoverageTrack(lCT=nlCTS.tumour,
                                    lSe=res$lSe,
                                    lGCT=res$lGCT,
                                    lNormals=list(ref_i))
        names(lCTS) <- names(nlCTS.tumour)
        lCTS
    }, mc.cores=mc.cores)
    names(out) <- targets
    out
}

## Build a full, plottable ASCAT.sc result whose per-cell tracks are LOO-normalized.
## For each target the reference pseudobulk is S - records_i (full pool sum S for
## non-pool cells); the existing pipeline (getTrackForAll -> searchGrid ->
## fitProfile -> getProfile) is then run unchanged on the normalized track.
## Returns a clone of res with allTracks.processed / allSolutions / allProfiles
## replaced (subset to targets); every other field is kept so ascatsc_plot()
## works on the returned object as-is. targets may be a subset for a quick look,
## but reference_cells should stay the full diploid pool so S is correct.
build_loo_res <- function(res,
                          reference_cells=identify_diploid_cells(res),
                          targets=names(res$allTracks),
                          mc.cores=1)
{
    require(parallel)
    require(DNAcopy)
    allchr <- res$chr
    ## Segmentation boundaries, as in run_sc_sequencing().
    data("SBDRYs_precomputed", package="ASCAT.sc", envir=environment())
    sa <- as.character(res$segmentation_alpha)
    if(!sa%in%names(SBDRYs))
    {
        nperms <- 10000
        max.ones <- floor(nperms*res$segmentation_alpha)+1
        SBDRY <- DNAcopy::getbdry(eta=0.05, nperm=nperms, max.ones=max.ones)
    }
    else SBDRY <- SBDRYs[[sa]]
    ## Total diploid pool sum, computed ONCE.
    S <- make_one_pseudobulk(reference_cells, res)
    cellnames <- names(res$allTracks)
    pos <- match(targets, cellnames)
    getpurs     <- function(k) if(is.list(res$purs))     res$purs[[k]]     else res$purs
    getploidies <- function(k) if(is.list(res$ploidies)) res$ploidies[[k]] else res$ploidies
    processed <- mclapply(targets, function(i)
    {
        nlCTS.tumour <- res$allTracks[[i]]$nlCTS.tumour
        ref_i <- S
        if(i%in%reference_cells)
            for(j in 1:length(ref_i))
                ref_i[[j]]$records <- S[[j]]$records-nlCTS.tumour[[j]]$records
        getTrackForAll(bamfile=NULL,
                       window=NULL,
                       lCT=nlCTS.tumour,
                       lSe=res$lSe,
                       lGCT=res$lGCT,
                       lNormals=list(ref_i),
                       allchr=allchr,
                       sdNormalise=0,
                       SBDRY=SBDRY,
                       segmentation_alpha=res$segmentation_alpha)
    }, mc.cores=mc.cores)
    names(processed) <- targets
    sols <- mclapply(seq_along(targets), function(t)
        try(searchGrid(processed[[t]],
                       purs=getpurs(pos[t]),
                       ploidies=getploidies(pos[t]),
                       maxTumourPhi=res$maxtumourpsi,
                       ismale=isTRUE(res$sex[pos[t]]=="male"),
                       isPON=res$isPON), silent=FALSE), mc.cores=mc.cores)
    names(sols) <- targets
    profs <- mclapply(seq_along(targets), function(t)
        try(getProfile(fitProfile(processed[[t]],
                                  purity=sols[[t]]$purity,
                                  ploidy=sols[[t]]$ploidy,
                                  ismale=isTRUE(res$sex[pos[t]]=="male")),
                       CHRS=allchr), silent=FALSE), mc.cores=mc.cores)
    names(profs) <- targets
    ## Fail loudly rather than return a silently corrupted object: a broken
    ## segmentation/fit (e.g. DNAcopy not on the library path) leaves try-error
    ## placeholders that only surface later as "$ operator is invalid for atomic
    ## vectors" when plotting.
    failed <- targets[sapply(processed, function(x) inherits(x, "try-error") || !is.list(x)) |
                      sapply(profs,     function(x) inherits(x, "try-error") || is.null(x))]
    if(length(failed)>0)
        warning(length(failed), "/", length(targets),
                " cells failed to segment/fit (e.g. ",
                paste(head(failed, 3), collapse=", "),
                "). Are DNAcopy/copynumber available? Run inside the ASCAT.sc container.")
    res_loo <- res
    res_loo$allTracks.processed <- processed
    res_loo$allSolutions <- sols
    res_loo$allProfiles <- profs
    res_loo
}
