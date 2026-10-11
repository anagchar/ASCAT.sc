## Bin-level read-length covariate pooled across cells ("pseudobulk"), aligned with the bins of the
## input tracks (nlCTS.tumour after treatTrack / sc_excludeBadBins, i.e. the same bins as res$nlGCT).
## Used by smoothCoverageTrack(readlen="pseudobulk") in place of the per-cell mean read length, which
## is far too noisy at single-cell depth (few, heavy-tailed read lengths per bin).
##
## lTracks:           named list over cells; each element is the per-chromosome list of coverage
##                    data.frames of one cell (res$allTracks[[cell]]$nlCTS.tumour)
## stat:              "short_fraction": pooled number of reads shorter than short_threshold / pooled
##                                      records (needs the n_len_* band columns written by
##                                      getCoverageTrack.10XBAM; falls back to "mean_length" otherwise)
##                    "mean_length":    pooled nucleotides / pooled records
## short_threshold:   band boundary (bp) defining "short"; NULL = the boundary closest to the pooled
##                    median read length (the fraction then behaves like a sign test)
## groups:            optional named vector cell -> group (e.g. res$clustering$clusters); one track per
##                    group is returned (a list), a group falling back to the global track when its
##                    pooled depth is below 5*min_reads_per_bin; attr(, "global") holds the global track
## min_reads_per_bin: bins with fewer pooled reads get the genome-wide median value (imputed)
##
## Returns a list per chromosome of numeric vectors (same shape as lGCT) with attributes "stat",
## "short_threshold", "n_imputed", "median_depth", "n_cells"; with groups, a list of such tracks.
getReadLengthTrack <- function(lTracks,
                               stat=c("short_fraction","mean_length"),
                               short_threshold=NULL,
                               groups=NULL,
                               min_reads_per_bin=10,
                               verbose=TRUE)
{
    stat <- match.arg(stat)
    if(length(lTracks)==0) stop("getReadLengthTrack: no tracks provided")
    if(is.null(names(lTracks))) names(lTracks) <- paste0("cell", seq_along(lTracks))
    cells <- names(lTracks)
    nchr <- length(lTracks[[1]])
    chrnames <- names(lTracks[[1]])
    band_cols <- grep("^n_len_", colnames(lTracks[[1]][[1]]), value=TRUE)
    if(stat=="short_fraction" & length(band_cols)==0)
    {
        if(verbose) print("getReadLengthTrack: no read-length band columns (n_len_*) in the tracks - falling back to stat='mean_length' (re-read the BAMs to get the bands)")
        stat <- "mean_length"
    }
    poolcol <- function(cc, cells)
    {
        lapply(seq_len(nchr), function(chr)
            Reduce(`+`, lapply(cells, function(cl) as.numeric(lTracks[[cl]][[chr]][[cc]]))))
    }
    short_cols <- NULL
    if(stat=="short_fraction")
    {
        upper <- as.numeric(sub("^n_len_[0-9.]+_", "", band_cols))
        ord <- order(upper)
        band_cols <- band_cols[ord]
        upper <- upper[ord]
        bounds <- upper[is.finite(upper)]
        if(is.null(short_threshold))
        {
            tot <- sapply(band_cols, function(cc) sum(unlist(poolcol(cc, cells)), na.rm=TRUE))
            cumfrac <- cumsum(tot)/sum(tot)
            short_threshold <- bounds[which.min(abs(cumfrac[seq_along(bounds)]-0.5))]
            if(verbose) print(paste0("getReadLengthTrack: short_threshold set to ", short_threshold,
                                     " bp (band boundary closest to the pooled median; ",
                                     round(100*cumfrac[match(short_threshold, bounds)]),
                                     "% of reads are shorter)"))
        }
        if(!short_threshold %in% bounds)
            stop(paste0("getReadLengthTrack: short_threshold must be one of the band boundaries: ",
                        paste(bounds, collapse=", ")))
        short_cols <- band_cols[upper <= short_threshold]
    }
    pool <- function(cells)
    {
        rec <- poolcol("records", cells)
        if(stat=="short_fraction")
        {
            s <- lapply(short_cols, function(cc) poolcol(cc, cells))
            num <- lapply(seq_len(nchr), function(chr) Reduce(`+`, lapply(s, `[[`, chr)))
        }
        else
            num <- poolcol("nucleotides", cells)
        out <- lapply(seq_len(nchr), function(chr)
            ifelse(rec[[chr]]>=min_reads_per_bin, num[[chr]]/rec[[chr]], NA))
        allx <- unlist(out)
        med <- median(allx, na.rm=TRUE)
        out <- lapply(out, function(x) { x[!is.finite(x)] <- med; x })
        names(out) <- chrnames
        attr(out, "stat") <- stat
        attr(out, "short_threshold") <- short_threshold
        attr(out, "n_imputed") <- sum(!is.finite(allx))
        attr(out, "median_depth") <- median(unlist(rec))
        attr(out, "n_cells") <- length(cells)
        out
    }
    glob <- pool(cells)
    if(verbose) print(paste0("getReadLengthTrack: ", stat, " pooled over ", length(cells),
                             " cells; median pooled reads per bin = ", attr(glob, "median_depth"),
                             "; imputed bins = ", attr(glob, "n_imputed")))
    if(is.null(groups)) return(glob)
    grp <- as.character(groups[cells])
    grp[is.na(grp)] <- "ungrouped"
    out <- lapply(split(cells, grp), function(cg)
    {
        depth <- median(unlist(poolcol("records", cg)))
        if(depth < 5*min_reads_per_bin)
        {
            if(verbose) print(paste0("getReadLengthTrack: group of ", length(cg),
                                     " cells too sparse (median pooled reads per bin = ", depth,
                                     ") - using the global track for it"))
            return(glob)
        }
        pool(cg)
    })
    if(verbose) print(paste0("getReadLengthTrack: ", length(out), " group track(s): ",
                             paste(names(out), " (", sapply(split(cells, grp), length), " cells)", sep="", collapse=", ")))
    attr(out, "global") <- glob
    attr(out, "groups") <- setNames(grp, cells)
    attr(out, "stat") <- stat
    attr(out, "short_threshold") <- short_threshold
    out
}

## Pick the read-length track of one cell from the output of getReadLengthTrack():
## the global track, or the cell's group track when groups were used (global if the cell has no group).
pickReadLengthTrack <- function(nlRLT, cell)
{
    if(is.null(nlRLT)) return(NULL)
    grp <- attr(nlRLT, "groups")
    if(is.null(grp)) return(nlRLT)
    g <- if(is.null(cell) || is.na(cell)) NA else grp[cell]
    if(is.na(g) || is.null(nlRLT[[g]])) return(attr(nlRLT, "global"))
    nlRLT[[g]]
}
