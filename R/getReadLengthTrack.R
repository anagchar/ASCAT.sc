## Bin-level read-length track pooled across cells ("pseudobulk"), aligned with the bins of the input
## tracks (nlCTS.tumour after treatTrack / sc_excludeBadBins, i.e. the same bins as res$nlGCT).
##
## Why: Tn5 cuts more in open chromatin, so an accessible bin gets more reads, and shorter ones
## (nucleosome-free fragments), in every cell of that cell type. The extra reads look like a copy-number
## gain. The read length of a bin marks the accessibility without carrying the copy number: a gain
## multiplies the short and the long reads of a bin alike, so their median length does not move. One
## cell has too few reads per bin for its own median to be informative, so the median is pooled over
## the cells, and smoothCoverageTrack(readlen="pseudobulk") uses it next to GC in each cell's loess.
## Pooled read COUNTS are not used as the covariate: they also contain the clonal copy-number changes
## shared by the tumour cells, which the loess would then remove from every cell.
##
## lTracks:           named list over cells; each element is the per-chromosome list of coverage
##                    data.frames of one cell (res$allTracks[[cell]]$nlCTS.tumour)
## stat:              "median": records-weighted median over cells of the per-cell bin medians
##                               (median_qwidth column); falls back to "mean" when the column is missing
##                    "mean":   pooled nucleotides / pooled records
## groups:            optional named vector cell -> group (e.g. res$clustering$clusters); one track per
##                    group is returned, a group falling back to the global track when its median pooled
##                    reads per bin is below 5*min_reads_per_bin; attr(, "global") holds the global track
## min_reads_per_bin: bins with fewer pooled reads get the genome-wide median value of the track
##
## Returns a list per chromosome of numeric vectors (same shape as lGCT) with attributes "stat",
## "n_imputed", "median_depth", "n_cells"; with groups, a list of such tracks (see pickReadLengthTrack).
getReadLengthTrack <- function(lTracks,
                               stat=c("median","mean"),
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
    if(stat=="median" & !all(sapply(lTracks, function(x) "median_qwidth" %in% colnames(x[[1]]))))
    {
        if(verbose) print("getReadLengthTrack: no median_qwidth column in the tracks - using stat='mean' (pooled nucleotides/records)")
        stat <- "mean"
    }
    ## bins x cells matrix of one column of the tracks, for one chromosome
    colmat <- function(cc, chr, cg)
        matrix(sapply(cg, function(cl) as.numeric(lTracks[[cl]][[chr]][[cc]])), ncol=length(cg))
    ## records-weighted median of each row (NA values and zero weights dropped);
    ## equals median(rep(x, w)) for integer weights w, i.e. the median of the pooled per-cell values
    wmedianRows <- function(X, Wt)
    {
        out <- rep(NA_real_, nrow(X))
        ok <- is.finite(X) & Wt > 0
        if(!any(ok)) return(out)
        r <- row(X)[ok]; x <- X[ok]; w <- Wt[ok]
        o <- order(r, x); r <- r[o]; x <- x[o]; w <- w[o]
        cw <- ave(w, r, FUN=cumsum)
        n <- ave(w, r, FUN=sum)
        k1 <- floor((n+1)/2); k2 <- floor(n/2)+1   # the two middle positions (equal when n is odd)
        rows <- unique(r)
        first <- match(rows, r)                     # first entry of each row in the sorted vectors
        i1 <- first + as.vector(rowsum(as.numeric(cw < k1), r))
        i2 <- first + as.vector(rowsum(as.numeric(cw < k2), r))
        out[rows] <- (x[i1] + x[i2])/2
        out
    }
    pool <- function(cg)
    {
        out <- lapply(seq_len(nchr), function(chr)
        {
            R <- colmat("records", chr, cg)
            tot <- rowSums(R)
            x <- if(stat=="median") wmedianRows(colmat("median_qwidth", chr, cg), R)
                 else rowSums(colmat("nucleotides", chr, cg))/tot
            x[tot < min_reads_per_bin] <- NA
            attr(x, "depth") <- tot
            x
        })
        depth <- unlist(lapply(out, attr, "depth"))
        out <- lapply(out, function(x) { attr(x, "depth") <- NULL; x })
        allx <- unlist(out)
        med <- median(allx, na.rm=TRUE)
        out <- lapply(out, function(x) { x[!is.finite(x)] <- med; x })
        names(out) <- chrnames
        attr(out, "stat") <- stat
        attr(out, "n_imputed") <- sum(!is.finite(allx))
        attr(out, "median_depth") <- median(depth)
        attr(out, "n_cells") <- length(cg)
        out
    }
    glob <- pool(cells)
    if(verbose) print(paste0("getReadLengthTrack: ", stat, " read length pooled over ", length(cells),
                             " cells; median pooled reads per bin = ", attr(glob, "median_depth"),
                             "; imputed bins = ", attr(glob, "n_imputed")))
    if(is.null(groups)) return(glob)
    grp <- as.character(groups[cells])
    grp[is.na(grp)] <- "ungrouped"
    out <- lapply(split(cells, grp), function(cg)
    {
        tr <- pool(cg)
        if(attr(tr, "median_depth") < 5*min_reads_per_bin)
        {
            if(verbose) print(paste0("getReadLengthTrack: group of ", length(cg),
                                     " cells too sparse (median pooled reads per bin = ", attr(tr, "median_depth"),
                                     ") - using the global track for it"))
            return(glob)
        }
        tr
    })
    if(verbose) print(paste0("getReadLengthTrack: ", length(out), " group track(s): ",
                             paste(names(out), " (", sapply(split(cells, grp), length), " cells)", sep="", collapse=", ")))
    attr(out, "global") <- glob
    attr(out, "groups") <- setNames(grp, cells)
    attr(out, "stat") <- stat
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
