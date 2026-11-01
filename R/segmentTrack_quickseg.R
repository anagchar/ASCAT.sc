## QuickSeg segmentation (Schlotmann et al., bioRxiv 2026, doi:10.64898/2026.08.31.747916): exact minimiser of
##     sum(|covtrack - segment value|) + penalty*sigma*(number of breakpoints)
## by dynamic programming over candidate segment values. Under the absolute error the best value of a segment
## is its median, so the observed values are the only candidates. Uses quickseg::segment() (C++,
## https://github.com/Balthasar-eu/Rquickseg) when that package is installed, otherwise .quickseg() below,
## which follows the same recursion and tie-breaking.
##
## penalty is in units of the noise sd sigma (see getTrackNoise), so that one value fits cells of different
## coverage, as copynumber::pcf scales its gamma by the noise variance. getTrackForAll passes one sigma per
## cell; sigma=NULL estimates it from covtrack alone and sigma=1 gives the raw penalty of the paper. The
## default 7 gives about as many false breakpoints as CBS at alpha=0.01 on simulated tracks
## (diagnostics/sim_segmentation_methods.R). Segments shorter than min.width bins (as min.width of CBS and
## kmin of PCF) are merged into a neighbour, which keeps single outlier bins, such as the wider last window
## of each chromosome (getstartends), from becoming segments; min.width=1 gives the exact QuickSeg optimum.
## The default 3 still keeps the 4-bin focal events that CBS (min.width=5) calls on simulated cells
## (diagnostics/compare_segmentations.R). Returns the list(data, output) of segmentTrack(). Non-finite
## values are skipped, as CBS does.
segmentTrack_quickseg <- function(covtrack,
                                  chr,
                                  starts,
                                  ends=NA,
                                  penalty=NULL,
                                  sigma=NULL,
                                  min.width=3)
{
    if(is.null(penalty)) penalty <- 7
    if(is.null(sigma)) sigma <- getTrackNoise(covtrack)
    ok <- which(is.finite(covtrack))
    y <- covtrack[ok]
    if(length(y)==0)
        return(list(data=NULL,
                    output=data.frame(ID="Sample.ID",
                                      chrom=chr,
                                      loc.start=starts[1],
                                      loc.end=starts[length(starts)],
                                      num.mark=0,
                                      seg.mean=NA)))
    if(requireNamespace("quickseg", quietly=TRUE))
        segs <- as.matrix(quickseg::segment(y, penalty*sigma)[, c("Startidx","Endidx")])+1L
    else
        segs <- .quickseg(y, penalty*sigma)
    segs <- .mergeShort(segs, y, min.width)
    list(data=NULL,
         output=data.frame(ID="Sample.ID",
                           chrom=chr,
                           loc.start=starts[ok][segs[,1]],
                           loc.end=starts[ok][segs[,2]],
                           num.mark=segs[,2]-segs[,1]+1,
                           seg.mean=sapply(seq_len(nrow(segs)), function(i) mean(y[segs[i,1]:segs[i,2]]))))
}

## Noise sd of a logR track (a vector, or a list of per-chromosome vectors): MAD of the differences between
## neighbouring bins / sqrt(2), which a few breakpoints do not move. Zero differences (runs of empty bins,
## set to 0 by smoothCoverageTrack) are left out. Falls back to 1 when there is nothing to estimate from.
getTrackNoise <- function(x)
{
    if(!is.list(x)) x <- list(x)
    d <- unlist(lapply(x, function(v) diff(v[is.finite(v)])))
    d <- d[d!=0]
    sigma <- if(length(d)>1) mad(d)/sqrt(2) else NA
    if(!is.finite(sigma) || sigma<=0) 1 else sigma
}

## R version of the QuickSeg recursion. score[j] is the cost of the best segmentation of y[1:i] whose last
## segment takes values[j]; state j starts a new segment at i when the best score at i-1 plus the penalty is
## lower than its own. Instead of the n x |values| decision bits of the C++ version, this keeps the start of
## the last segment of the best segmentation at each position (O(n) memory). Returns 1-based start/end indices.
.quickseg <- function(y, penalty, values=sort(unique(y)))
{
    n <- length(y)
    score <- abs(y[1]-values)
    segstart <- rep(1L, length(values))
    best <- bestStart <- integer(n)
    best[1] <- which.min(score)
    bestStart[1] <- 1L
    for(i in seq_len(n)[-1])
    {
        restart <- score[best[i-1]]+penalty
        new <- restart<score
        score[new] <- restart
        segstart[new] <- i
        score <- score+abs(y[i]-values)
        best[i] <- which.min(score)
        bestStart[i] <- segstart[best[i]]
    }
    ends <- n
    while(bestStart[ends[1]]>1) ends <- c(bestStart[ends[1]]-1L, ends)
    cbind(start=c(1L, ends[-length(ends)]+1L), end=ends)
}

## Merge segments (rows of start/end indices into y) shorter than min.width bins, the shortest first, each into
## the neighbour with the closest median.
.mergeShort <- function(segs, y, min.width)
{
    while(nrow(segs)>1)
    {
        len <- segs[,2]-segs[,1]+1
        i <- which.min(len)
        if(len[i]>=min.width) break
        nbr <- c(i-1, i+1)
        nbr <- nbr[nbr>=1 & nbr<=nrow(segs)]
        med <- sapply(c(i, nbr), function(k) median(y[segs[k,1]:segs[k,2]]))
        j <- nbr[which.min(abs(med[-1]-med[1]))]
        segs[min(i, j), 2] <- segs[max(i, j), 2]
        segs <- segs[-max(i, j), , drop=FALSE]
    }
    segs
}
