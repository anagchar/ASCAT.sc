smoothCoverageTrack <- function(lCT,
                                lSe,
                                lGCT,
                                lNormals=NULL,
                                lRLT=NULL,
                                readlen=c("cell","pseudobulk","none"),
                                method=c("loess",
                                         "lowess"))
{
    ## readlen: read-length correction of the coverage on top of the GC loess
    ##   "cell":       per-track mean read length per bin (nucleotides/records) as 2nd covariate;
    ##                 fine for bulk samples, too noisy for single cells (few reads per bin)
    ##   "pseudobulk": bin-level covariate lRLT (list per chromosome, same bins as lGCT) pooled
    ##                 across cells by getReadLengthTrack(); the loess is still fitted per track
    ##   "none":       GC-only loess (upstream behaviour)
    readlen <- match.arg(readlen)
    if(method[1]!="loess") warning("smoothCoverageTrack: only method='loess' is supported, using loess")
    allRec <- unlist(lapply(lCT,function(x) log2(x$records+1)))
    if(!is.null(lNormals))
    {
        if(!"records"%in%names(lNormals[[1]]))
            allRec <- smoothNormals(allRec, lNormals)
        else
            allRec <- allRec-unlist(lapply(lNormals,function(x) log2(x$records+1)))
    }
    allGC <- unlist(lGCT)
    starts <- c(0,cumsum(sapply(lCT,nrow)[-c(length(lCT))]))+1
    ends <- cumsum(sapply(lCT,nrow))
    if(readlen=="pseudobulk")
    {
        if(is.null(lRLT))
        {
            warning("smoothCoverageTrack: readlen='pseudobulk' but no lRLT given - using GC-only correction")
            readlen <- "none"
        }
        else if(length(unlist(lRLT))!=length(allRec))
            stop("smoothCoverageTrack: lRLT does not match the bins of lCT")
    }
    ## loess fitted on the bins with finite response and covariates; the other bins (e.g. bins without
    ## reads, hence without a read length in 'cell' mode) get the prediction at the track's median
    ## covariate value, instead of residual 0 (which would report them as average coverage)
    .fitLoess <- function(df, form, ...)
    {
        ok <- apply(is.finite(as.matrix(df)), 1, all)
        fit <- myloess(LL=sum(ok), form, data=df[ok,,drop=FALSE], ...)
        fitted <- rep(NA_real_, nrow(df))
        fitted[ok] <- fit$fitted
        if(any(!ok))
        {
            nd <- df[!ok, setdiff(colnames(df), "y"), drop=FALSE]
            for(cc in colnames(nd)) nd[[cc]][!is.finite(nd[[cc]])] <- median(df[[cc]][ok])
            p <- predict(fit, newdata=nd)
            p[!is.finite(p)] <- median(fitted, na.rm=TRUE)
            fitted[!ok] <- p
        }
        fitted
    }
    df <- data.frame(y=allRec, gc=allGC)
    if(readlen=="none")
        fitted <- .fitLoess(df, y ~ gc)
    if(readlen=="pseudobulk")
    {
        df$rl <- unlist(lRLT)
        fitted <- .fitLoess(df, y ~ gc + rl, degree=2, normalize=TRUE)
    }
    if(readlen=="cell")
    {
        df$rl <- unlist(lapply(lCT,function(x) ifelse(x$records>0, x$nucleotides/x$records, NA)))
        fitted <- .fitLoess(df, y ~ gc + rl, degree=2, normalize=TRUE)
    }
    residuals <- df$y-fitted
    residuals[!is.finite(residuals)] <- 0
    smoothT <- list(fitted=fitted, residuals=residuals)
    for(i in 1:length(lCT))
    {
        lCT[[i]] <- cbind(lCT[[i]],smoothT$fitted[starts[i]:ends[i]],
                          smoothT$residuals[starts[i]:ends[i]])
        colnames(lCT[[i]])[(ncol(lCT[[i]])-1):ncol(lCT[[i]])] <- c("fitted","smoothed")
    }
    return(lCT)
}
