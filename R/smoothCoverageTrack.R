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
    ##   "cell":       the track's own median read length per bin (median_qwidth; nucleotides/records when
    ##                 the column is missing) as 2nd covariate; fine for bulk samples, noisy for single
    ##                 cells (few reads per bin)
    ##   "pseudobulk": bin-level read length lRLT (list per chromosome, same bins as lGCT) pooled across
    ##                 cells by getReadLengthTrack(); the loess is still fitted per track
    ##   "none":       GC-only loess (upstream behaviour)
    readlen <- match.arg(readlen)
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
    d <- data.frame(y=allRec, gc=allGC)
    if(readlen=="pseudobulk")
        d$rl <- unlist(lRLT)
    if(readlen=="cell")
    {
        # Read-length covariate: median read length per bin (median_qwidth)
        if (all(sapply(lCT, function(x) "median_qwidth" %in% colnames(x))))
            d$rl <- unlist(lapply(lCT, function(x) x$median_qwidth))
        else
        {
            message("smoothCoverageTrack: no median_qwidth column, using nucleotides/records as read length covariate")
            d$rl <- unlist(lapply(lCT, function(x) {
                ifelse(x$records > 0, x$nucleotides / x$records, NA) # set to NA if no reads are found
            }))
        }
    }

    # Loess on the bins with a coverage value, a GC value and (if used) a read length
    valid_idx <- apply(is.finite(as.matrix(d)), 1, all)
    form <- if(readlen=="none") y ~ gc else y ~ gc + rl
    smoothT <- myloess(LL=sum(valid_idx),
                       form,
                       data=d[valid_idx,],
                       degree=2,
                       normalize=TRUE)

    fitted_vals <- rep(NA, length(allRec))
    residual_vals <- rep(NA, length(allRec))
    fitted_vals[valid_idx] <- smoothT$fitted
    residual_vals[valid_idx] <- smoothT$residuals

    # Bins without reads ("cell" mode: no read length): predict at their GC (within the fitted GC range)
    # and the track's median read length, so they report lower coverage than expected rather than 0
    noread_idx = !valid_idx & is.finite(d$y) & is.finite(d$gc)
    if (any(noread_idx))
    {
        nd <- d[noread_idx, setdiff(colnames(d), "y"), drop=FALSE]
        for(cc in colnames(nd))
        {
            rng <- range(d[[cc]][valid_idx])
            nd[[cc]][!is.finite(nd[[cc]])] <- median(d[[cc]][valid_idx])
            nd[[cc]] <- pmin(pmax(nd[[cc]], rng[1]), rng[2])
        }
        pred = predict(smoothT, newdata=nd)
        pred[!is.finite(pred)] <- median(fitted_vals, na.rm=TRUE)
        fitted_vals[noread_idx] <- pred
        residual_vals[noread_idx] <- d$y[noread_idx] - pred
    }

    # Other invalid bins (no GC or no coverage value): uncorrected values, as before
    other_idx = !valid_idx & !noread_idx
    fitted_vals[other_idx] <- allRec[other_idx]
    residual_vals[other_idx] <- 0

    smoothT <- list(fitted=fitted_vals, residuals=residual_vals)

    cat("\n--- Applying correction to coverage tracks ---\n")
    for(i in 1:length(lCT))
    {
        lCT[[i]] <- cbind(lCT[[i]],smoothT$fitted[starts[i]:ends[i]],
                          smoothT$residuals[starts[i]:ends[i]])
        colnames(lCT[[i]])[(ncol(lCT[[i]])-1):ncol(lCT[[i]])] <- c("fitted","smoothed")
    }
    return(lCT)
}
