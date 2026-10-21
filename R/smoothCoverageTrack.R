smoothCoverageTrack <- function(lCT,
                                lSe,
                                lGCT,
                                lNormals=NULL,
                                method=c("loess",
                                         "lowess"))
{
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

    # Read-length covariate: median read length per bin (median_qwidth)
    if (all(sapply(lCT, function(x) "median_qwidth" %in% colnames(x))))
        allReadLen <- unlist(lapply(lCT, function(x) x$median_qwidth))
    else
    {
        message("smoothCoverageTrack: no median_qwidth column, using nucleotides/records as read length covariate")
        allReadLen <- unlist(lapply(lCT, function(x) {
            ifelse(x$records > 0, x$nucleotides / x$records, NA) # set to NA if no reads are found
        }))
    }

    # Bivariate correction with GC and read length
    valid_idx <- !is.na(allReadLen) & !is.na(allGC) & !is.na(allRec)
    d <- data.frame(y=allRec, gc=allGC, rl=allReadLen)
    smoothT <- myloess(LL=sum(valid_idx),
                       y ~ gc + rl,
                       data=d[valid_idx,],
                       degree=2,
                       normalize=TRUE)

    # Initialize with NA for invalid bins
    fitted_vals <- rep(NA, length(allRec))
    residual_vals <- rep(NA, length(allRec))

    fitted_vals[valid_idx] <- smoothT$fitted
    residual_vals[valid_idx] <- smoothT$residuals

    # For invalid bins, use uncorrected values
    fitted_vals[!valid_idx] <- allRec[!valid_idx]
    residual_vals[!valid_idx] <- 0

    # Bins without reads: predict at their GC and the cells's median read length
    noread_idx = is.na(allReadLen) & !is.na(allGC) & !is.na(allRec)
    if (any(noread_idx)) 
    {
        gcrange = range(allGC[valid_idx])
        pred = predict(smoothT, newdata=data.frame(gc=pmin(pmax(allGC[noread_idx], gcrange[1]), gcrange[2]), rl=median(allReadLen[valid_idx])))
        fitted_vals[noread_idx] <- pred
        residual_vals[noread_idx] <- allRec[noread_idx] - pred
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