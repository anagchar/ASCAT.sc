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

    # Average read length correction
    # Calculate average read length per bin from nucleotides/records 
    allReadLen <- unlist(lapply(lCT, function(x) {
        ifelse(x$records > 0, x$nucleotides / x$records, NA) # set to NA if no reads are found
    }))

    # Bivariate correction with GC and read length
    valid_idx <- !is.na(allReadLen) & !is.na(allGC) & !is.na(allRec)

    smoothT <- myloess(LL=sum(valid_idx),
                       allRec[valid_idx] ~ allGC[valid_idx] + allReadLen[valid_idx],
                       degree=2,
                       normalize=TRUE)

    # Initialize with NA for invalid bins
    fitted_vals <- rep(NA, length(allRec))
    residual_vals <- rep(NA, length(allRec))

    fitted_vals[valid_idx] <- smoothT$fitted
    residual_vals[valid_idx] <- smoothT$residuals

    # For invalid bins, use uncorrected values
    residual_vals[!valid_idx] <- allRec[!valid_idx]
    fitted_vals[!valid_idx] <- 0

    smoothT <- list(fitted=fitted_vals, residuals=residual_vals)

    for(i in 1:length(lCT))
    {
        lCT[[i]] <- cbind(lCT[[i]],smoothT$fitted[starts[i]:ends[i]],
                          smoothT$residuals[starts[i]:ends[i]])
        colnames(lCT[[i]])[(ncol(lCT[[i]])-1):ncol(lCT[[i]])] <- c("fitted","smoothed")
    }
    return(lCT)
}
