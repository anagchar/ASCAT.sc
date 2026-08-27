getReadN50 <- function(bamPath,
                       isDuplicate=FALSE,
                       isSecondaryAlignment=FALSE,
                       isSupplementaryAlignment=FALSE,
                       isUnmappedQuery=FALSE)
{
    require(Rsamtools)
    sbp <- ScanBamParam(flag=scanBamFlag(isDuplicate=isDuplicate,
                                         isSecondaryAlignment=isSecondaryAlignment,
                                         isSupplementaryAlignment=isSupplementaryAlignment,
                                         isUnmappedQuery=isUnmappedQuery),
                        what="qwidth",
                        mapqFilter=30)
    qwidth <- scanBam(bamPath, param=sbp)[[1]]$qwidth
    qwidth <- sort(as.numeric(qwidth[!is.na(qwidth)]), decreasing=TRUE)
    qwidth[which(cumsum(qwidth)>=sum(qwidth)/2)[1]]
}
