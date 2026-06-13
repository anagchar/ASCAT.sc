treatGCT <- function(lGCT, window=ceiling(as.numeric(WINDOW)/as.numeric("10000")), groups=NULL)
{
    nlGCT <- lapply(seq_along(lGCT),function(i)
    {
        gc <- lGCT[[i]]
        l <- length(gc)
        se <- if(is.null(groups)) getstartends(end=l,window=window) else groups[[i]]
        sapply(seq_along(se$starts),function(x) mean(gc[se$starts[x]:se$ends[x]]))
    })
    names(nlGCT) <- names(lGCT)
    nlGCT
}
