treatlSe <- function(lSe, window=ceiling(as.numeric(WINDOW)/as.numeric("10000")), groups=NULL)
{
    nlSe <- lapply(seq_along(lSe),function(i)
    {
        x <- lSe[[i]]
        l <- length(x$starts)
        se <- if(is.null(groups)) getstartends(end=l,window=window) else groups[[i]]
        list(starts=x$starts[se$starts],
             ends=x$ends[se$ends])
    })
    names(nlSe) <- names(lSe)
    nlSe
}
