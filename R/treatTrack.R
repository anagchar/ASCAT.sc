treatTrack <- function(lCTS, window=NULL, groups=NULL)
{
    nlCTS <- lapply(seq_along(lCTS),function(i)
    {
        df <- lCTS[[i]]
        nr <- nrow(df)
        se <- if(is.null(groups)) getstartends(end=nr,window=window) else groups[[i]]
        starts <- se$starts
        ends <- se$ends
        w <- if(is.null(groups)) df[starts,"width"]*window else (ends-starts+1)*df[starts,"width"]
        ndf <- data.frame(space=df[starts,"space"],
                          start=df[starts,"start"],
                          end=df[ends,"end"],
                          width=w,
                          file=df[starts,"file"],
                          records=sapply(seq_along(starts),function(x) sum(df[starts[x]:ends[x],"records"])),
                          nucleotides=sapply(seq_along(starts),function(x) sum(df[starts[x]:ends[x],"nucleotides"])))
        ndf
    })
    names(nlCTS) <- names(lCTS)
    nlCTS
}
