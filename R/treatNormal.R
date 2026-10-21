treatNormal <- function(lCTS, window)
{
    nlCTS <- lapply(lCTS,function(df)
    {
        cat(".")
        nr <- nrow(df)
        starts <- getstartends(end=nr,window=window)$starts
        ends <- getstartends(end=nr,window=window)$ends
        ndf <- data.frame(space=df[starts,"space"],
                          start=df[starts,"start"],
                          end=df[ends,"end"],
                          width=df[starts,"width"]*window,
                          file=df[starts,"file"],
                          records=sapply(1:length(starts),function(x) sum(df[starts[x]:ends[x],"records"])),
                          nucleotides=sapply(1:length(starts),function(x) sum(df[starts[x]:ends[x],"nucleotides"])))
        if ("median_qwidth" %in% colnames(df))
            ndf$median_qwidth <- sapply(1:length(starts), function(x) median(rep(df[starts[x]:ends[x], "median_qwidth"], df[starts[x]:ends[x], "records"]), na.rm=TRUE))
        ## extra count columns are summed; median_qwidth (a median, set above) is not a count
        for(cc in setdiff(colnames(df)[-(1:7)], "median_qwidth"))
            ndf[[cc]] <- sapply(1:length(starts),function(x) sum(df[starts[x]:ends[x],cc]))
        ndf
    })
    names(nlCTS) <- names(lCTS)
    nlCTS
}
