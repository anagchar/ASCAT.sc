treatTrack <- function(lCTS, window)
{
    nlCTS <- lapply(lCTS,function(df)
    {
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
                          if (all(c("hap1", "hap2") %in% colnames(df))) {
                            ndf$hap1 <- sapply(1:length(starts),function(x) sum(df[starts[x]:ends[x],"hap1"]))
                            ndf$hap2 <- sapply(1:length(starts),function(x) sum(df[starts[x]:ends[x],"hap2"]))
                          }
        ndf
    })
    names(nlCTS) <- names(lCTS)
    nlCTS
}
