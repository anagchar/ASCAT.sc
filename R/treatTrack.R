treatTrack <- function(lCTS, window)
{
    nlCTS <- lapply(lCTS,function(df)
    {
        nr <- nrow(df)
        starts <- getstartends(end=nr,window=window)$starts
        ends <- getstartends(end=nr,window=window)$ends
        sumCol <- function(col) sapply(1:length(starts),function(x) sum(df[starts[x]:ends[x],col]))
        ndf <- data.frame(space=df[starts,"space"],
                          start=df[starts,"start"],
                          end=df[ends,"end"],
                          width=df[starts,"width"]*window,
                          file=df[starts,"file"],
                          records=sumCol("records"),
                          nucleotides=sumCol("nucleotides"))
        ## haplotype counts from the HP tag (only present when hp_tag was set)
        if(all(c("hap1","hap2") %in% colnames(df)))
        {
            ndf$hap1 <- sumCol("hap1")
            ndf$hap2 <- sumCol("hap2")
        }
        ndf
    })
    names(nlCTS) <- names(lCTS)
    nlCTS
}
