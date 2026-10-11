treatTrack <- function(lCTS, window)
{
    nlCTS <- lapply(lCTS,function(df)
    {
        nr <- nrow(df)
        se <- getstartends(end=nr,window=window)
        starts <- se$starts
        ends <- se$ends
        grp <- rep.int(seq_along(starts), ends-starts+1)
        sumcol <- function(cc) as.vector(rowsum(as.numeric(df[[cc]]), grp, reorder=TRUE))
        ndf <- data.frame(space=df[starts,"space"],
                          start=df[starts,"start"],
                          end=df[ends,"end"],
                          width=df[starts,"width"]*window,
                          file=df[starts,"file"],
                          records=sumcol("records"),
                          nucleotides=sumcol("nucleotides"))
        ## extra per-bin count columns (read-length bands n_len_*) are summed the same way
        extra <- setdiff(colnames(df), colnames(ndf))
        extra <- extra[vapply(extra, function(cc) is.numeric(df[[cc]]), logical(1))]
        for(cc in extra) ndf[[cc]] <- sumcol(cc)
        ndf
    })
    names(nlCTS) <- names(lCTS)
    nlCTS
}
