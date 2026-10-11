getstartends <- function(start=1, end, window)
{
    starts <- seq(start,end,window)
    starts <- starts[-c(length(starts))]
    ## fewer bins than one window: a single window covering everything
    ## (previously starts was empty and callers crashed on starts[1]:ends[1])
    if(length(starts)==0) return(list(starts=start,ends=end))
    ends <- c(starts[-c(1)]-1,end)
    list(starts=starts,ends=ends)
}
