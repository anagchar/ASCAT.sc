## Define variable-width bins from the pooled (pseudobulk) coverage.
##
## Pools the Level-1 (small-window) read counts across all cells per
## chromosome and walks along the windows, cutting a bin boundary each time
## the accumulated pooled count reaches target_reads. Boundaries therefore
## snap to existing window edges, so the resulting bins inherit the mappability
## filtering already baked into lSe and stay aligned with lGCT and every cell's
## coverage track. The returned per-chromosome index groups are consumed by
## treatTrack / treatlSe / treatGCT (via their `groups` argument).
##
## lCTS_list: list over cells, each a per-chromosome list of small-window
##            coverage data.frames (the `lCTS.tumour` objects).
getPseudobulkGroups <- function(lCTS_list, target_reads=1e6)
{
    nchr <- length(lCTS_list[[1]])
    lapply(seq_len(nchr), function(i)
    {
        pooled <- Reduce(`+`, lapply(lCTS_list, function(cell) cell[[i]]$records))
        n <- length(pooled)
        starts <- integer(0)
        ends <- integer(0)
        acc <- 0
        bin_start <- 1L
        for(j in seq_len(n))
        {
            acc <- acc + pooled[j]
            if(acc >= target_reads)
            {
                starts <- c(starts, bin_start)
                ends <- c(ends, j)
                bin_start <- j + 1L
                acc <- 0
            }
        }
        ## Trailing windows that never reach the target are merged into the last
        ## bin (or form a single bin if the whole chromosome is below target).
        ## This also absorbs zero-read stretches, so a homozygous deletion can
        ## not collapse into its own giant bin.
        if(bin_start <= n)
        {
            if(length(ends) == 0)
            {
                starts <- 1L
                ends <- n
            }
            else
            {
                ends[length(ends)] <- n
            }
        }
        list(starts=starts, ends=ends)
    })
}
