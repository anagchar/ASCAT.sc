combineDiploid <- function(lNormals)
{
    nlNormals <- lNormals[[1]]
    if(length(lNormals)>1)
        for(i in 2:length(lNormals))
        {
            for(j in 1:length(lNormals[[1]]))
            {
                ## sum every per-bin count column (records, nucleotides, read-length bands n_len_*)
                cols <- intersect(c("records","nucleotides",grep("^n_len_",colnames(nlNormals[[j]]),value=TRUE)),
                                  colnames(lNormals[[i]][[j]]))
                for(cc in cols)
                    nlNormals[[j]][[cc]] <- nlNormals[[j]][[cc]]+lNormals[[i]][[j]][[cc]]
            }
        }
    nlNormals
}
