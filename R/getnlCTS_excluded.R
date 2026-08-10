getnlCTS_excluded <- function(nlCTS.tumour, lInds)
{
    nnlCTS.tumour <- lapply(1:length(nlCTS.tumour),function(chr)
    {
        if(length(lInds[[chr]])>0)
            ## drop=FALSE: PCA panel members are single-column data.frames
            ## (records only) and would otherwise collapse to a plain vector,
            ## breaking y$records in smoothNormals for every cell.
            return(nlCTS.tumour[[chr]][-lInds[[chr]],,drop=FALSE])
        nlCTS.tumour[[chr]]
    })
    names(nnlCTS.tumour) <- names(nlCTS.tumour)
    nnlCTS.tumour
}

