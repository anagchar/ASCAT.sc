sc_excludeBadBins <- function(res)
{
    get_nbins_lSe <- function(x)
    {
        if(is.null(x)) return(NA_integer_)
        sum(sapply(x, function(chr) length(chr$starts)))
    }
    get_nbins_lGCT <- function(x)
    {
        if(is.null(x)) return(NA_integer_)
        sum(sapply(x, length))
    }
    get_nbins_nlCTS <- function(x)
    {
        if(is.null(x)) return(NA_integer_)
        sum(sapply(x, nrow))
    }

    ## pick the Se/GCT set that matches the actual nlCTS bins, otherwise default to lSe/lGCT
    nBins.ref <- NA_integer_
    if(!is.null(res$allTracks) && length(res$allTracks) > 0 && !is.null(res$allTracks[[1]]$nlCTS.tumour))
        nBins.ref <- get_nbins_nlCTS(res$allTracks[[1]]$nlCTS.tumour)
    else if(!is.null(res$lNormals) && length(res$lNormals) > 0)
        nBins.ref <- get_nbins_nlCTS(res$lNormals[[1]])
    else if(!is.null(res$lCTS.normal) && length(res$lCTS.normal) > 0 && !is.null(res$lCTS.normal[[1]]$nlCTS.normal))
        nBins.ref <- get_nbins_nlCTS(res$lCTS.normal[[1]]$nlCTS.normal)

    nBins.nlSe <- get_nbins_lSe(res$nlSe)
    nBins.lSe <- get_nbins_lSe(res$lSe)
    if(is.null(res$nlSe) || (!is.na(nBins.ref) && nBins.nlSe != nBins.ref && !is.null(res$lSe) && nBins.lSe == nBins.ref))
    {
        if(is.null(res$lSe))
            stop("sc_excludeBadBins(): missing both res$nlSe and res$lSe")
        if(!is.null(res$nlSe) && !is.na(nBins.ref) && nBins.nlSe != nBins.ref)
            print("sc_excludeBadBins(): res$nlSe does not match nlCTS bins, using res$lSe")
        res$nlSe <- res$lSe
    }

    nBins.target <- get_nbins_lSe(res$nlSe)
    nBins.nlGCT <- get_nbins_lGCT(res$nlGCT)
    nBins.lGCT <- get_nbins_lGCT(res$lGCT)
    if(is.null(res$nlGCT) || nBins.nlGCT != nBins.target)
    {
        if(!is.null(res$lGCT) && nBins.lGCT == nBins.target)
        {
            if(!is.null(res$nlGCT) && nBins.nlGCT != nBins.target)
                print("sc_excludeBadBins(): res$nlGCT does not match selected bins, using res$lGCT")
            res$nlGCT <- res$lGCT
        }
        else if(is.null(res$nlGCT))
        {
            stop("sc_excludeBadBins(): missing both res$nlGCT and res$lGCT")
        }
        else
        {
            stop("sc_excludeBadBins(): neither res$nlGCT nor res$lGCT matches selected bin definition")
        }
    }
    ## Non-10X: normal BAMs stored under res$lCTS.normal
    ## 10X: normal barcodes stored under res$lNormals (set via normal_barcodes)
    has_bam_normals <- !is.null(res$lCTS.normal) &&
        length(res$lCTS.normal) > 0 &&
        !is.null(res$lCTS.normal[[1]]$nlCTS.normal)
    has_barcode_normals <- !has_bam_normals && !is.null(res$lNormals) && length(res$lNormals) > 0

    allTracks_nms <- names(res$allTracks)

    if(has_bam_normals)
    {
        print("Using normal BAM samples for removal of bad bins")
        allTracks.normal <- lapply(res$lCTS.normal, function(x) list(lCTS = x$nlCTS.normal))
        lInds <- filterBins(allTracks=allTracks.normal, logr=NULL, lSe=res$nlSe)
        res$nlGCT <- getlGCT_excluded(res$nlGCT, lInds)
        res$nlSe <- getlSe_excluded(res$nlSe, lInds)
        if(!is.null(res$lNormals))
            res$lNormals <- lapply(res$lNormals, function(x) getnlCTS_excluded(x, lInds))
        res$allTracks <- lapply(names(res$allTracks),function(x)
        {
            res$allTracks[[x]]$nlCTS.tumour = getnlCTS_excluded(res$allTracks[[x]]$nlCTS.tumour, lInds)
            res$allTracks[[x]]
        })
    }
    else if(has_barcode_normals)
    {
        print(paste0("Using ", length(res$lNormals), " normal barcodes for removal of bad bins"))
        allTracks.normal <- lapply(res$lNormals, function(x) list(lCTS = x))
        lInds <- filterBins(allTracks=allTracks.normal, logr=NULL, lSe=res$nlSe)
        res$nlGCT <- getlGCT_excluded(res$nlGCT, lInds)
        res$nlSe <- getlSe_excluded(res$nlSe, lInds)
        res$lNormals <- lapply(res$lNormals, function(x) getnlCTS_excluded(x, lInds))
        res$allTracks <- lapply(names(res$allTracks),function(x)
        {
            res$allTracks[[x]]$nlCTS.tumour = getnlCTS_excluded(res$allTracks[[x]]$nlCTS.tumour, lInds)
            res$allTracks[[x]]
        })
    }
    else
    {
        print("Using all samples for removal of bad bins")
        allTracks <- lapply(res$allTracks, function(x) {
            list(lCTS = x$nlCTS.tumour)})
        lInds <- filterBins(allTracks=allTracks, logr=NULL, lSe=res$nlSe, fraction=.1, k=11)
        res$nlGCT <- getlGCT_excluded(res$nlGCT, lInds)
        res$nlSe <- getlSe_excluded(res$nlSe, lInds)
        if(!is.null(res$lNormals))
            res$lNormals <- lapply(names(res$lNormals), function(x) getnlCTS_excluded(res$lNormals[[x]], lInds))
        res$allTracks <- lapply(names(res$allTracks),function(x)
        {
            res$allTracks[[x]]$nlCTS.tumour = getnlCTS_excluded(res$allTracks[[x]]$nlCTS.tumour, lInds)
            res$allTracks[[x]]
        })
    }
    names(res$allTracks) <- allTracks_nms
    res
}
