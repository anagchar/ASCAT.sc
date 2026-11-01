segmentTrack <- function (covtrack,
                          chr,
                          starts,
                          ends = NA,
                          sd = 0,
                          transform=FALSE,
                          min.width = 5,
                          ALPHA=0.01,
                          SBDRY=NULL,
                          smooth=TRUE,
                          segmentation=c("cbs","pcf","quickseg"),
                          penalty=NULL,
                          sigma=NULL)
{
    ## segmentation: "cbs" (DNAcopy, below), "pcf" (copynumber::pcf, gamma=penalty in units of the noise variance
    ## sigma^2) or "quickseg" (penalty in units of the noise sd sigma, default 7, see segmentTrack_quickseg).
    ## sigma=NULL: estimated from covtrack. The pcf default gamma=20 gives at most as many false breakpoints as
    ## CBS at alpha=0.01 on simulated single-cell tracks; 1/ALPHA (=100, the multipcf translation) leaves most
    ## CNAs of a single cell unsegmented (diagnostics/sim_segmentation_methods.R).
    segmentation <- match.arg(segmentation)
    if(segmentation=="pcf")
        return(segmentTrack_pcf(covtrack,
                                chr=chr,
                                starts=starts,
                                ends=ends,
                                sd=sd,
                                transform=transform,
                                min.width=min.width,
                                ALPHA=ALPHA,
                                gamma=if(is.null(penalty)) 20 else penalty,
                                sigma=sigma))
    ## QuickSeg keeps its own default min.width (3) unless one is given
    if(segmentation=="quickseg")
        return(segmentTrack_quickseg(covtrack,
                                     chr=chr,
                                     starts=starts,
                                     ends=ends,
                                     penalty=penalty,
                                     sigma=sigma,
                                     min.width=if(missing(min.width)) 3 else min.width))
    require(DNAcopy)
    covtrack <- covtrack * rnorm(length(covtrack),
                                 mean = 1,
                                 sd = sd)
    if(!transform)
        cna <- CNA(covtrack,
                   chr = rep(chr, length(covtrack)),
                   maploc = starts,
                   data.type = "logratio")
    if(transform)
        cna <- CNA(sqrt(2^covtrack +3/8),
                   chr = rep(chr, length(covtrack)),
                   maploc = starts,
                   data.type = "logratio")
    if(smooth)
        cna <- smooth.CNA(cna)
    if(is.null(SBDRY))
        capture.output(cna <- DNAcopy::segment(cna, min.width = min.width,alpha=ALPHA))
    else
        capture.output(cna <- DNAcopy::segment(cna, min.width = min.width,alpha=ALPHA, sbdry=SBDRY))
    if(transform)
    {
        cna. <- cna$output
        indices <- lapply(1:nrow(cna.),function(x)
        {
            if(x==1)
            {
                return(1:cna.[x,"num.mark"])
            }
            else
            {
                offset <- cumsum(cna.[1:(x-1),"num.mark"])
                return((offset+1):(offset+cna.[x,"num.mark"]))
            }
        })
        for(i in 1:nrow(cna.))
        {
            cna$output[i,"seg.mean"] <- mean(covtrack[indices[[i]]],na.rm=T)
        }
    }
    cna
}
