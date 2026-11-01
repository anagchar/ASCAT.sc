segmentTrack_pcf <- function (covtrack,
                              chr,
                              starts,
                              ends = NA,
                              sd = 0,
                              transform=FALSE,
                              ismedian=FALSE,
                              min.width = 5,
                              ALPHA=0.01,
                              SBDRY=NULL,
                              smooth=TRUE,
                              gamma=NULL,
                              sigma=NULL)
{
    if(is.null(gamma)) gamma <- 1/ALPHA ## rough translation from CBS pval to pcf penalty
    ## pcf scales gamma by the noise variance of the data it gets (here one chromosome); a given sigma (e.g. of
    ## the whole cell, see getTrackNoise) is used instead, which is more stable on short chromosomes
    if(!is.null(sigma)) gamma <- gamma*sigma^2
    suppressPackageStartupMessages(require(copynumber))
    covtrack <- covtrack * rnorm(length(covtrack),
                                 mean = 1,
                                 sd = sd)
    if(transform)
        covtrack <- sqrt(2^covtrack +3/8)
    x <- pcf(data.frame(chr=rep(chr,length(starts)),
                        positions=starts,
                        sample1=covtrack),
             gamma=gamma,
             normalize=is.null(sigma),
             verbose=F,
             kmin=min.width)
    cna <- list(data=NULL,
                output=data.frame(ID="Sample.ID",
                                  chrom=x[,"chrom"],
                                  loc.start=as.numeric(as.character(x[,"start.pos"])),
                                  loc.end=as.numeric(as.character(x[,"end.pos"])),
                                  num.mark=as.numeric(as.character(x[,"n.probes"])),
                                  seg.mean=as.numeric(as.character(x[,"mean"]))))
    if(ismedian)
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
                offset <- sum(cna.[1:(x-1),"num.mark"])
                return((offset+1):(offset+cna.[x,"num.mark"]))
            }
        })
        for(i in 1:nrow(cna.))
        {
            cna$output[i,"seg.mean"] <- median(covtrack[indices[[i]]],na.rm=T)
        }
    }
    cna
}
