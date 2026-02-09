getCoverageTrack.10XBAM <- function (bamPath,
                                     chr,
                                     starts,
                                     ends,
                                     chrstring = "",
                                     isDuplicate=F,
                                     isSecondaryAlignment=F,
                                     isNotPassingQualityControls=NA,
                                     isUnmappedQuery=NA,
                                     mapqFilter=30,
                                     barcodes=NULL,
                                     pcchromosome=0.8,
                                     tag=c("CB","CR"))
{
    suppressPackageStartupMessages(require(Rsamtools))
    suppressPackageStartupMessages(require(GenomicRanges))
    sbp <- ScanBamParam(flag = scanBamFlag(isDuplicate=isDuplicate,
                                           isSecondaryAlignment=isSecondaryAlignment,
                                           isNotPassingQualityControls=isNotPassingQualityControls,
                                           isUnmappedQuery=isUnmappedQuery),
                        which = GRanges(paste0(chrstring, chr),
                                        IRanges(starts,ends)),
                        mapqFilter=mapqFilter,
                        tag=tag[1],
                        what="qwidth")
    coverageTrack <- scanBam(bamPath, param = sbp)
    .guessBarcodes <- function(cT, pcchromosome=pcchromosome)
    {
        counts <- table(unlist(lapply(cT,function(x) unique(x$tag[[tag[1]]]))))
        names(counts)[counts>length(cT)*pcchromosome]
    }
    if(is.null(barcodes)) barcodes <- .guessBarcodes(coverageTrack, pcchromosome = pcchromosome)
    .countBarcodesAndNucleotides <- function(binData, barcodesUNIQUE)
    {
        bc <- binData$tag[[tag[1]]]
        qw <- binData$qwidth
        counts <- sapply(barcodesUNIQUE, function(b) sum(bc == b, na.rm=TRUE))
        nucleotides <- sapply(barcodesUNIQUE, function(b) sum(qw[bc == b], na.rm=TRUE))
        list(counts=counts, nucleotides=nucleotides)
    }
    .formatCounts <- function(counts, nucleotides, chr, start, end, bamPath, bc)
    {
        samplename = paste0(basename(bamPath), "_", bc)
        data.frame(space=NA,
                   start=start,
                   end=end,
                   width=NA,
                   file=samplename,
                   records=counts,
                   nucleotides=nucleotides)
    }
    .formatOutput <- function(cT, chr, starts, ends)
    {
        out <- lapply(cT, .countBarcodesAndNucleotides, barcodes)
        out <- lapply(barcodes,
                      function(bc)
                      {
                          df <- do.call("rbind",
                                        lapply(1:length(out),
                                               function(i)
                                                   .formatCounts(out[[i]]$counts[bc],
                                                                 out[[i]]$nucleotides[bc],
                                                                 chr,starts[i],ends[i], bamPath, bc)))
                          df$width <- df$end-df$start
                          df$space <- chr
                          df
                      })
        names(out) <- barcodes
        out
    }
    coverageTrack <- .formatOutput(coverageTrack,chr,starts,ends)
    return(list(coverageTrack=coverageTrack,
                barcodes=barcodes))
}
