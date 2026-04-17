getCoverageTrack.Fix <-
function(bamPath,
                                 chr,
                                 lengthChr,
                                 step,
                                 CHRSTRING="")
{
    require(Rsamtools)
    require(GenomicRanges)
    divideChr <- seq(0, lengthChr, step)
    starts <- divideChr[-c(length(divideChr))] + 1
    ends <- divideChr[-c(1)]
    sbp <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE),
        which = GRanges(paste0(CHRSTRING,chr), IRanges(starts, ends)), what = "qwidth")
    coverageTrack <- countBam(bamPath, param = sbp)
    
    counts <- sapply(coverageTrack, function(x) length(x$qwidth))
    nucleotides <- sapply(coverageTrack, function(x) sum(x$qwidth))

    return(list(
        counts = counts,
        nucleotides = nucleotides))
    return(coverageTrack)
  }
