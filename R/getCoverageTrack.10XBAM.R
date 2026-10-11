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
                                     tag=c("CB","CR"),
                                     length_breaks=c(150, 250, 500, 1000))
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
                        what=c("qwidth"))
    coverageTrack <- scanBam(bamPath, param = sbp)
    .guessBarcodes <- function(cT, pcchromosome=pcchromosome)
    {
        counts <- table(unlist(lapply(cT,function(x) unique(x$tag[[tag[1]]]))))
        names(counts)[counts>length(cT)*pcchromosome]
    }
    if(is.null(barcodes)) barcodes <- .guessBarcodes(coverageTrack, pcchromosome = pcchromosome)
    ## Read-length (qwidth) bands [0,b1) [b1,b2) ... [bk,Inf): one read-count column per band and
    ## barcode (n_len_<lower>_<upper>), so that a short-read fraction can be computed later at any
    ## band boundary (getReadLengthTrack) without re-reading the BAM.
    length_breaks <- sort(unique(as.numeric(length_breaks)))
    band_breaks <- c(0, length_breaks, Inf)
    band_names <- paste0("n_len_", band_breaks[-length(band_breaks)], "_", band_breaks[-1])
    .countBarcodes <- function(binData, barcodesUNIQUE)
    {
        f <- factor(binData$tag[[tag[1]]], levels=barcodesUNIQUE)   # reads of other barcodes -> NA, ignored
        qw <- as.numeric(binData$qwidth)
        counts <- tabulate(f, nbins=length(barcodesUNIQUE))
        names(counts) <- barcodesUNIQUE
        nucleotides <- vapply(split(qw, f), sum, numeric(1))        # split() keeps empty levels (sum 0)
        band <- cut(qw, breaks=band_breaks, right=FALSE, labels=band_names)
        bands <- unclass(table(f, band))[barcodesUNIQUE, band_names, drop=FALSE]
        list(counts=counts, nucleotides=nucleotides, bands=bands)
    }
    .formatOutput <- function(cT, chr, starts, ends)
    {
        out <- lapply(cT, .countBarcodes, barcodes)
        nb <- length(barcodes)
        asmat <- function(x) matrix(x, nrow=nb, dimnames=list(barcodes, NULL))   # barcodes x bins
        recs <- asmat(unlist(lapply(out, `[[`, "counts")))
        nucs <- asmat(unlist(lapply(out, `[[`, "nucleotides")))
        bands <- lapply(band_names, function(bn) asmat(unlist(lapply(out, function(o) o$bands[, bn]))))
        names(bands) <- band_names
        res <- lapply(barcodes, function(bc)
        {
            df <- data.frame(space=chr,
                             start=starts,
                             end=ends,
                             width=ends-starts,
                             file=paste0(basename(bamPath), "_", bc),
                             records=as.vector(recs[bc, ]),
                             nucleotides=as.vector(nucs[bc, ]))
            for(bn in band_names) df[[bn]] <- as.vector(bands[[bn]][bc, ])
            df
        })
        names(res) <- barcodes
        res
    }
    ## qwidth quantiles of the reads of the kept barcodes (diagnostic to choose short_threshold)
    .qwidthQuantiles <- function(cT)
    {
        qw <- unlist(lapply(cT, function(x) x$qwidth[x$tag[[tag[1]]] %in% barcodes]))
        probs <- c(.05,.25,.5,.75,.95)
        if(length(qw)==0) return(setNames(rep(NA_real_, length(probs)), paste0(100*probs, "%")))
        quantile(as.numeric(qw), probs, na.rm=TRUE)
    }
    qwidth_quantiles <- .qwidthQuantiles(coverageTrack)
    coverageTrack <- .formatOutput(coverageTrack,chr,starts,ends)
    return(list(coverageTrack=coverageTrack,
                barcodes=barcodes,
                qwidth_quantiles=qwidth_quantiles,
                length_breaks=length_breaks))
}
