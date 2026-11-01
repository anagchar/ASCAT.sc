## QuickSeg segmentation (segmentTrack_quickseg / segmentation="quickseg") and the segmentation option.

## exact optimum of sum(|y - segment median|) + p*(segments-1), O(n^2) optimal partitioning
l1_optimum <- function(y, p)
{
    n <- length(y); F <- c(-p, rep(Inf, n))
    for(t in 1:n) for(s in 1:t) { v <- y[s:t]; F[t+1] <- min(F[t+1], F[s]+sum(abs(v-median(v)))+p) }
    F[n+1]
}
l1_cost <- function(y, segs, p)
    sum(apply(segs, 1, function(s) { v <- y[s[1]:s[2]]; sum(abs(v-median(v))) }))+p*(nrow(segs)-1)

test_that(".quickseg reaches the exact L1 optimum and returns contiguous segments", {
    set.seed(1)
    for(k in 1:150)
    {
        n <- sample(1:30, 1)
        y <- round(rnorm(n)+rep(rnorm(3, 0, 2), length.out=n)[sort(sample(n))], sample(c(1, 8), 1))
        p <- sample(c(0, 0.1, 0.5, 1, 3, 10), 1)
        segs <- unname(.quickseg(y, p))
        expect_equal(segs[1, 1], 1L)
        expect_equal(segs[nrow(segs), 2], n)
        expect_equal(segs[-1, 1], segs[-nrow(segs), 2]+1L)
        expect_equal(l1_cost(y, segs, p), l1_optimum(y, p), tolerance=1e-10)
    }
    ## example of the quickseg package
    expect_equal(unname(.quickseg(c(1.1, 0.9, 1, 1, 3, 3, 3, 3), 1)), cbind(c(1L, 5L), c(4L, 8L)))
})

test_that("segmentTrack(segmentation='quickseg') returns the DNAcopy-like output used downstream", {
    set.seed(2)
    y <- c(rnorm(60), rnorm(30, 3), rnorm(60))
    starts <- seq(1, by=5e5, length.out=150)
    out <- segmentTrack(y, chr="3", starts=starts, segmentation="quickseg")
    expect_null(out$data)
    o <- out$output
    expect_equal(colnames(o), c("ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean"))
    expect_equal(nrow(o), 3)
    expect_true(all(o$chrom=="3"))
    expect_equal(sum(o$num.mark), 150)
    expect_equal(o$loc.start[1], starts[1])
    expect_equal(o$loc.end[3], starts[150])
    expect_true(all(abs(cumsum(o$num.mark)[1:2]-c(60, 90))<=2))
    expect_equal(o$seg.mean[2], mean(y[(o$num.mark[1]+1):(o$num.mark[1]+o$num.mark[2])]))
    ## CBS and PCF still go through the same function
    expect_equal(colnames(segmentTrack(y, chr="3", starts=starts)$output), colnames(o))
    expect_equal(colnames(segmentTrack(y, chr="3", starts=starts, segmentation="pcf")$output), colnames(o))
    expect_error(segmentTrack(y, chr="3", starts=starts, segmentation="hmm"))
})

test_that("the QuickSeg penalty is in noise units: same breakpoints after scaling and shifting the track", {
    set.seed(3)
    y <- c(rnorm(80), rnorm(15, -3), rnorm(80), rnorm(40, 2))
    starts <- seq_along(y)*1e5
    o1 <- segmentTrack_quickseg(y, chr="1", starts=starts)$output
    expect_equal(nrow(o1), 4)
    o2 <- segmentTrack_quickseg(3*y+10, chr="1", starts=starts)$output
    expect_equal(o2$num.mark, o1$num.mark)
    ## a higher penalty gives no more segments; sigma=1 makes it a raw penalty
    expect_lte(nrow(segmentTrack_quickseg(y, chr="1", starts=starts, penalty=50)$output), nrow(o1))
    expect_equal(segmentTrack_quickseg(y, chr="1", starts=starts, penalty=3)$output$num.mark,
                 segmentTrack_quickseg(y, chr="1", starts=starts, penalty=3*getTrackNoise(y), sigma=1)$output$num.mark)
})

test_that("segments shorter than min.width are merged into the neighbour with the closest median", {
    set.seed(8)
    y <- c(rnorm(100), 15)                                       # outlier bin at the end of the chromosome
    starts <- seq_along(y)*1e5
    expect_equal(segmentTrack_quickseg(y, chr="1", starts=starts, min.width=1)$output$num.mark, c(100, 1))
    expect_equal(segmentTrack_quickseg(y, chr="1", starts=starts)$output$num.mark, 101)
    z <- c(rnorm(50), rnorm(3, 12), rnorm(50, 2))                # 3-bin segment, closer to the right neighbour
    expect_equal(segmentTrack_quickseg(z, chr="1", starts=seq_along(z)*1e5)$output$num.mark, c(50, 3, 50))
    expect_equal(segmentTrack_quickseg(z, chr="1", starts=seq_along(z)*1e5, min.width=5)$output$num.mark, c(50, 53))
    ## segmentTrack passes min.width only when it is given (its default 5 is for CBS and PCF)
    expect_equal(segmentTrack(z, chr="1", starts=seq_along(z)*1e5, segmentation="quickseg")$output$num.mark, c(50, 3, 50))
    expect_equal(segmentTrack(z, chr="1", starts=seq_along(z)*1e5, segmentation="quickseg", min.width=5)$output$num.mark,
                 c(50, 53))
    expect_equal(unname(.mergeShort(cbind(c(1, 3, 6), c(2, 5, 10)), c(0, 0, 5, 5, 5, 9, 9, 9, 9, 9), 3)),
                 cbind(c(1, 6), c(5, 10)))
    expect_equal(nrow(.mergeShort(cbind(1, 4), 1:4, 5)), 1)    # a chromosome shorter than min.width
})

test_that("getTrackNoise is robust to breakpoints and runs of empty bins", {
    set.seed(4)
    y <- rnorm(5000, sd=0.3)
    y[2001:3000] <- y[2001:3000]+2
    expect_equal(getTrackNoise(y), 0.3, tolerance=0.1)
    expect_equal(getTrackNoise(list(y[1:2500], y[2501:5000])), 0.3, tolerance=0.1)
    expect_equal(getTrackNoise(c(rep(0, 4000), y[1:1000])), 0.3, tolerance=0.15)
    expect_equal(getTrackNoise(rep(0, 10)), 1)
    expect_equal(getTrackNoise(3), 1)
})

test_that("segmentTrack_quickseg skips non-finite values as CBS does", {
    set.seed(5)
    y <- c(rnorm(50), rnorm(50, 2))
    y[c(10, 75)] <- c(NA, Inf)
    starts <- seq_along(y)*1e5
    o <- segmentTrack_quickseg(y, chr="1", starts=starts)$output
    expect_equal(sum(o$num.mark), 98)
    expect_false(any(is.na(o$seg.mean)))
    empty <- segmentTrack_quickseg(rep(NA_real_, 5), chr="1", starts=1:5)$output
    expect_equal(nrow(empty), 1)
    expect_equal(empty$num.mark, 0)
})

test_that("the built-in QuickSeg gives the same segments as quickseg::segment", {
    skip_if_not_installed("quickseg")
    set.seed(6)
    for(k in 1:50)
    {
        n <- sample(1:400, 1)
        y <- rnorm(n)+rep(c(0, sample(c(-1.5, 1.5), 1)), c(n%/%2, n-n%/%2))
        if(k%%2==0) y <- round(y, 1)                                  # ties between candidate values
        p <- runif(1, 0, 8)
        cpp <- quickseg::segment(y, p)
        expect_equal(unname(.quickseg(y, p)), unname(cbind(cpp$Startidx, cpp$Endidx)+1L))
    }
})

## ---------- run_sc_sequencing and fitFromLogR on simulated cells (real hg38 chr21/22 bins, no BAMs) ----------
sim_res <- function(ncell=6, seed=7)
{
    data("lSe_filtered_30000.hg38", package="ASCAT.sc", envir=environment())
    set.seed(seed)
    allchr <- c("21", "22"); binsize <- 1e6; window <- ceiling(binsize/30000)
    win <- lapply(paste0("chr", allchr), function(chr) lSe.hg38.filtered[[chr]])
    cells <- paste0("BC", seq_len(ncell))
    allTracks <- lapply(cells, function(bc)
    {
        lCTS <- lapply(seq_along(win), function(j)
        {
            se <- win[[j]]; n <- length(se$starts)
            cn <- rep(2, n)
            if(j==2) cn[round(n/2):n] <- 3                          # gain of distal chr22
            N <- rpois(n, 12*cn/2)
            len <- round(rlnorm(sum(N), log(400), 0.9))
            f <- factor(rep.int(seq_len(n), N), levels=seq_len(n))
            data.frame(space=allchr[j], start=se$starts, end=se$ends, width=se$ends-se$starts, file=bc,
                       records=N, nucleotides=vapply(split(len, f), sum, numeric(1)))
        })
        names(lCTS) <- allchr
        list(lCTS.tumour=lCTS, nlCTS.tumour=treatTrack(lCTS, window=window))
    })
    names(allTracks) <- cells
    list(stub=list(allTracks=allTracks, allTracks.processed=list(stale=TRUE), allSolutions=list(), allProfiles=list(),
                   purs=1, ploidies=seq(1.7, 3, 0.1), maxtumourpsi=5, chr=allchr, sex="female",
                   segmentation_alpha=0.01, binsize=binsize),
         allchr=allchr, binsize=binsize, cells=cells)
}

test_that("run_sc_sequencing segments with the chosen method and re-segments when the setting changes", {
    skip_on_cran()
    s <- sim_res()
    run <- function(res, ...) suppressWarnings(run_sc_sequencing(tumour_bams="sim.bam", res=res, allchr=s$allchr,
        sex="female", build="hg38", binsize=s$binsize, barcodes_10x=s$cells, purs=1, ploidies=seq(1.7, 3, 0.1),
        predict_refit=FALSE, print_results=FALSE, outdir=tempdir(), ...))
    nsegs <- function(r) sapply(r$allTracks.processed, function(tr) sum(sapply(tr$lSegs, function(x) nrow(x$output))))
    qs <- run(s$stub, segmentation="quickseg")
    expect_equal(qs$segmentation, "quickseg")
    expect_null(qs$segmentation_penalty)
    expect_equal(names(qs$allTracks.processed), s$cells)
    expect_true(all(nsegs(qs)>=3))                                   # chr21, and chr22 split at the gain
    expect_true(all(sapply(qs$allSolutions, function(x) !inherits(x, "try-error"))))
    ## same settings: segments reused; other method or penalty: segmented again
    expect_identical(run(qs, segmentation="quickseg")$allTracks.processed, qs$allTracks.processed)
    hi <- run(qs, segmentation="quickseg", segmentation_penalty=1e6)
    expect_equal(hi$segmentation_penalty, 1e6)
    expect_true(all(nsegs(hi)==2))                                   # one segment per chromosome
    pcf <- run(qs, segmentation="pcf")
    expect_equal(pcf$segmentation, "pcf")
    expect_true(all(nsegs(pcf)>=3))
    ## results saved before the option existed count as CBS
    old <- qs; old$segmentation <- NULL; old$segmentation_penalty <- NULL
    expect_identical(run(old)$allTracks.processed, qs$allTracks.processed)
    expect_error(run(s$stub, segmentation="quickseg", multipcf=TRUE), "multipcf")
})

test_that("fitFromLogR takes the segmentation from res unless given", {
    skip_on_cran()
    s <- sim_res(ncell=3)
    res <- suppressWarnings(run_sc_sequencing(tumour_bams="sim.bam", res=s$stub, allchr=s$allchr, sex="female",
        build="hg38", binsize=s$binsize, barcodes_10x=s$cells, purs=1, ploidies=seq(1.7, 3, 0.1),
        predict_refit=FALSE, print_results=FALSE, outdir=tempdir(), segmentation="quickseg", segmentation_penalty=3))
    logr <- getProcessedLogR(res)
    r1 <- fitFromLogR(res, logr)
    expect_equal(r1$segmentation, "quickseg")
    expect_equal(r1$segmentation_penalty, 3)
    expect_equal(r1$allTracks.processed$BC1$lSegs, res$allTracks.processed$BC1$lSegs)
    r2 <- fitFromLogR(res, logr, segmentation="cbs")
    expect_equal(r2$segmentation, "cbs")
    expect_null(r2$segmentation_penalty)
})
