## Simulation: CBS vs PCF vs QuickSeg as ASCAT.sc runs them (segmentTrack with segmentation=...), on
## single-cell-like logR tracks (unit noise sd). Used to pick the per-cell defaults against CBS at the default
## segmentation_alpha=0.01: QuickSeg penalty 7 (in units of the noise sd) and PCF gamma 20 (gamma=1/alpha=100,
## the multipcf translation, finds almost nothing in a single cell). Results: sim_segmentation_methods_output.txt
##   null: chromosomes without CNA -> false breakpoints per chromosome
##   alt : one central CNA of L bins shifted by delta noise sds -> share of the 2 true breakpoints found
##         within +-2 bins, and false breakpoints per chromosome
## Noise: gaussian, or gaussian with 2% of bins shifted by +-6 sd (outliers).
## Usage: Rscript diagnostics/sim_segmentation_methods.R [path/to/ASCAT.sc]
args <- commandArgs(TRUE)
repo <- if(length(args)) args[1] else "."
suppressPackageStartupMessages({ library(DNAcopy); library(copynumber) })
for(f in c("segmentTrack.R", "segmentTrack_pcf.R", "segmentTrack_quickseg.R")) source(file.path(repo, "R", f))
load(file.path(repo, "data", "SBDRYs_precomputed.rda"))
set.seed(20261028)
cat(sprintf("QuickSeg engine: %s\n", if(requireNamespace("quickseg", quietly=TRUE)) "quickseg::segment (C++)" else "built-in R"))

pcf_gamma <- c(100, 40, 20, 10)   # 100 = 1/segmentation_alpha, the default translation (also used by multipcf)
qs_penalty <- c(4, 5, 6, 7, 8)
methods <- c(list(cbs=list("cbs", NULL)),
             setNames(lapply(pcf_gamma, function(g) list("pcf", g)), paste0("pcf", pcf_gamma)),
             setNames(lapply(qs_penalty, function(g) list("quickseg", g)), paste0("qs", qs_penalty)))
## breakpoints = index of the last bin before each break. 100 kb spacing keeps the track inside the p-arm of
## chr1, as copynumber::pcf always breaks at the centromere.
breaks <- function(y, m)
{
    out <- segmentTrack(y, chr="1", starts=seq_along(y)*1e5, SBDRY=SBDRYs[["0.01"]], ALPHA=0.01,
                        segmentation=m[[1]], penalty=m[[2]])$output
    cumsum(out$num.mark)[-nrow(out)]
}
noisy <- function(n, outliers)
{
    y <- rnorm(n)
    if(outliers) { k <- sample.int(n, round(0.02*n)); y[k] <- y[k]+sample(c(-6, 6), length(k), TRUE) }
    y
}
timing <- setNames(numeric(length(methods)), names(methods))
run <- function(y, truth)
{
    sapply(names(methods), function(mm)
    {
        t0 <- proc.time()[[3]]
        b <- breaks(y, methods[[mm]])
        timing[mm] <<- timing[mm]+proc.time()[[3]]-t0
        near <- function(x, to) if(length(to)==0) rep(FALSE, length(x)) else sapply(x, function(v) any(abs(v-to)<=2))
        c(found=if(length(truth)) mean(near(truth, b)) else NA, false=sum(!near(b, truth)))
    })
}

## ---------- null: false breakpoints per chromosome ----------
null <- lapply(c(FALSE, TRUE), function(outl) lapply(c(50, 150, 300), function(n)
{
    r <- replicate(200, run(noisy(n, outl), integer(0))["false", ])
    lab <- data.frame(noise=if(outl) "outliers" else "gauss", n=n)
    list(mean=cbind(lab, t(round(rowMeans(r), 3))), any=cbind(lab, t(round(rowMeans(r>0), 3))))
}))
null <- unlist(null, recursive=FALSE)
cat("\n## NULL - mean false breakpoints per chromosome (200 chromosomes per row)\n")
print(do.call(rbind, lapply(null, `[[`, "mean")), row.names=FALSE)
cat("\n## NULL - share of chromosomes with >=1 false breakpoint\n")
print(do.call(rbind, lapply(null, `[[`, "any")), row.names=FALSE)

## ---------- alt: one CNA of L bins, shift delta ----------
alt <- unlist(lapply(c(FALSE, TRUE), function(outl) unlist(lapply(c(5, 10, 20, 40), function(L)
    lapply(c(1, 1.5, 2, 3), function(delta)
    {
        truth <- c(100, 100+L)
        r <- replicate(100, { y <- noisy(200, outl); y[101:(100+L)] <- y[101:(100+L)]+delta; run(y, truth) }, simplify=FALSE)
        lab <- data.frame(noise=if(outl) "outliers" else "gauss", L=L, delta=delta)
        list(found=cbind(lab, t(round(rowMeans(sapply(r, function(x) x["found", ])), 2))),
             false=cbind(lab, t(round(rowMeans(sapply(r, function(x) x["false", ])), 2))))
    })), recursive=FALSE)), recursive=FALSE)
cat("\n## ALT - share of the 2 true breakpoints found within +-2 bins (n=200, 100 chromosomes per row)\n")
print(do.call(rbind, lapply(alt, `[[`, "found")), row.names=FALSE)
cat("\n## ALT - breakpoints not within +-2 bins of a true one, per chromosome (includes misplaced true ones)\n")
print(do.call(rbind, lapply(alt, `[[`, "false")), row.names=FALSE)

cat("\n## total segmentation time per method (s)\n")
print(round(timing, 2))
