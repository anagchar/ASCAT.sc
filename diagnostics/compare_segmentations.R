## Compare per-cell segmentations (CBS, PCF, QuickSeg) of the same cells, starting from a saved
## run_sc_sequencing() result. Every method segments the same logR (res$allTracks.processed, via
## getProcessedLogR) and is fitted with fitFromLogR(), so the differences come from the segmentation only.
## BAMs are not read again.
##
## Usage (cluster):
##   Rscript compare_segmentations.R res.rds outdir [key=value ...]
## keys (comma-separated values give one run per value, "none" leaves the method out):
##   cbs=0.01         segmentation_alpha (default: the one stored in res)
##   pcf=20           gamma, in units of the cell's noise variance
##   quickseg=7       penalty, in units of the cell's noise sd
##   cells=500        random subset of 500 cells, or a file with one cell name per line (default: all cells)
##   tol=2            breakpoints of two methods at most tol bins apart count as the same breakpoint
##   examples=8       cells drawn in the pdf: the ones where the methods disagree most, plus 2 random ones
##   save_res=TRUE    also save the fitted res of each run as res_<run>.rds (default FALSE)
##   src=path/to/ASCAT.sc/R   source segmentTrack*.R and customLogR.R from a checkout, for an installed
##                    ASCAT.sc that does not have the segmentation option yet
## e.g. Rscript compare_segmentations.R res.rds cmp quickseg=5,7,10 pcf=none src=$HOME/ASCAT.sc/R
## Cores: SLURM_CPUS_PER_TASK.
##
## Outputs in outdir:
##   segmentation_per_cell.tsv   run x cell: noise sd, segments, breakpoints, cn_changes (breakpoints where the
##                               copy number changes; PCF also breaks at every centromere), purity, ploidy,
##                               frac_altered (bins off the modal copy number), rmse_sigma (rmse of logR around the
##                               segment means / noise sd: ~1 the segments explain the signal, >1 changes are
##                               missed), time_segmentation (s), failed
##   segmentation_summary.tsv    medians per run
##   segmentation_pairwise.tsv   pair of runs (a, b) x cell: share of a's breakpoints found by b (precision, b as
##                               truth), of b's found by a (recall), F1, F1 of the copy number changes only
##                               (cnbp_f1), share of bins with the same copy number, ploidy_b - ploidy_a
##   segmentation_comparison.pdf distributions, agreement with the first run, example cells

suppressPackageStartupMessages(library(ASCAT.sc))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2) stop("usage: Rscript compare_segmentations.R res.rds outdir [key=value ...]")
res <- readRDS(args[1])
outdir <- args[2]
kv <- strsplit(args[-(1:2)], "=", fixed=TRUE)
opt <- setNames(lapply(kv, function(x) paste(x[-1], collapse="=")), sapply(kv, `[`, 1))
keys <- c("cbs","pcf","quickseg","cells","tol","examples","save_res","src")
if(any(!names(opt)%in%keys))
    stop(paste0("unknown option(s): ", paste(setdiff(names(opt), keys), collapse=", "), " - use ", paste(keys, collapse=", ")))
getOpt <- function(key, default) if(is.null(opt[[key]])) default else opt[[key]]
if(!is.null(opt[["src"]]))
    for(f in c("segmentTrack.R", "segmentTrack_pcf.R", "segmentTrack_quickseg.R", "customLogR.R"))
        source(file.path(path.expand(opt[["src"]]), f))
if(!exists("fitFromLogR") || !"segmentation"%in%names(formals(fitFromLogR)))
    stop("this ASCAT.sc has no segmentation option: install the updated package or add src=path/to/ASCAT.sc/R")
MC.CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
set.seed(1)
print(paste0("## QuickSeg engine: ", if(requireNamespace("quickseg", quietly=TRUE)) "quickseg::segment (C++)" else "built-in R"))

## 1. runs: one per method and parameter value
values <- function(key, default)
{
    v <- getOpt(key, default)
    if(identical(v, "none")) return(numeric(0))
    x <- suppressWarnings(as.numeric(strsplit(as.character(v), ",")[[1]]))
    if(any(is.na(x))) stop(paste0(key, "= takes comma-separated numbers or none, not ", v))
    x
}
alpha <- if(is.null(res$segmentation_alpha)) 0.01 else res$segmentation_alpha
runs <- c(lapply(values("cbs", alpha), function(a) list(segmentation="cbs", alpha=a, penalty=NULL)),
          lapply(values("pcf", 20), function(g) list(segmentation="pcf", alpha=alpha, penalty=g)),
          lapply(values("quickseg", 7), function(p) list(segmentation="quickseg", alpha=alpha, penalty=p)))
names(runs) <- sapply(runs, function(r) paste0(r$segmentation, "_", if(r$segmentation=="cbs") r$alpha else r$penalty))
if(length(runs)<2) stop("need at least two runs to compare")
print(paste0("## runs: ", paste(names(runs), collapse=", ")))
tol <- as.numeric(getOpt("tol", 2))
save_res <- as.logical(getOpt("save_res", FALSE))

## 2. logR (bins x cells) of the cells to compare
logr <- getProcessedLogR(res)
cells <- colnames(logr)
if(!is.null(opt[["cells"]]))
{
    if(file.exists(opt[["cells"]]))
        cells <- intersect(readLines(opt[["cells"]]), cells)
    else if(!is.na(suppressWarnings(as.integer(opt[["cells"]]))))
        cells <- sample(cells, min(as.integer(opt[["cells"]]), length(cells)))
    else
        stop("cells= takes a number of cells or a file with one cell name per line")
}
ok <- apply(is.finite(logr[, cells, drop=FALSE]), 2, all)
if(any(!ok)) print(paste0("## leaving out ", sum(!ok), " cells with NA/Inf logR"))
cells <- cells[ok]
if(length(cells)==0) stop("no cells left to compare (check cells=)")
nb <- sapply(res$allTracks.processed[[cells[1]]]$lCTS, nrow)
chr <- rep(seq_along(nb), nb)
chrnames <- if(length(res$chr)==length(nb)) res$chr else names(nb)
sigma <- vapply(cells, function(cc) getTrackNoise(split(logr[, cc], chr)), numeric(1), USE.NAMES=FALSE)
print(paste0("## ", length(cells), " cells, ", nrow(logr), " bins; median noise sd ", signif(median(sigma), 3)))

## segments of a fitted profile (getProfile): bins per segment, copy number, breakpoints (index of the last bin
## before each new segment within a chromosome) and the breakpoints where the copy number changes (PCF also
## breaks at every centromere, which usually leaves the copy number as it is)
segsOf <- function(prof)
{
    if(!is.data.frame(prof) || sum(prof$num.mark)!=nrow(logr)) return(NULL)
    k <- which(duplicated(prof$chromosome, fromLast=TRUE))
    bp <- cumsum(prof$num.mark)[k]
    change <- prof$total_copy_number[k]!=prof$total_copy_number[k+1]
    list(nm=prof$num.mark,
         cn=prof$total_copy_number,
         bp=bp,
         cnbp=bp[!is.na(change) & change])
}
modalCN <- function(cn) as.numeric(names(which.max(table(cn))))
cellStats <- function(s, y, sd)
{
    if(is.null(s)) return(c(n_segments=NA, n_breakpoints=NA, n_cn_changes=NA, frac_altered=NA, rmse_sigma=NA))
    cn <- rep(s$cn, s$nm)
    fitted <- ave(y, rep(seq_along(s$nm), s$nm))
    c(n_segments=length(s$nm),
      n_breakpoints=length(s$bp),
      n_cn_changes=length(s$cnbp),
      frac_altered=mean(cn!=modalCN(cn), na.rm=TRUE),
      rmse_sigma=sqrt(mean((y-fitted)^2))/sd)
}

## 3. segment and fit every run (slow part); purs/ploidies/maxtumourpsi/isPON as in res
mysafe <- function(x, f) if(is.list(x) && !is.null(x[[f]])) x[[f]] else NA
segs <- list()
percell <- NULL
for(m in names(runs))
{
    print(paste0("## ", m))
    r <- fitFromLogR(res, logr, cells=cells, segmentation=runs[[m]]$segmentation,
                     segmentation_penalty=runs[[m]]$penalty, segmentation_alpha=runs[[m]]$alpha,
                     MC.CORES=MC.CORES)
    if(save_res) saveRDS(r, file.path(outdir, paste0("res_", m, ".rds")))
    segs[[m]] <- lapply(r$allProfiles[cells], segsOf)
    st <- t(sapply(seq_along(cells), function(i) cellStats(segs[[m]][[i]], logr[, cells[i]], sigma[i])))
    percell <- rbind(percell, data.frame(run=m,
                                         cell=cells,
                                         sigma=sigma,
                                         st,
                                         purity=sapply(r$allSolutions[cells], mysafe, "purity"),
                                         ploidy=sapply(r$allSolutions[cells], mysafe, "ploidy"),
                                         time_segmentation=if(is.null(r$timetosegment)) NA else r$timetosegment[cells],
                                         failed=sapply(segs[[m]], is.null),
                                         row.names=NULL,
                                         stringsAsFactors=FALSE))
    rm(r); gc()
}
percell$run <- factor(percell$run, levels=names(runs))
write.table(percell, file.path(outdir, "segmentation_per_cell.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
summ <- do.call(rbind, lapply(split(percell, percell$run), function(d)
    data.frame(run=d$run[1],
               cells=nrow(d),
               failed=sum(d$failed),
               breakpoints=median(d$n_breakpoints, na.rm=TRUE),
               cn_changes=median(d$n_cn_changes, na.rm=TRUE),
               rmse_sigma=round(median(d$rmse_sigma, na.rm=TRUE), 3),
               frac_altered=round(median(d$frac_altered, na.rm=TRUE), 3),
               ploidy=round(median(d$ploidy, na.rm=TRUE), 2),
               seconds_per_cell=signif(median(d$time_segmentation, na.rm=TRUE), 3),
               total_seconds=round(sum(d$time_segmentation, na.rm=TRUE), 1))))
write.table(summ, file.path(outdir, "segmentation_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
print("## medians per run (seconds: segmentation only, summed over cells for total_seconds)")
print(summ, row.names=FALSE)

## 4. agreement between every pair of runs, per cell
near <- function(a, b) vapply(a, function(x) any(abs(x-b)<=tol), logical(1))
f1 <- function(a, b)
{
    p <- if(length(a)) mean(near(a, b)) else NA
    r <- if(length(b)) mean(near(b, a)) else NA
    c(precision=p, recall=r, f1=if(!length(a) && !length(b)) 1 else if(!length(a) || !length(b)) 0 else
        if(p+r==0) 0 else 2*p*r/(p+r))
}
pairStats <- function(sa, sb)
{
    if(is.null(sa) || is.null(sb))
        return(c(bp_precision=NA, bp_recall=NA, bp_f1=NA, cnbp_f1=NA, cn_agreement=NA))
    b <- f1(sa$bp, sb$bp)
    c(bp_precision=b[["precision"]],
      bp_recall=b[["recall"]],
      bp_f1=b[["f1"]],
      cnbp_f1=f1(sa$cnbp, sb$cnbp)[["f1"]],
      cn_agreement=mean(rep(sa$cn, sa$nm)==rep(sb$cn, sb$nm)))
}
pairs <- combn(names(runs), 2, simplify=FALSE)
pairwise <- do.call(rbind, lapply(pairs, function(pr)
{
    st <- t(sapply(seq_along(cells), function(i) pairStats(segs[[pr[1]]][[i]], segs[[pr[2]]][[i]])))
    pa <- percell$ploidy[percell$run==pr[1]]
    pb <- percell$ploidy[percell$run==pr[2]]
    data.frame(run_a=pr[1], run_b=pr[2], cell=cells, st, ploidy_diff=pb-pa, row.names=NULL, stringsAsFactors=FALSE)
}))
write.table(pairwise, file.path(outdir, "segmentation_pairwise.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
psumm <- do.call(rbind, lapply(split(pairwise, paste(pairwise$run_a, pairwise$run_b)), function(d)
    data.frame(run_a=d$run_a[1],
               run_b=d$run_b[1],
               bp_f1=round(median(d$bp_f1, na.rm=TRUE), 3),
               cnbp_f1=round(median(d$cnbp_f1, na.rm=TRUE), 3),
               cn_agreement=round(median(d$cn_agreement, na.rm=TRUE), 3),
               cells_ploidy_diff_gt_0.5=round(mean(abs(d$ploidy_diff)>0.5, na.rm=TRUE), 3))))
print(paste0("## agreement between runs (medians over cells; breakpoints within ", tol, " bins)"))
print(psumm, row.names=FALSE)

## 5. plots
cols <- setNames(hcl.colors(length(runs), "Dark 3"), names(runs))
pdf(file.path(outdir, "segmentation_comparison.pdf"), width=11, height=8.5)
par(mfrow=c(2, 3), mar=c(8, 4, 3, 1))
boxplot(n_breakpoints ~ run, data=percell, las=2, col=cols, xlab="", main="breakpoints per cell", outline=FALSE)
boxplot(n_cn_changes ~ run, data=percell, las=2, col=cols, xlab="", main="copy number changes per cell", outline=FALSE)
boxplot(rmse_sigma ~ run, data=percell, las=2, col=cols, xlab="", main="rmse around segments / noise sd", outline=FALSE)
abline(h=1, lty=2)
boxplot(ploidy ~ run, data=percell, las=2, col=cols, xlab="", main="ploidy")
if(any(is.finite(percell$time_segmentation)))
    boxplot(pmax(time_segmentation, 1e-3) ~ run, data=percell, las=2, col=cols, xlab="", ylab="", log="y",
            main="segmentation time per cell (s, log scale)")
## agreement of every other run with the first one
ref <- names(runs)[1]
par(mfrow=c(min(3, length(runs)-1), 4), mar=c(4.5, 4.5, 3, 1))
for(m in names(runs)[-1])
{
    d <- pairwise[pairwise$run_a==ref & pairwise$run_b==m, ]
    plot(jitter(percell$ploidy[percell$run==ref], amount=0.02), jitter(percell$ploidy[percell$run==m], amount=0.02),
         pch=16, cex=0.5, col=cols[m], xlab=paste("ploidy", ref), ylab=paste("ploidy", m), main=paste(m, "vs", ref))
    abline(0, 1, lty=2)
    hist(d$cn_agreement, breaks=seq(0, 1, 0.05), col=cols[m], xlab="share of bins with the same copy number",
         main=sprintf("CN agreement (median %.2f)", median(d$cn_agreement, na.rm=TRUE)))
    hist(d$bp_f1, breaks=seq(0, 1, 0.05), col=cols[m], xlab=paste0("F1, all breakpoints (within ", tol, " bins)"),
         main=sprintf("breakpoints (median %.2f)", median(d$bp_f1, na.rm=TRUE)))
    hist(d$cnbp_f1, breaks=seq(0, 1, 0.05), col=cols[m], xlab=paste0("F1, CN changes (within ", tol, " bins)"),
         main=sprintf("CN changes (median %.2f)", median(d$cnbp_f1, na.rm=TRUE)))
}
## example cells: lowest copy number agreement over all pairs, plus 2 random cells
nex <- min(as.integer(getOpt("examples", 8)), length(cells))
worst <- tapply(pairwise$cn_agreement, factor(pairwise$cell, levels=cells), function(x) suppressWarnings(min(x, na.rm=TRUE)))
ex <- head(names(sort(worst[is.finite(worst)])), max(nex-2, 0))
ex <- c(ex, head(sample(setdiff(cells, ex)), nex-length(ex)))
plotCell <- function(cc)
{
    y <- logr[, cc]
    yl <- quantile(y, c(0.002, 0.998))
    yl <- yl+c(-1, 1)*0.1*diff(yl)
    par(mfrow=c(length(runs), 1), mar=c(1, 4, 1.6, 1), oma=c(2, 0, 2.5, 0))
    for(m in names(runs))
    {
        plot(seq_along(y), pmin(pmax(y, yl[1]), yl[2]), pch=16, cex=0.3, col="grey70", ylim=yl, xaxs="i", xaxt="n",
             xlab="", ylab="logR")
        abline(v=cumsum(nb)[-length(nb)]+0.5, col="grey85")
        s <- segs[[m]][[match(cc, cells)]]
        if(!is.null(s))
        {
            e <- cumsum(s$nm)
            b <- e-s$nm+1
            mu <- sapply(seq_along(e), function(i) mean(y[b[i]:e[i]]))
            mode <- modalCN(rep(s$cn, s$nm))
            segments(b-0.5, mu, e+0.5, mu, lwd=2.5, col=ifelse(s$cn>mode, "#B2182B", ifelse(s$cn<mode, "#2166AC", "black")))
        }
        row <- percell[percell$run==m & percell$cell==cc, ]
        title(main=sprintf("%s   ploidy %.2f   purity %.2f   %s breakpoints   rmse/sd %.2f", m, row$ploidy, row$purity,
                           row$n_breakpoints, row$rmse_sigma), adj=0, font.main=1, cex.main=0.9)
    }
    axis(1, at=cumsum(nb)-nb/2, labels=chrnames, tick=FALSE, cex.axis=0.7, line=-0.5)
    mtext(paste0(cc, "   (segments: red above, blue below the modal copy number)"), outer=TRUE, line=0.8, font=2)
}
for(cc in ex) plotCell(cc)
invisible(dev.off())
print(paste0("## done: ", normalizePath(outdir)))
