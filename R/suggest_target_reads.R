## suggest_target_reads() ----------------------------------------------------
##
## Inspect, from the actual pooled depth and the number of cells, which
## TARGET_READS gives equal-count (pseudobulk) bins that are (a) fine enough to
## localise CNAs, (b) deep enough per cell to beat Poisson shot noise, and
## (c) free enough of residual GC / read-length structure that a separate
## correction would not change much.
##
## Why scan instead of guess: equal-count binning already absorbs the *shared*
## technical bias into the bin boundaries, so the only knob left is how many
## reads to put in each bin. As TARGET_READS grows:
##   - per-cell shot noise (1/sqrt(reads_per_cell_per_bin))      goes DOWN  (good)
##   - residual GC/read-length structure (fraction of bin var)   goes UP    (bad)
##   - number of bins / resolution                               goes DOWN  (bad)
## The useful TARGET_READS is the smallest one whose per-cell shot noise is low
## enough for detection while GC/RL structure and bins-per-chromosome are still
## acceptable. This function measures all three on the real data.
##
## Inputs
##   res             : a run_sc_sequencing() result (needs res$allTracks with
##                     per-cell $lCTS.tumour small-window tracks, res$build, res$chr)
##   reads_per_cell  : grid of target reads-per-cell-per-bin R; TARGET_READS = R * n_cells
##   diploid_cells   : optional vector of normal/diploid barcodes to use for the
##                     residual-GC test (their logR should be flat). Defaults to a
##                     random sample of all cells.
##   n_gc_cells      : how many cells to use for the GC/RL regression (speed)
##   gc_r2_threshold : max fraction of bin-to-bin variance GC+RL may explain for
##                     "correction not needed" (default 0.02 = 2%)
##   shotnoise_cv_threshold : max acceptable per-cell shot-noise CV (default 0.25)
##   min_bins_per_chr: minimum bins on the smallest chromosome (>=3 avoids the
##                     smooth.CNA degeneracy)
##   lGCT_smallwindow: optional small-window GC list (per chr); if NULL it is
##                     reloaded from the package data for res$build (hg19/hg38)
##
## Returns (invisibly) a list: $grid (one row per candidate), $recommended_target_reads,
## $total_pooled, $n_cells, $mean_cell_depth. Also prints a readable summary.

suggest_target_reads <- function(res,
                                 reads_per_cell        = c(10, 20, 30, 50, 75, 100, 150),
                                 diploid_cells         = NULL,
                                 n_gc_cells            = 100,
                                 gc_r2_threshold       = 0.02,
                                 shotnoise_cv_threshold= 0.25,
                                 min_bins_per_chr      = 3,
                                 lGCT_smallwindow      = NULL,
                                 plot                  = FALSE)
{
    if(!exists("getPseudobulkGroups"))
        stop("getPseudobulkGroups() not found - source it (package R/) before calling.")
    if(is.null(res$allTracks) || is.null(res$allTracks[[1]]$lCTS.tumour))
        stop("res$allTracks[[cell]]$lCTS.tumour (small-window tracks) not found.")

    cells     <- names(res$allTracks)
    n_cells   <- length(cells)
    lCTS_list <- lapply(res$allTracks, function(x) x$lCTS.tumour)
    chrs      <- names(lCTS_list[[1]])
    nchr      <- length(chrs)

    ## ---- pooled small-window coverage (reads + nucleotides) -----------------
    pooled_rec <- lapply(seq_len(nchr), function(i)
        Reduce(`+`, lapply(lCTS_list, function(cl) cl[[i]]$records)))
    pooled_nuc <- lapply(seq_len(nchr), function(i)
        Reduce(`+`, lapply(lCTS_list, function(cl) cl[[i]]$nucleotides)))
    names(pooled_rec) <- names(pooled_nuc) <- chrs
    total_pooled    <- sum(sapply(pooled_rec, sum, na.rm = TRUE))
    mean_cell_depth <- total_pooled / n_cells

    ## ---- small-window GC ----------------------------------------------------
    lGCT_sw <- if(!is.null(lGCT_smallwindow)) lGCT_smallwindow
               else .load_smallwindow_gc(res$build, chrs)
    have_gc <- !is.null(lGCT_sw) &&
               all(sapply(seq_len(nchr), function(i)
                   length(lGCT_sw[[i]]) == length(pooled_rec[[i]])))
    if(!is.null(lGCT_sw) && !have_gc)
        warning("small-window GC track does not align with coverage; skipping GC test.")

    ## ---- cells used for the residual-GC test --------------------------------
    gc_cells <- if(!is.null(diploid_cells)) intersect(diploid_cells, cells) else cells
    if(length(gc_cells) > n_gc_cells) gc_cells <- sample(gc_cells, n_gc_cells)

    ## fast contiguous-group rebinning via cumulative sums
    rebin <- function(x, g) { cs <- c(0, cumsum(x)); cs[g$ends + 1] - cs[g$starts] }

    ## ---- scan the grid ------------------------------------------------------
    rows <- lapply(reads_per_cell, function(R)
    {
        target <- R * n_cells
        groups <- getPseudobulkGroups(lCTS_list, target_reads = target)

        binned_rec <- lapply(seq_len(nchr), function(i) rebin(pooled_rec[[i]], groups[[i]]))
        binned_nuc <- lapply(seq_len(nchr), function(i) rebin(pooled_nuc[[i]], groups[[i]]))
        nbins_chr  <- sapply(groups, function(g) length(g$starts))
        widths_mb  <- unlist(lapply(seq_len(nchr), function(i){
            st <- lCTS_list[[1]][[i]]$start; en <- lCTS_list[[1]][[i]]$end
            (en[groups[[i]]$ends] - st[groups[[i]]$starts]) / 1e6
        }))

        n_bins   <- sum(nbins_chr)
        reads_pb_pool <- unlist(binned_rec)
        reads_pb_cell <- median(reads_pb_pool, na.rm = TRUE) / n_cells   # ~reads per cell per bin
        shot_cv  <- 1 / sqrt(max(reads_pb_cell, 1e-9))

        ## residual GC / read-length structure on the normal cells
        gc_r2 <- NA; gcrl_r2 <- NA
        if(have_gc && length(gc_cells) > 0)
        {
            ## per-bin GC = unweighted mean of window GCs (matches treatGCT, i.e.
            ## exactly the GC the pipeline's correction would regress against)
            binGC <- unlist(lapply(seq_len(nchr), function(i){
                g <- groups[[i]]
                rebin(lGCT_sw[[i]], g) / (g$ends - g$starts + 1)
            }))
            binRL <- unlist(binned_nuc) / pmax(unlist(binned_rec), 1)
            logR  <- unlist(lapply(gc_cells, function(cell){
                cl <- lCTS_list[[cell]]
                cnt <- unlist(lapply(seq_len(nchr), function(i) rebin(cl[[i]]$records, groups[[i]])))
                log2((cnt + 1) / median(cnt + 1, na.rm = TRUE))
            }))
            GCv <- rep(binGC, length(gc_cells)); RLv <- rep(binRL, length(gc_cells))
            ok  <- is.finite(logR) & is.finite(GCv) & is.finite(RLv)
            gc_r2   <- summary(lm(logR[ok] ~ GCv[ok]))$r.squared
            gcrl_r2 <- summary(lm(logR[ok] ~ GCv[ok] + RLv[ok]))$r.squared
        }

        data.frame(reads_per_cell      = R,
                   target_reads        = target,
                   n_bins              = n_bins,
                   min_bins_per_chr    = min(nbins_chr),
                   median_bin_mb       = round(median(widths_mb, na.rm = TRUE), 3),
                   reads_per_cell_bin  = round(reads_pb_cell, 1),
                   shotnoise_cv        = round(shot_cv, 3),
                   gc_r2               = round(gc_r2, 4),
                   gc_readlen_r2       = round(gcrl_r2, 4))
    })
    grid <- do.call(rbind, rows)

    ## ---- recommendation -----------------------------------------------------
    feas <- grid$shotnoise_cv <= shotnoise_cv_threshold &
            grid$min_bins_per_chr >= min_bins_per_chr &
            (is.na(grid$gc_readlen_r2) | grid$gc_readlen_r2 <= gc_r2_threshold)
    rec <- if(any(feas)) grid$target_reads[which(feas)[1]] else NA_real_

    cat("\n=== suggest_target_reads ===\n")
    cat(sprintf("cells = %d | total pooled reads = %s | mean depth/cell = %s\n",
                n_cells, format(total_pooled, big.mark=","),
                format(round(mean_cell_depth), big.mark=",")))
    if(!have_gc) cat("(GC/read-length test skipped - no aligned small-window GC)\n")
    print(grid, row.names = FALSE)
    cat("\nReading the table: shotnoise_cv DOWN as target grows (good); ",
        "gc_readlen_r2 = fraction of bin variance a GC/RL correction would remove ",
        "(want it small => correction unnecessary); min_bins_per_chr must be >= ",
        min_bins_per_chr, ".\n", sep="")
    if(!is.na(rec))
        cat(sprintf("\n>> Recommended TARGET_READS = %s  (R=%d reads/cell/bin, ~%d bins)\n",
                    format(rec, big.mark=","),
                    grid$reads_per_cell[which(feas)[1]],
                    grid$n_bins[which(feas)[1]]))
    else
        cat("\n>> No candidate satisfies all thresholds: at the depth needed to beat\n",
            "   shot noise, GC/RL structure or resolution is already too large. This\n",
            "   data likely DOES need correction (or more depth). Relax thresholds or\n",
            "   inspect the table.\n", sep="")

    if(plot && have_gc) .plot_target_scan(grid)

    invisible(list(grid = grid,
                   recommended_target_reads = rec,
                   total_pooled = total_pooled,
                   n_cells = n_cells,
                   mean_cell_depth = mean_cell_depth))
}

## reload the small-window GC track the pipeline uses, aligned to `chrs`.
## Mirrors run_sc_sequencing: load both the GC and start/end objects and copy
## the lSe names onto the GC list before indexing by chromosome.
.load_smallwindow_gc <- function(build, chrs)
{
    cfg <- switch(build,
        hg38 = list(gc="lGCT_filtered_30000.hg38", se="lSe_filtered_30000.hg38",
                    gcobj="lGCT.hg38.filtered", seobj="lSe.hg38.filtered",
                    key=paste0("chr", gsub("chr", "", chrs))),
        hg19 = list(gc="lGCT_filtered_30000.hg19", se="lSe_filtered_30000.hg19",
                    gcobj="lGCT.hg19.filtered", seobj="lSe.hg19.filtered",
                    key=gsub("chr", "", chrs)),
        NULL)
    if(is.null(cfg)) { warning("auto GC load supports hg19/hg38 only; pass lGCT_smallwindow="); return(NULL) }
    e  <- new.env()
    ok <- tryCatch({ utils::data(list = c(cfg$gc, cfg$se), package = "ASCAT.sc", envir = e); TRUE },
                   error = function(err) FALSE)
    if(!ok) { warning("could not load ", cfg$gc, " / ", cfg$se, " - pass lGCT_smallwindow="); return(NULL) }
    lg <- get(cfg$gcobj, envir = e)
    ls <- get(cfg$seobj, envir = e)
    names(lg) <- names(ls)
    out <- lapply(cfg$key, function(k) lg[[k]]); names(out) <- chrs
    out
}

## simple trade-off plot (target_reads vs shot noise / GC-R2 / n_bins)
.plot_target_scan <- function(grid)
{
    op <- par(mar = c(4,4,2,4)); on.exit(par(op))
    x <- grid$target_reads
    plot(x, grid$shotnoise_cv, type = "b", pch = 19, log = "x",
         xlab = "TARGET_READS (pooled reads/bin)", ylab = "shot-noise CV (per cell)",
         main = "Choosing TARGET_READS")
    par(new = TRUE)
    plot(x, grid$gc_readlen_r2, type = "b", pch = 17, col = "red", log = "x",
         axes = FALSE, xlab = "", ylab = "")
    axis(4, col = "red", col.axis = "red")
    mtext("GC+readlen R^2 (correction would remove)", side = 4, line = 2.5, col = "red")
    legend("top", c("shot-noise CV", "GC+RL R^2"), pch = c(19,17),
           col = c("black","red"), bty = "n")
}
