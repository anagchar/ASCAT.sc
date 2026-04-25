identify_diploid_cells <- function(res,
                                   ploidy_range = c(1.9, 2.1),
                                   exclude_ambiguous = TRUE,
                                   max_focal_fraction = 0.03,
                                   mapd_threshold = NULL,
                                   mapd_mad_k = 3,
                                   min_bin_coverage = 0.5,
                                   diploid_cn = 2)
{
    allSols <- res$allSolutions
    allProfiles <- res$allProfiles
    allTracks <- res$allTracks.processed
    cell_sex <- res$sex

    is_valid <- sapply(allSols, function(x) !is.null(x) && is.list(x) && !inherits(x, "try-error"))
    valid_cells <- names(allSols)[is_valid]
    if(length(valid_cells) == 0) {
        warning("No valid solutions in res$allSolutions")
        return(character(0))
    }

    ploidies <- sapply(allSols[valid_cells], function(x) x$ploidy)
    pass_ploidy <- ploidies >= ploidy_range[1] & ploidies <= ploidy_range[2]

    pass_ambig <- rep(TRUE, length(valid_cells))
    if(exclude_ambiguous) {
        pass_ambig <- !sapply(allSols[valid_cells], function(x) {
            if(!is.null(x$ambiguous)) isTRUE(x$ambiguous) else FALSE
        })
    }

    ## Per-cell metrics: max_focal_fraction (segment-based) and MAPD (raw-based)
    metrics <- lapply(valid_cells, function(cell) {
        prof <- allProfiles[[cell]]
        if(is.null(prof) || inherits(prof, "try-error") || !is.data.frame(prof))
            return(list(focal = NA_real_, mapd = NA_real_, coverage = NA_real_))
        if(!"total_copy_number" %in% colnames(prof))
            return(list(focal = NA_real_, mapd = NA_real_, coverage = NA_real_))

        cn <- prof$total_copy_number
        cell_idx <- match(cell, names(allSols))
        this_sex <- if(!is.null(cell_sex) && length(cell_sex) >= cell_idx) cell_sex[cell_idx] else NA
        chr_clean <- gsub("^chr", "", as.character(prof$chromosome))
        expected <- ifelse(!is.na(this_sex) && this_sex == "male" & chr_clean %in% c("X", "Y"),
                           1L, diploid_cn)

        if("num.mark" %in% colnames(prof)) {
            weights <- as.numeric(prof$num.mark)
        } else if(all(c("start", "end") %in% colnames(prof))) {
            weights <- prof$end - prof$start + 1
        } else {
            weights <- rep(1, nrow(prof))
        }

        ## max_focal_fraction: largest contiguous run of non-diploid segments,
        ## per chromosome, weighted by num.mark, divided by total weight.
        is_nondip <- !is.na(cn) & cn != expected
        total_w <- sum(weights, na.rm=TRUE)
        if(total_w <= 0) {
            focal <- NA_real_
        } else {
            ## per-chromosome run-length
            run_weights <- numeric(0)
            for(ch in unique(chr_clean)) {
                idx <- which(chr_clean == ch)
                if(length(idx) == 0) next
                r <- rle(is_nondip[idx])
                ends <- cumsum(r$lengths)
                starts <- ends - r$lengths + 1
                for(j in seq_along(r$values)) {
                    if(isTRUE(r$values[j])) {
                        seg_idx <- idx[starts[j]:ends[j]]
                        run_weights <- c(run_weights, sum(weights[seg_idx], na.rm=TRUE))
                    }
                }
            }
            focal <- if(length(run_weights) == 0) 0 else max(run_weights) / total_w
        }

        ## MAPD from raw smoothed bins, matching the interactive QC formulation:
        ## trimmed-median baseline -> log2 ratio -> median |adjacent diff|.
        mapd <- NA_real_
        coverage <- NA_real_
        if(!is.null(allTracks) && !is.null(allTracks[[cell]]) &&
           !is.null(allTracks[[cell]]$lCTS)) {
            lcts <- allTracks[[cell]]$lCTS
            sm <- tryCatch(unlist(lapply(lcts, function(x) x$smoothed)),
                           error = function(e) NULL)
            if(!is.null(sm) && length(sm) > 0) {
                ploidy_cell <- allSols[[cell]]$ploidy
                raw_vals <- ploidy_cell * 2^sm
                coverage <- mean(!is.na(raw_vals) & raw_vals > 0)
                valid <- !is.na(raw_vals) & raw_vals > 0
                rv <- raw_vals[valid]
                if(length(rv) >= 10) {
                    sv <- sort(rv)
                    lo <- ceiling(0.05 * length(sv))
                    hi <- floor(0.95 * length(sv))
                    if(lo < 1) lo <- 1
                    if(hi > length(sv)) hi <- length(sv)
                    baseline <- median(sv[lo:hi])
                    if(!is.na(baseline) && baseline > 0) {
                        log2r <- log2(rv / baseline)
                        mapd <- median(abs(diff(log2r)), na.rm=TRUE)
                    }
                }
            }
        }
        list(focal = focal, mapd = mapd, coverage = coverage)
    })

    focal_fracs <- sapply(metrics, `[[`, "focal")
    mapds       <- sapply(metrics, `[[`, "mapd")
    coverages   <- sapply(metrics, `[[`, "coverage")

    ## Resolve MAPD threshold (data-driven if NULL).
    mapd_for_calc <- mapds[!is.na(mapds)]
    if(is.null(mapd_threshold)) {
        if(length(mapd_for_calc) >= 5) {
            mapd_threshold_used <- median(mapd_for_calc) + mapd_mad_k * mad(mapd_for_calc)
        } else {
            mapd_threshold_used <- Inf
        }
    } else {
        mapd_threshold_used <- mapd_threshold
    }

    pass_focal    <- !is.na(focal_fracs) & focal_fracs <= max_focal_fraction
    pass_mapd     <- is.na(mapds) | mapds <= mapd_threshold_used
    pass_coverage <- is.na(coverages) | coverages >= min_bin_coverage

    is_diploid <- pass_ploidy & pass_ambig & pass_focal & pass_mapd & pass_coverage
    diploid_cells <- valid_cells[is_diploid]

    cat("Diploid cell identification:\n")
    cat("  Total cells:                     ", length(allSols), "\n")
    cat("  Valid solutions:                 ", length(valid_cells), "\n")
    cat("  Pass ploidy [", ploidy_range[1], "-", ploidy_range[2], "]:    ",
        sum(pass_ploidy), "\n", sep="")
    if(exclude_ambiguous)
        cat("  Pass non-ambiguous:              ", sum(pass_ambig), "\n")
    cat("  Pass max_focal_fraction <=", max_focal_fraction, ":   ",
        sum(pass_focal, na.rm=TRUE), "\n")
    cat("  Pass MAPD <=", round(mapd_threshold_used, 3),
        if(is.null(mapd_threshold)) " (auto)" else " (user)", ":  ",
        sum(pass_mapd, na.rm=TRUE), "\n", sep="")
    cat("  Pass bin coverage >=", min_bin_coverage, ":         ",
        sum(pass_coverage, na.rm=TRUE), "\n")
    if(any(!is.na(focal_fracs[is_diploid]))) {
        cat("  Focal fraction range (kept):    ",
            round(min(focal_fracs[is_diploid], na.rm=TRUE), 4), "-",
            round(max(focal_fracs[is_diploid], na.rm=TRUE), 4), "\n")
    }
    if(any(!is.na(mapds[is_diploid]))) {
        cat("  MAPD range (kept):              ",
            round(min(mapds[is_diploid], na.rm=TRUE), 3), "-",
            round(max(mapds[is_diploid], na.rm=TRUE), 3), "\n")
    }
    cat("  Final diploid reference cells:   ", length(diploid_cells), "\n")

    if(length(diploid_cells) == 0)
        warning("No diploid cells found. Consider relaxing ploidy_range, ",
                "raising max_focal_fraction, or raising mapd_threshold.")

    attr(diploid_cells, "metrics") <- data.frame(
        cell = valid_cells,
        ploidy = ploidies,
        max_focal_fraction = focal_fracs,
        mapd = mapds,
        bin_coverage = coverages,
        passes = is_diploid,
        stringsAsFactors = FALSE
    )
    attr(diploid_cells, "mapd_threshold") <- mapd_threshold_used

    diploid_cells
}
