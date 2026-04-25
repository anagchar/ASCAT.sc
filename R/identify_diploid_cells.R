identify_diploid_cells <- function(res,
                                   ploidy_range = c(1.9, 2.1),
                                   exclude_ambiguous = TRUE,
                                   min_diploid_fraction = 0.95,
                                   diploid_cn = 2)
{
    allSols <- res$allSolutions
    allProfiles <- res$allProfiles
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

    ## Compute per-cell fraction of segments (weighted by number of bins) that are diploid.
    ## For males: chrX/Y expected CN = 1, all autosomes = diploid_cn.
    ## A cell passes if >= min_diploid_fraction of its bins are at the expected diploid CN.
    diploid_fractions <- sapply(seq_along(valid_cells), function(i) {
        cell <- valid_cells[i]
        prof <- allProfiles[[cell]]
        if(is.null(prof) || inherits(prof, "try-error") || !is.data.frame(prof))
            return(NA_real_)
        if(!"total_copy_number" %in% colnames(prof))
            return(NA_real_)
        cn <- prof$total_copy_number
        if(any(is.na(cn))) return(NA_real_)

        cell_idx <- match(cell, names(allSols))
        this_sex <- if(!is.null(cell_sex) && length(cell_sex) >= cell_idx) cell_sex[cell_idx] else NA
        chr_clean <- gsub("^chr", "", as.character(prof$chromosome))
        expected <- ifelse(!is.na(this_sex) && this_sex == "male" & chr_clean %in% c("X", "Y"),
                           1L, diploid_cn)

        ## Weight each segment by its number of bins if available, else count segments equally
        if("n_bins" %in% colnames(prof)) {
            weights <- prof$n_bins
        } else if(all(c("start", "end") %in% colnames(prof))) {
            weights <- prof$end - prof$start + 1
        } else {
            weights <- rep(1, nrow(prof))
        }
        sum(weights[cn == expected]) / sum(weights)
    })

    pass_segments <- !is.na(diploid_fractions) & diploid_fractions >= min_diploid_fraction

    is_diploid <- pass_ploidy & pass_ambig & pass_segments
    diploid_cells <- valid_cells[is_diploid]

    cat("Diploid cell identification:\n")
    cat("  Total cells:                    ", length(allSols), "\n")
    cat("  Valid solutions:                ", length(valid_cells), "\n")
    cat("  Pass ploidy [", ploidy_range[1], "-", ploidy_range[2], "]:    ",
        sum(pass_ploidy), "\n", sep="")
    if(exclude_ambiguous)
        cat("  Pass non-ambiguous:             ", sum(pass_ambig), "\n")
    cat("  Pass diploid fraction >=", min_diploid_fraction, ":",
        sum(pass_segments, na.rm=TRUE), "\n")
    if(sum(pass_segments, na.rm=TRUE) > 0) {
        frac_pass <- diploid_fractions[pass_segments & !is.na(pass_segments)]
        cat("  Diploid fraction range (passing):", round(min(frac_pass),3),
            "-", round(max(frac_pass),3), "\n")
    }
    cat("  Final diploid reference cells:  ", length(diploid_cells), "\n")

    if(length(diploid_cells) == 0)
        warning("No diploid cells found. Consider relaxing ploidy_range or lowering min_diploid_fraction.")

    diploid_cells
}
