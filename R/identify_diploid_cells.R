identify_diploid_cells <- function(res,
                                   ploidy_range = c(1.9, 2.1),
                                   exclude_ambiguous = TRUE,
                                   require_all_diploid_segments = TRUE,
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

    pass_segments <- rep(TRUE, length(valid_cells))
    if(require_all_diploid_segments) {
        pass_segments <- sapply(seq_along(valid_cells), function(i) {
            cell <- valid_cells[i]
            prof <- allProfiles[[cell]]
            if(is.null(prof) || inherits(prof, "try-error") || !is.data.frame(prof))
                return(FALSE)
            if(!"total_copy_number" %in% colnames(prof))
                return(FALSE)
            ## per-cell expected CN: chrX/Y = 1 in males, otherwise diploid_cn
            cell_idx <- match(cell, names(allSols))
            this_sex <- if(!is.null(cell_sex) && length(cell_sex) >= cell_idx) cell_sex[cell_idx] else NA
            chr_clean <- gsub("^chr", "", as.character(prof$chromosome))
            expected <- ifelse(!is.na(this_sex) && this_sex == "male" & chr_clean %in% c("X", "Y"),
                               1, diploid_cn)
            cn <- prof$total_copy_number
            !any(is.na(cn)) && all(cn == expected)
        })
    }

    is_diploid <- pass_ploidy & pass_ambig & pass_segments
    diploid_cells <- valid_cells[is_diploid]

    cat("Diploid cell identification:\n")
    cat("  Total cells:                    ", length(allSols), "\n")
    cat("  Valid solutions:                ", length(valid_cells), "\n")
    cat("  Pass ploidy [", ploidy_range[1], "-", ploidy_range[2], "]:    ",
        sum(pass_ploidy), "\n", sep="")
    if(exclude_ambiguous)
        cat("  Pass non-ambiguous:             ", sum(pass_ambig), "\n")
    if(require_all_diploid_segments)
        cat("  Pass all-segments-diploid:      ", sum(pass_segments), "\n")
    cat("  Final diploid reference cells:  ", length(diploid_cells), "\n")

    if(length(diploid_cells) == 0)
        warning("No diploid cells found. Consider relaxing ploidy_range or set require_all_diploid_segments=FALSE.")

    diploid_cells
}
