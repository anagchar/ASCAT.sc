## ------------------------------------------------------------------------
## Count-space PCA panel of normals for ASCAT.sc::smoothNormals
##
## Return an lNormals object:
##      lNormals[[j]][[chr]]$records
## with j = 1 ... M+1 panel members (mean profile + M PC directions).
##
## ------------------------------------------------------------------------

X <- 50L # number of cells to sample from each lCT in lCT_list
N <- 1000L # number of randomly sampled pseudobulks
M <- 50L # number of PCs to return

create_pca_pon <- function(lCT_list,
                           X = 50L,
                           N = 1000L,
                           M = 50L,
                           replace = FALSE,
                           seed = 1L,
                           verbose = TRUE,
                           plot_knee = TRUE,
                           save_knee_plot = TRUE,
                           knee_plot_file = "pca_pon_knee_plot.pdf",
                           cumulative_thresholds = c(0.8, 0.9, 0.95, 0.99)) {
    stopifnot(is.list(lCT_list), length(lCT_list) >= 2L)

    ## -- 1. flatten to bins x cells --------------------------------
    template <- lCT_list[[1]]
    bins_per_chr <- sapply(template, nrow)
    chr_names <- names(template)

    counts <- sapply(lCT_list, function(x) {
        unlist(lapply(chr_names, function(y) x[[y]]$records))
    })
    counts <- as.matrix(counts)
    storage.mode(counts) <- "double" # ensure double for PCA

    if (nrow(counts) != sum(bins_per_chr)) {
        stop("bin count mismatch between cells - do all cells share the same lSe?")
    }

    if (anyNA(counts)) {
        stop("counts contains NA")
    }

    if (!replace && X > ncol(counts)) {
        stop("X (", X, ") exceeds the number of cells (", ncol(counts), ")")
    }

    n_bins <- nrow(counts)
    n_cells <- ncol(counts)

    ## rank ceiling: resampling cannot create more directions that cells
    rank_ceiling <- min(
        n_cells - 1L,
        N - 1L,
        n_bins
    )
    if (M > rank_ceiling) {
        warning(
            "M = ", M, " exceeds the rank ceiling set by ", n_cells,
            " cells; reducing to ", rank_ceiling
        )
        M <- rank_ceiling
    }

    ## -- 2. pseudo-bulks via one sparse matrix product -------------
    set.seed(seed)
    idx <- replicate(N, sample.int(n_cells, X, replace = replace))

    S <- Matrix::sparseMatrix(
        i = as.vector(idx),
        j = rep(seq_len(N), each = X),
        x = 1,
        dims = c(n_cells, N)
    )

    # bins x N pseudobulks, containing summed counts
    pb <- as.matrix(counts %*% S) # bins x N, summed

    ## -- 3. library size normalization to the median --------------- (sweep counts between bins of different pseudobulks)
    lib <- colSums(pb)
    if (any(lib <= 0)) stop("at least one pseudobulk has zero total counts")
    pb <- sweep(pb, 2L, lib, "/") * median(lib)

    ## -- 4. PCA ------------------------------------------------------
    mu <- rowMeans(pb)
    Zc <- pb - mu

    keep <- apply(Zc, 1L, function(z) any(z != 0))
    if (sum(keep) < M + 1L) {
        stop("fewer variable bins (", sum(keep), ") than requested components (", M + 1L, ")")
    }

    sv <- svd(Zc, nu = M, nv = M)
    U <- sv$u # bins x M
    s <- sv$d[seq_len(M)] / sqrt(N - 1) # SD of scores on each PC

    ev <- sv$d^2 / sum(sv$d^2)

    ## Checking the "knee" of the variance explained curve
    pc_variance <- sv$d^2
    ## Proportion of total variants explained by every PC
    explained_variance <- pc_variance / sum(pc_variance)

    # Eigenvalues
    eigenvalues <- sv$d^2 / (N - 1)

    cumulative_variance <- cumsum(explained_variance)

    n_available_pcs <- length(explained_variance)

    ## -- 4a. plot the variance explained curve --------------------------
    pc_number <- seq_len(n_available_pcs)

    if (n_available_pcs >= 3L) {
        x_norm <- (
            pc_number - min(pc_number)
        ) / (
            max(pc_number) - min(pc_number)
        )

        y_norm <- (
            eigenvalues - min(eigenvalues)
        ) / (
            max(eigenvalues) - min(eigenvalues)
        )

        ## Straight line from (0,1) to (1,0).
        reference_line <- 1 - x_norm

        ## For a convex decreasing scree curve, the curve is
        ## generally above the reference line before flattening.
        knee_distance <- reference_line - y_norm

        estimated_knee <- which.max(knee_distance)
    } else {
        estimated_knee <- NA_integer_
        knee_distance <- rep(NA_real_, n_available_pcs)
    }

    ## PCs required to reach specified cumulative thresholds.
    pcs_for_threshold <- sapply(
        cumulative_thresholds,
        function(threshold) {
            hit <- which(
                cumulative_variance >= threshold
            )[1]

            if (length(hit) == 0L || is.na(hit)) {
                return(NA_integer_)
            }

            hit
        }
    )

    names(pcs_for_threshold) <- paste0(
        cumulative_thresholds * 100,
        "%"
    )

    pca_diagnostics <- data.frame(
        PC = pc_number,
        singular_value = sv$d,
        variance = pc_variance,
        eigenvalue = eigenvalues,
        PC_SD = sqrt(eigenvalues),
        explained_variance = explained_variance,
        explained_variance_percent =
            100 * explained_variance,
        cumulative_variance = cumulative_variance,
        cumulative_variance_percent =
            100 * cumulative_variance,
        knee_distance = knee_distance
    )

    ## ------------------------------------------------------------
    ## 4B. Plot PCA knee diagnostics
    ## ------------------------------------------------------------

    draw_knee_plot <- function() {
        old_par <- par(no.readonly = TRUE)
        on.exit(par(old_par), add = TRUE)

        par(
            mfrow = c(1, 2),
            mar = c(5, 5, 4, 2) + 0.1
        )

        ## Limit the displayed scree plot for readability.
        ## All PCs are still used in the calculations.
        n_plot <- min(
            n_available_pcs,
            max(100L, M + 10L)
        )

        plot_pcs <- seq_len(n_plot)

        ## Panel 1: singular-value scree plot
        plot(
            plot_pcs,
            eigenvalues[plot_pcs],
            type = "b",
            pch = 16,
            cex = 0.55,
            lwd = 1,
            col = "grey30",
            xlab = "Principal component",
            ylab = "Eigenvalue",
            main = "PCA scree / knee plot"
        )

        grid(
            col = "grey90",
            lty = 1
        )

        if (!is.na(estimated_knee) &&
            estimated_knee <= n_plot) {
            abline(
                v = estimated_knee,
                col = "#0072B2",
                lwd = 2,
                lty = 2
            )

            points(
                estimated_knee,
                eigenvalues[estimated_knee],
                pch = 19,
                cex = 1.2,
                col = "#0072B2"
            )

            text(
                estimated_knee,
                eigenvalues[estimated_knee],
                labels = paste0(
                    " Estimated knee: PC",
                    estimated_knee
                ),
                pos = 4,
                col = "#0072B2",
                cex = 0.85
            )
        }

        if (M <= n_plot) {
            abline(
                v = M,
                col = "#D55E00",
                lwd = 2,
                lty = 3
            )

            legend(
                "topright",
                legend = c(
                    paste0(
                        "Estimated knee: PC",
                        estimated_knee
                    ),
                    paste0(
                        "Requested M: ",
                        M
                    )
                ),
                col = c(
                    "#0072B2",
                    "#D55E00"
                ),
                lty = c(2, 3),
                lwd = 2,
                bty = "n",
                cex = 0.85
            )
        } else {
            legend(
                "topright",
                legend = paste0(
                    "Estimated knee: PC",
                    estimated_knee
                ),
                col = "#0072B2",
                lty = 2,
                lwd = 2,
                bty = "n",
                cex = 0.85
            )
        }

        ## Panel 2: cumulative explained variance
        plot(
            plot_pcs,
            100 * cumulative_variance[plot_pcs],
            type = "l",
            lwd = 2,
            col = "black",
            ylim = c(
                0,
                min(
                    100,
                    max(
                        100 * cumulative_variance[plot_pcs]
                    ) * 1.05
                )
            ),
            xlab = "Number of principal components",
            ylab = "Cumulative variance explained (%)",
            main = "Cumulative explained variance"
        )

        grid(
            col = "grey90",
            lty = 1
        )

        threshold_colours <- c(
            "#009E73",
            "#56B4E9",
            "#CC79A7",
            "#E69F00"
        )

        threshold_colours <- rep(
            threshold_colours,
            length.out = length(cumulative_thresholds)
        )

        for (i in seq_along(cumulative_thresholds)) {
            threshold <- cumulative_thresholds[i]
            threshold_pc <- pcs_for_threshold[i]

            abline(
                h = 100 * threshold,
                col = threshold_colours[i],
                lty = 3
            )

            if (!is.na(threshold_pc) &&
                threshold_pc <= n_plot) {
                abline(
                    v = threshold_pc,
                    col = threshold_colours[i],
                    lty = 3
                )

                points(
                    threshold_pc,
                    100 * cumulative_variance[threshold_pc],
                    pch = 19,
                    col = threshold_colours[i]
                )
            }
        }

        if (M <= n_plot) {
            abline(
                v = M,
                col = "#D55E00",
                lwd = 2,
                lty = 2
            )
        }

        legend_labels <- paste0(
            names(pcs_for_threshold),
            " variance: ",
            ifelse(
                is.na(pcs_for_threshold),
                "not reached",
                paste0(
                    "PC",
                    pcs_for_threshold
                )
            )
        )

        legend(
            "bottomright",
            legend = c(
                legend_labels,
                paste0(
                    "Requested M: ",
                    M,
                    " (",
                    round(
                        100 * cumulative_variance[M],
                        1
                    ),
                    "%)"
                )
            ),
            col = c(
                threshold_colours,
                "#D55E00"
            ),
            lty = c(
                rep(3, length(cumulative_thresholds)),
                2
            ),
            lwd = c(
                rep(1, length(cumulative_thresholds)),
                2
            ),
            bty = "n",
            cex = 0.72
        )
    }

    if (plot_knee) {
        draw_knee_plot()
    }

    if (save_knee_plot) {
        pdf(
            knee_plot_file,
            width = 12,
            height = 5.5
        )

        draw_knee_plot()

        dev.off()

        if (verbose) {
            message(
                "PCA knee plot saved to: ",
                normalizePath(
                    knee_plot_file,
                    mustWork = FALSE
                )
            )
        }
    }

    ## -- 5. coverage-shaped panel columns ----------------------------
    ## Any non-zero scale spans the same subspace; s is chosen so the columns look like plausible profiles.
    panel <- cbind(mu, sweep(U, 2L, s, "*") + mu)
    n_neg <- sum(panel < 0)
    if (n_neg > 0) {
        if (verbose) {
            message(
                "clamping ", n_neg, " negative entries (",
                round(100 * n_neg / length(panel), 2), "%)"
            )
        }
        panel[panel < 0] <- 0
    }
    colnames(panel) <- c("mean", paste0("PC", seq_len(M)))

    ## collinearity check: if library-size normalization leaked, PC1 ~ mu
    r1 <- suppressWarnings(cor(mu, U[, 1L]))
    if (!is.na(r1) && abs(r1) > 0.95) {
        warning(
            "PC1 is highly correlated with the mean profile (r = ", round(r1, 3),
            "); library-size normalization may have leaked into PCA"
        )
    }

    ## -- 6. reconstitute lNormals ------------------------------------
    split_idx <- rep(seq_along(bins_per_chr), bins_per_chr)

    lNormals <- lapply(seq_len(ncol(panel)), function(k) {
        v <- split(panel[, k], split_idx)
        out <- lapply(
            seq_along(bins_per_chr),
            function(i) {
                data.frame(records = as.numeric(v[[i]]))
            }
        )
        names(out) <- chr_names
        out
    })

    ## -- sanity checks before returning ------------------------------
    if ("records" %in% names(lNormals[[1]])) {
        stop(
            "lNormals[[1]] has a 'records' name - smoothCoverageTrack ",
            "would take the single-normal branch"
        )
    }

    len_panel <- length(unlist(lapply(lNormals[[1]], function(y) y$records)))
    if (len_panel != n_bins) {
        stop("panel member length (", len_panel, ") != number of bins (", n_bins, ")")
    }

    if (verbose) {
        message("panel members : ", length(lNormals), " (mean + ", M, " PCs)")
        message("bins          : ", n_bins)
        message("cells         : ", n_cells, "   rank ceiling = ", n_cells)
        message(
            "var explained : PC1-PC", M, " = ",
            round(100 * sum(ev[seq_len(M)]), 1), "%"
        )
        message("cor(mu, PC1)  : ", round(r1, 3))
    }

    attr(lNormals, "explained_variance") <- explained_variance
    attr(lNormals, "eigenvalues") <- eigenvalues
    attr(lNormals, "cumulative_variance") <- cumulative_variance
    attr(lNormals, "pcs_for_threshold") <- pcs_for_threshold
    attr(lNormals, "estimated_knee") <- estimated_knee
    attr(lNormals, "singular_values") <- sv$d
    attr(lNormals, "panel_matrix") <- panel
    attr(lNormals, "sampled_cells") <- idx
    lNormals
}
