## ===========================================================================
## Negative binomial significance test for chromosome-level CNA in single cells
##
## ASCAT.sc fits copy number to a continuous logR and has no notion of
## statistical significance: no p-value, no false-positive control. In
## single-cell WGS/ATAC on normal tissue, where real alterations are rare and
## every cell is noisy, that makes a real event indistinguishable from
## background. This adds a calibrated null.
##
##   X_cr ~ NegBinom(mean = alpha_c * beta_r * gamma_cr, dispersion = phi)
##   Var(X) = mu + phi*mu^2   ->   R size parameter = 1/phi
##
## alpha_c is the cell's library size, beta_r the chromosome's share of reads,
## gamma_cr the copy ratio being tested. sqrt(phi) is the irreducible
## cell-to-cell CV: how much two identical healthy cells differ, which is what
## sets how far a cell must deviate before the deviation is believable.
## ===========================================================================


## --- internals --------------------------------------------------------------

## res$allTracks can be unnamed (only positionally aligned with allProfiles).
.nbcna_tracks <- function(res) {
    tracks <- res$allTracks
    if (is.null(tracks)) stop("res$allTracks not found")
    if (is.null(names(tracks))) {
        nm <- names(res$allProfiles)
        if (is.null(nm)) nm <- res$barcodes
        if (length(nm) != length(tracks)) {
            stop("res$allTracks is unnamed and cannot be matched to allProfiles/barcodes")
        }
        names(tracks) <- nm
    }
    tracks
}

## Null expectation with leave-one-out library size: rest_cr = lib_c - x_cr
## estimates alpha_c * (1 - beta_r), so mu_cr = rest_cr * beta_r/(1 - beta_r).
## Using the plain library size would put x_cr inside its own denominator,
## which shrinks the apparent ratio of a real event on a large chromosome by
## several percent - the same order as the CV we are trying to get under 10%.
.nbcna_expected <- function(counts, libsize) {
    beta <- apply(counts / libsize, 2, median, na.rm = TRUE)
    beta <- beta / sum(beta, na.rm = TRUE)
    mu <- t(t(libsize - counts) * (beta / (1 - beta)))
    dimnames(mu) <- dimnames(counts)
    list(mu = mu, beta = beta)
}

.nbcna_qnb <- function(p, mu, phi) {
    qnbinom(p, size = 1 / pmax(phi, 1e-8), mu = mu)
}

## merge adjacent bins that are too thin (only bites when quantile breaks tie)
.nbcna_merge_bins <- function(idx, min_per_bin) {
    ub <- sort(unique(idx))
    tab <- table(idx)
    grp <- integer(length(ub))
    g <- 1L
    n <- 0L
    for (i in seq_along(ub)) {
        grp[i] <- g
        n <- n + as.integer(tab[[as.character(ub[i])]])
        if (n >= min_per_bin) {
            g <- g + 1L
            n <- 0L
        }
    }
    last <- grp[length(grp)]
    if (sum(idx %in% ub[grp == last]) < min_per_bin && last > 1L) {
        grp[grp == last] <- last - 1L
    }
    grp[match(idx, ub)]
}

## variance lost to symmetric trimming, normal theory (0.623 at trim=0.05)
.nbcna_trim_factor <- function(a) {
    if (a <= 0) return(1)
    z <- qnorm(1 - a)
    1 - 2 * z * dnorm(z) / (1 - 2 * a)
}

## two-sided tail probability, Poisson limit when phi collapses to 0
.nbcna_pval <- function(x, mu, phi) {
    p <- numeric(length(x))
    ip <- !is.finite(phi) | phi <= 0
    if (any(ip)) {
        lo <- ppois(x[ip], mu[ip])
        hi <- ppois(x[ip] - 1, mu[ip], lower.tail = FALSE)
        p[ip] <- pmin(1, 2 * pmin(lo, hi))
    }
    if (any(!ip)) {
        sz <- 1 / phi[!ip]
        lo <- pnbinom(x[!ip], size = sz, mu = mu[!ip])
        hi <- pnbinom(x[!ip] - 1, size = sz, mu = mu[!ip], lower.tail = FALSE)
        p[!ip] <- pmin(1, 2 * pmin(lo, hi))
    }
    p
}


## --- 1. cell x chromosome raw counts ----------------------------------------

#' @title Cell x chromosome read counts from an ASCAT.sc result object
#' @description Uses \code{res$allTracks[[cell]]$nlCTS.tumour}, the rebinned but
#'   UNCORRECTED counts. The 'smoothed' track is not usable here: the negative
#'   binomial models counts, not log-ratios.
#' @param res An ASCAT.sc result object
#' @param chrs Chromosomes to keep (default: autosomes)
#' @param cells Cells to keep (default: all in \code{res$allTracks})
#' @return List with \code{$counts} (cells x chromosomes) and \code{$libsize}
#' @export
nbcna_chrom_counts <- function(res, chrs = NULL, cells = NULL) {
    tracks <- .nbcna_tracks(res)
    if (is.null(cells)) cells <- names(tracks)

    .nl <- function(cell) {
        tr <- tracks[[cell]]
        if (!is.null(tr$nlCTS.tumour)) tr$nlCTS.tumour else tr[[2]]
    }

    all_chrs <- names(.nl(cells[1]))
    if (is.null(chrs)) {
        clean <- gsub("^chr", "", all_chrs)
        chrs <- all_chrs[clean %in% as.character(1:22)]
        if (length(chrs) == 0) chrs <- all_chrs
    }

    counts <- t(sapply(cells, function(cell) {
        nl <- .nl(cell)
        sapply(chrs, function(ch) sum(nl[[ch]]$records, na.rm = TRUE))
    }))
    dimnames(counts) <- list(cells, chrs)
    list(counts = counts, libsize = rowSums(counts, na.rm = TRUE))
}


## --- 2. drop broken barcodes ------------------------------------------------

#' @title Drop cells whose chromosomes disagree with each other
#' @description Empty droplets, doublets and read pile-ups are extreme on
#'   several chromosomes at once and inflate the dispersion estimate. Scores
#'   each cell by the MAD of its standardised residuals
#'   \code{(x - mu)/sqrt(mu + phi0*mu^2)} across chromosomes: MAD has a 50\%
#'   breakdown point, so one genuinely altered chromosome out of 22 does not
#'   get the cell excluded. The denominator is the full negative binomial sd,
#'   using a robust pilot \code{phi0}, so that the score does not track
#'   coverage: a raw ratio spread would drop low-coverage cells and a
#'   Poisson-only \code{sqrt(mu)} would drop high-coverage ones, in both cases
#'   instead of the broken barcodes.
#' @param m Output of \code{\link{nbcna_chrom_counts}}
#' @param cell_types Optional named vector cell -> type. Scores are computed
#'   within type, so genuine between-type accessibility differences do not
#'   count as within-cell disagreement.
#' @param nmad Cells above \code{median + nmad * mad} of the score are dropped
#' @param min_libsize Library size floor
#' @return List with \code{$keep} (cell names) and \code{$stats} (per-cell
#'   diagnostics, including \code{mad_ratio} on the raw ratio scale)
#' @export
nbcna_filter_cells <- function(m, cell_types = NULL, nmad = 3, min_libsize = 5000) {
    counts <- m$counts
    libsize <- m$libsize
    cells <- rownames(counts)

    grp <- if (is.null(cell_types)) {
        setNames(rep("all", length(cells)), cells)
    } else {
        setNames(as.character(cell_types[cells]), cells)
    }
    grp[is.na(grp)] <- "_unlabelled"

    stats <- do.call("rbind", lapply(sort(unique(grp)), function(g) {
        gc <- cells[grp == g]
        cnt <- counts[gc, , drop = FALSE]
        lib <- libsize[gc]
        mu <- .nbcna_expected(cnt, lib)$mu
        ratio <- cnt / mu
        ## robust pilot dispersion: mad() not var(), so the broken cells we are
        ## about to remove do not set the scale used to find them
        phi0 <- max(0, mad(as.vector(ratio), na.rm = TRUE)^2 -
                       mean(1 / as.vector(mu), na.rm = TRUE))
        mad_z <- apply((cnt - mu) / sqrt(mu + phi0 * mu^2), 1, mad, na.rm = TRUE)
        mad_ratio <- apply(ratio, 1, mad, na.rm = TRUE)
        thr <- median(mad_z, na.rm = TRUE) + nmad * mad(mad_z, na.rm = TRUE)
        data.frame(
            cell = gc, cell_type = g, libsize = lib,
            mad_z = mad_z, mad_ratio = mad_ratio, threshold = thr,
            keep = is.finite(mad_z) & mad_z < thr & lib >= min_libsize,
            stringsAsFactors = FALSE
        )
    }))
    rownames(stats) <- NULL

    n_mad <- sum(!is.finite(stats$mad_z) | stats$mad_z >= stats$threshold)
    n_lib <- sum(stats$libsize < min_libsize)
    message(
        "kept ", sum(stats$keep), "/", nrow(stats), " cells (",
        n_mad, " with disagreeing chromosomes, ", n_lib, " below library floor)"
    )
    list(keep = stats$cell[stats$keep], stats = stats)
}


## --- 3. cell type labels ----------------------------------------------------

#' @title Map RNA cell-type labels onto ASCAT.sc track names
#' @description Coverage tracks are keyed by ATAC barcode (wrapped in the BAM
#'   name), cell types are called on the GEX side. This walks
#'   track name -> ATAC barcode -> GEX barcode -> label and returns a vector
#'   named by \code{names(res$allTracks)}, ready to pass to
#'   \code{\link{nbcna_test}} as \code{cell_types}. Different cell types have
#'   genuinely different chromosome accessibility, and pooling them folds that
#'   difference into the noise estimate. Even 2-4 coarse groups recover most of
#'   the benefit: the labels need to be correct partitions, not biologically
#'   meaningful ones.
#' @param res An ASCAT.sc result object
#' @param annotation data.frame or tab-separated file with a GEX barcode column
#'   and a label column
#' @param barcode_map data.frame or CSV file with columns \code{ATAC} and \code{GEX}
#' @param barcode_col Name of the GEX barcode column in \code{annotation}
#' @param type_col Name of the label column in \code{annotation}
#' @return Named character vector: track name -> cell type, matched cells only
#' @export
nbcna_map_celltypes <- function(res,
                                annotation,
                                barcode_map,
                                barcode_col = "barcode",
                                type_col = "cell_type") {
    ann <- if (is.character(annotation)) {
        read.delim(annotation, stringsAsFactors = FALSE)
    } else annotation
    map <- if (is.character(barcode_map)) {
        read.csv(barcode_map, stringsAsFactors = FALSE)
    } else barcode_map

    if (!all(c("ATAC", "GEX") %in% colnames(map))) {
        stop("barcode_map needs columns ATAC and GEX")
    }
    if (!all(c(barcode_col, type_col) %in% colnames(ann))) {
        stop("annotation needs columns ", barcode_col, " and ", type_col)
    }

    .strip <- function(b) sub("-\\d+$", "", b)
    tracks <- names(.nbcna_tracks(res))
    atac <- .strip(sub("^.*\\.bam_", "", tracks))

    atac2gex <- setNames(.strip(map$GEX), .strip(map$ATAC))
    gex2type <- setNames(as.character(ann[[type_col]]), .strip(ann[[barcode_col]]))

    types <- setNames(unname(gex2type[atac2gex[atac]]), tracks)
    types <- types[!is.na(types)]
    message(
        "labelled ", length(types), "/", length(tracks), " tracks across ",
        length(unique(types)), " types"
    )
    print(table(types))
    types
}


## --- 4. dispersion, estimated in quantile bins of the mean ------------------

#' @title Fit the negative binomial dispersion as a function of expected count
#' @description Bins are quantiles of \code{mu}, not equal-width. Coverage spans
#'   ~1000x, so equal-width bins collapse: nearly every observation lands in the
#'   first bin and one dispersion figure gets applied across the whole range.
#' @param mu Vector of expected counts
#' @param x Vector of observed counts
#' @param method \code{"ratio"} estimates \code{phi = Var(x/mu) - E[1/mu]}, which
#'   stays unbiased when \code{mu} varies within a bin; \code{"count"} is the
#'   textbook \code{(Var(x) - mean(x))/mean(x)^2}, valid only at constant mu
#' @param nbins Number of quantile bins
#' @param min_per_bin Minimum observations per bin after merging
#' @param trim Symmetric trimming fraction, so real events do not enter the null
#' @return data.frame, one row per bin, with \code{phi} and \code{cv}
#' @export
nbcna_fit_dispersion <- function(mu, x,
                                 method = c("ratio", "count"),
                                 nbins = 20,
                                 min_per_bin = 200,
                                 trim = 0.05) {
    method <- match.arg(method)
    ok <- is.finite(mu) & is.finite(x) & mu > 0
    mu <- mu[ok]
    x <- x[ok]
    if (length(mu) < min_per_bin) {
        stop("too few observations (", length(mu), ") to estimate dispersion")
    }

    tfac <- .nbcna_trim_factor(trim)
    breaks <- unique(quantile(mu, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE))
    idx <- cut(mu, breaks, include.lowest = TRUE, labels = FALSE)
    grp <- .nbcna_merge_bins(idx, min_per_bin)

    out <- do.call("rbind", lapply(sort(unique(grp)), function(g) {
        keep <- grp == g
        xs <- x[keep]
        ms <- mu[keep]
        v <- if (method == "ratio") xs / ms else xs
        if (trim > 0 && length(v) >= 20) {
            qs <- quantile(v, c(trim, 1 - trim), na.rm = TRUE)
            inb <- v >= qs[1] & v <= qs[2]
            xs <- xs[inb]
            ms <- ms[inb]
            v <- v[inb]
        }
        if (length(xs) < 3) return(NULL)
        mbar <- mean(xs)
        if (method == "ratio") {
            phi <- var(v) / tfac - mean(1 / ms)
            vbar <- var(xs) / tfac
        } else {
            vbar <- var(xs)
            phi <- (vbar - mbar) / mbar^2
        }
        data.frame(
            mu_lo = min(mu[keep]), mu_hi = max(mu[keep]),
            n = sum(keep), n_used = length(xs),
            mean = mbar, var = vbar, var_poisson = mbar,
            phi_raw = phi, phi = max(phi, 0), cv = sqrt(max(phi, 0))
        )
    }))
    rownames(out) <- NULL
    if (any(out$phi_raw < 0)) {
        warning(
            sum(out$phi_raw < 0), "/", nrow(out), " bins gave negative ",
            "dispersion; floored at 0 (treated as Poisson)"
        )
    }
    out
}

#' @title Look up the fitted dispersion at a given expected count
#' @param fit Output of \code{\link{nbcna_fit_dispersion}}
#' @param mu Vector of expected counts
#' @export
nbcna_phi_at <- function(fit, mu) {
    i <- findInterval(mu, fit$mu_lo, rightmost.closed = TRUE, all.inside = TRUE)
    fit$phi[pmin(pmax(i, 1), nrow(fit))]
}


## --- 5. main entry point ----------------------------------------------------

#' @title Negative binomial test for chromosome-level gains and losses
#' @description Runs the whole pipeline: counts, broken-barcode filtering,
#'   per-cell-type dispersion, two-sided NB tail probabilities, BH correction.
#' @param res An ASCAT.sc result object
#' @param cell_types Named vector cell -> type, e.g. from
#'   \code{\link{nbcna_map_celltypes}}. Must be a vector, not a data.frame.
#'   Leaving this NULL pools unlike cell types and inflates phi.
#' @param chrs,cells Passed to \code{\link{nbcna_chrom_counts}}
#' @param filter Drop broken barcodes with \code{\link{nbcna_filter_cells}}
#' @param nmad,min_libsize Passed to \code{\link{nbcna_filter_cells}}
#' @param min_cells Cell types with fewer cells are dropped
#' @param method,nbins,min_per_bin,trim Passed to \code{\link{nbcna_fit_dispersion}}
#' @param alpha FDR threshold
#' @return List with \code{$results} (one row per cell x chromosome test),
#'   \code{$fits}, \code{$beta}, \code{$counts}, \code{$libsize},
#'   \code{$cell_types}, \code{$filter}
#' @export
nbcna_test <- function(res,
                       cell_types = NULL,
                       chrs = NULL,
                       cells = NULL,
                       filter = TRUE,
                       nmad = 3,
                       min_libsize = 5000,
                       min_cells = 30,
                       method = c("ratio", "count"),
                       nbins = 20,
                       min_per_bin = 200,
                       trim = 0.05,
                       alpha = 0.05) {
    method <- match.arg(method)
    if (!is.null(cell_types) && !is.null(dim(cell_types))) {
        stop("cell_types must be a named vector (cell -> type), not a ",
             class(cell_types)[1], "; see nbcna_map_celltypes()")
    }

    mat <- nbcna_chrom_counts(res, chrs = chrs, cells = cells)
    counts <- mat$counts
    libsize <- mat$libsize
    chrnames <- colnames(counts)

    if (is.null(cell_types)) {
        warning(
            "cell_types is NULL: treating all cells as one group. ",
            "This is the correct model if the sample is a single population. ",
            "If it is not, between-type accessibility differences are ",
            "absorbed into phi and cost sensitivity."
        )
        grp <- setNames(rep("all", nrow(counts)), rownames(counts))
    } else {
        grp <- setNames(as.character(cell_types[rownames(counts)]), rownames(counts))
        if (any(is.na(grp))) message(sum(is.na(grp)), " cells unlabelled, dropped")
    }

    ## restrict to labelled cells before filtering, so that beta and the MAD
    ## threshold are set by the cells actually being tested
    lab <- !is.na(grp) & libsize > 0
    counts <- counts[lab, , drop = FALSE]
    libsize <- libsize[lab]
    grp <- grp[lab]

    filt <- NULL
    if (filter) {
        filt <- nbcna_filter_cells(list(counts = counts, libsize = libsize),
                                   cell_types = grp, nmad = nmad,
                                   min_libsize = min_libsize)
        keep <- rownames(counts) %in% filt$keep
        counts <- counts[keep, , drop = FALSE]
        libsize <- libsize[keep]
        grp <- grp[keep]
    }

    groups <- names(which(table(grp) >= min_cells))
    if (length(groups) == 0) stop("no group has >= ", min_cells, " cells")
    if (length(setdiff(unique(grp), groups))) {
        message("dropping small group(s): ",
                paste(setdiff(unique(grp), groups), collapse = ", "))
    }

    fits <- betas <- res_list <- list()
    for (g in groups) {
        gcells <- names(grp)[grp == g]
        cnt <- counts[gcells, , drop = FALSE]
        lib <- libsize[gcells]

        expec <- .nbcna_expected(cnt, lib)
        mu <- expec$mu

        fit <- nbcna_fit_dispersion(as.vector(mu), as.vector(cnt),
                                    method = method, nbins = nbins,
                                    min_per_bin = min_per_bin, trim = trim)
        phi <- matrix(nbcna_phi_at(fit, as.vector(mu)), nrow = nrow(mu))
        pv <- matrix(.nbcna_pval(as.vector(cnt), as.vector(mu), as.vector(phi)),
                     nrow = nrow(mu))

        res_list[[g]] <- data.frame(
            cell = rep(gcells, times = ncol(cnt)),
            cell_type = g,
            chr = rep(chrnames, each = nrow(cnt)),
            observed = as.vector(cnt),
            libsize = rep(lib, times = ncol(cnt)),
            beta = rep(expec$beta, each = nrow(cnt)),
            expected = as.vector(mu),
            phi = as.vector(phi),
            ratio = as.vector(cnt) / as.vector(mu),
            cn = 2 * as.vector(cnt) / as.vector(mu),
            pvalue = as.vector(pv),
            stringsAsFactors = FALSE
        )
        fits[[g]] <- fit
        betas[[g]] <- expec$beta
    }

    out <- do.call("rbind", res_list)
    rownames(out) <- NULL
    out$padj <- p.adjust(out$pvalue, method = "BH")
    out$call <- ifelse(out$padj >= alpha, "neutral",
                       ifelse(out$ratio > 1, "gain", "loss"))
    message(
        sum(out$padj < alpha), "/", nrow(out), " tests significant at FDR ",
        alpha, " (", round(100 * mean(out$padj < alpha), 2), "%)"
    )

    list(results = out, fits = fits, beta = betas, counts = counts,
         libsize = libsize, cell_types = grp, filter = filt, alpha = alpha)
}


## --- 6. what the data can resolve -------------------------------------------

#' @title Report cell-to-cell variability and what it makes detectable
#' @description \code{cv_pct} is \code{sqrt(phi)}, the irreducible cell-to-cell
#'   CV: a clean single population runs 8-10\%. The \code{min_cn_*} columns
#'   turn that into the smallest whole-chromosome gain and loss that could
#'   clear the multiple-testing threshold, evaluated at BH's most stringent
#'   rung (\code{alpha/n_tests}) - that is, what it takes to find a single
#'   isolated event in an otherwise quiet dataset. \code{cv_for_cn3} inverts
#'   the question: the CV that would be needed to bring an ordinary trisomy
#'   into reach at that coverage. \code{NA} there means no amount of noise
#'   reduction is enough and the limit is coverage, not QC. All of this is a
#'   property of the data alone and can be read before doing any biology.
#' @param x Output of \code{\link{nbcna_test}}, or a single dispersion fit
#' @param n_tests Number of tests corrected for (default: as run)
#' @param alpha Significance level
#' @export
nbcna_variability_report <- function(x, n_tests = NULL, alpha = 0.05) {
    fits <- if (is.data.frame(x)) list(all = x) else x$fits
    if (is.null(n_tests)) {
        n_tests <- if (is.data.frame(x)) 1 else nrow(x$results)
    }
    a2 <- alpha / (2 * n_tests)

    out <- do.call("rbind", lapply(names(fits), function(g) {
        f <- fits[[g]]
        sd_nb <- sqrt(f$mean + f$phi * f$mean^2)
        sd_po <- sqrt(f$mean)
        gain <- 2 * (.nbcna_qnb(1 - a2, f$mean, f$phi) + 1) / f$mean
        loss <- 2 * (.nbcna_qnb(a2, f$mean, f$phi) - 1) / f$mean
        ## CV that would put a trisomy (1.5x) exactly at the threshold:
        ## 1.5*mu = mu + z*sqrt(mu + phi*mu^2)  ->  phi = 0.25/z^2 - 1/mu
        phi_need <- 0.25 / qnorm(1 - a2)^2 - 1 / f$mean
        data.frame(
            cell_type = g, mu = round(f$mean), n = f$n_used,
            phi = signif(f$phi, 3), cv_pct = round(100 * sqrt(f$phi), 1),
            sd_poisson = round(sd_po, 1), sd_nb = round(sd_nb, 1),
            widening = round(sd_nb / sd_po, 2),
            min_cn_gain = round(gain, 2), min_cn_loss = round(pmax(loss, 0), 2),
            cv_for_cn3 = ifelse(phi_need > 0,
                                round(100 * sqrt(pmax(phi_need, 0)), 1), NA),
            stringsAsFactors = FALSE
        )
    }))
    rownames(out) <- NULL
    cat("\n  phi         : overdispersion (0 = pure Poisson)\n",
        " cv_pct      : sqrt(phi) = irreducible cell-to-cell CV\n",
        " widening    : how many times wider the NB null is than Poisson\n",
        " min_cn_gain : smallest isolated gain, in copies, findable over ",
        n_tests, " tests\n",
        " cv_for_cn3  : CV needed to bring a trisomy into reach",
        " (NA = coverage-limited)\n\n", sep = "")
    print(out)
    cat("\nMedian CV across bins: ",
        round(100 * median(sqrt(unlist(lapply(fits, `[[`, "phi")))), 1), "%\n",
        sep = "")
    invisible(out)
}

#' @title Mean-variance and CV diagnostics for one cell type
#' @param x Output of \code{\link{nbcna_test}}, or a single dispersion fit
#' @param cell_type Which group to plot (default: the first)
#' @export
nbcna_plot_dispersion <- function(x, cell_type = NULL) {
    fits <- if (is.data.frame(x)) list(all = x) else x$fits
    g <- if (is.null(cell_type)) names(fits)[1] else cell_type
    f <- fits[[g]]
    op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
    on.exit(par(op))

    plot(f$mean, f$var, log = "xy", pch = 19,
         xlab = "expected count (mu)", ylab = "observed variance",
         main = paste0("mean-variance: ", g))
    xs <- exp(seq(log(min(f$mean)), log(max(f$mean)), length.out = 200))
    phis <- nbcna_phi_at(f, xs)
    lines(xs, xs, col = "grey50", lty = 2, lwd = 2)
    lines(xs, xs + phis * xs^2, col = "firebrick", lwd = 2)
    legend("topleft", bty = "n", lwd = 2, lty = c(2, 1),
           col = c("grey50", "firebrick"), legend = c("Poisson", "fitted NB"))

    cv <- sqrt(f$var) / f$mean
    plot(f$mean, 100 * cv, log = "x", pch = 19, ylim = c(0, max(100 * cv) * 1.1),
         xlab = "expected count (mu)", ylab = "CV (%)",
         main = "cell-to-cell variability")
    lines(xs, 100 / sqrt(xs), col = "grey50", lty = 2, lwd = 2)
    lines(xs, 100 * sqrt(1 / xs + phis), col = "firebrick", lwd = 2)
    abline(h = 100 * sqrt(median(f$phi)), col = "steelblue", lty = 3, lwd = 2)
    legend("topright", bty = "n", lwd = 2, lty = c(2, 1, 3),
           col = c("grey50", "firebrick", "steelblue"),
           legend = c("Poisson only", "fitted NB", "sqrt(phi) floor"))
    invisible(f)
}


## --- 7. is the null calibrated ----------------------------------------------

#' @title Check whether the fitted null is calibrated
#' @description Under a correct null the p-values are uniform. A left-heavy
#'   histogram means the null is too tight and the hits are mostly false; a
#'   histogram sagging in the middle means it is too loose. \code{pi0} is the
#'   Storey estimate of the fraction of true nulls: near 1 means nothing beyond
#'   noise is present. Counts are discrete and the two-sided p-value doubles
#'   the smaller tail, so a spike at p = 1 is expected and is not
#'   miscalibration.
#' @param x Output of \code{\link{nbcna_test}}
#' @param lambda Tail cut used for the pi0 estimate
#' @param plot Draw the p-value histogram
#' @export
nbcna_calibration <- function(x, lambda = 0.5, plot = TRUE) {
    p <- x$results$pvalue
    padj <- x$results$padj
    alpha <- x$alpha
    pi0 <- min(1, mean(p > lambda, na.rm = TRUE) / (1 - lambda))
    nsig <- sum(padj < alpha, na.rm = TRUE)

    if (plot) {
        hist(p, breaks = 50, col = "grey85", border = "white",
             xlab = "p-value", main = "null calibration")
        abline(h = pi0 * length(p) / 50, col = "firebrick", lwd = 2, lty = 2)
        legend("topright", bty = "n", lwd = 2, lty = 2, col = "firebrick",
               legend = paste0("expected under null (pi0 = ", round(pi0, 3), ")"))
    }

    cat("\n tests            : ", length(p),
        "\n significant      : ", nsig, " (", round(100 * nsig / length(p), 2), "%)",
        "\n pi0 (true nulls) : ", round(pi0, 3),
        "\n excess over null : ", round(100 * (1 - pi0), 2), "%\n\n", sep = "")
    invisible(list(pi0 = pi0, n_tests = length(p), n_sig = nsig))
}
