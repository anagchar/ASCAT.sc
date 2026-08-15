## Absolute-ploidy validation for single cells from allele-specific evidence.
##
## Coverage alone cannot fix the absolute scale: logR only constrains copy number
## up to a multiplicative constant. BAF adds the missing constraint, because a
## segment's allelic ratio must decompose into non-negative INTEGER major/minor
## copies. A wrong scale factor pushes those implied integers off the lattice.
##
## What this can and cannot decide
## -------------------------------
## Doubling maps the integer lattice onto itself: if (nA,nB) are integers then so
## are (2nA,2nB), and BAF is unchanged. So allelic evidence can NEVER exclude a
## whole-genome doubling. It can exclude HALVING, because odd copy numbers do not
## survive it (2+1 -> 1+0.5). This function therefore reports the smallest ploidy
## consistent with the data and states plainly that multiples remain possible;
## it does not pretend to a resolution the data cannot carry.
##
## Two limits measured on ATLAS/Skin/DNA/donor0
## -------------------------------------------
## 1. Only segments with imbalanced BAF carry scale information. A cell whose
##    segments are all balanced (BAF ~ 0.5) is genuinely unidentifiable and is
##    reported as such rather than silently accepted.
##
## 2. Within a cell, reads arrive in input DNA molecules and every read from one
##    molecule carries the same haplotype, so BAF precision is set by the
##    MOLECULE count, not the read count. Measured overdispersion on donor0 was
##    ~5 reads per independent molecule (variance 5x binomial; obs/exp |0.5-f|
##    ratio 2.23 over 1 Mb blocks). Binomial intervals from fitBinom.1dist are
##    therefore ~2.2x too narrow and are widened by sqrt(reads_per_molecule).
##    The same measurement implies per-cell BAF is only informative over large
##    segments (SE ~ 0.17/sqrt(Mb)), hence the 10 Mb min_seg_width default.

nbcna_absolutePloidy = function(nprof,
                                ploidy,
                                purity = 1,
                                scales = c(0.5, 1, 2),
                                min_seg_width = 1e7,
                                reads_per_molecule = 5,
                                min_imbalance = 0.05,
                                fit_cutoff = 3)
{
    stopifnot(all(c("logr", "BAF", "startpos", "endpos") %in% colnames(nprof)))
    scales = sort(scales)

    width = as.numeric(nprof[, "endpos"]) - as.numeric(nprof[, "startpos"])
    baf = as.numeric(nprof[, "BAF"])
    logr = as.numeric(nprof[, "logr"])

    ## BAF sd from the stored 90% interval, widened for molecule-level
    ## overdispersion (see header note 2)
    if(all(c("q05", "q95") %in% colnames(nprof)))
    {
        sdbaf = (as.numeric(nprof[, "q95"]) - as.numeric(nprof[, "q05"])) /
            (2 * qnorm(0.95))
        sdbaf = sdbaf * sqrt(reads_per_molecule)
    }
    else sdbaf = rep(0.05, nrow(nprof))
    sdbaf[is.na(sdbaf) | sdbaf <= 0] = NA

    sdlogr = if("logr.sd" %in% colnames(nprof))
        as.numeric(nprof[, "logr.sd"]) else rep(0, nrow(nprof))
    sdlogr[is.na(sdlogr)] = 0

    .fail = function(reason, ninf)
        list(verdict = "ambiguous", reason = reason,
             ploidy.called = ploidy, ploidy.min = NA, scale = NA,
             identifiable = FALSE, doubling.excluded = FALSE,
             smaller.scale.fits = FALSE, n.informative = ninf,
             scores = setNames(rep(NA_real_, length(scales)), scales),
             segments = NULL)

    keep = !is.na(baf) & !is.na(logr) & !is.na(sdbaf) & width >= min_seg_width
    if(sum(keep) < 2)
        return(.fail("fewer than 2 segments wide enough for a usable per-cell BAF", 0))

    ## only imbalanced segments carry scale information (see header)
    inf = keep & abs(baf - 0.5) > min_imbalance
    if(sum(inf) < 2)
        return(.fail(paste0("all usable segments are allelically balanced; ",
                            "ploidy and any multiple of it are indistinguishable"),
                     sum(inf)))

    ## Reduced chi-square of the implied (major,minor) against the integer
    ## lattice, propagating both BAF and logR error into the copy numbers.
    .score = function(s)
    {
        ntot = transform_bulk2tumour(logr[inf], purity, ploidy * s)
        nmaj = baf[inf] * (purity * ntot + 2 * (1 - purity)) - (1 - purity)
        nmin = ntot - nmaj
        varn = (ntot * log(2) * sdlogr[inf])^2            # from logR
        vmaj = baf[inf]^2 * varn + (ntot * sdbaf[inf])^2  # + from BAF
        vmin = (1 - baf[inf])^2 * varn + (ntot * sdbaf[inf])^2
        z = ((nmaj - round(nmaj))^2 / pmax(vmaj, 1e-8) +
             (nmin - round(nmin))^2 / pmax(vmin, 1e-8)) / 2
        ## a negative minor copy number is impossible, not merely a poor fit
        z[nmin < -0.5 | ntot < -0.5] = 1e6
        w = width[inf] / sum(width[inf])
        sum(z * w)
    }

    scores = sapply(scales, .score)
    names(scores) = scales
    fits = is.finite(scores) & scores <= fit_cutoff

    ## The job is to VALIDATE the coverage call, so the verdict turns on whether
    ## the called ploidy itself survives. Picking the smallest fitting scale
    ## instead would bias downward whenever BAF is noisy: wide intervals make
    ## every scale "fit", and the minimum would always win.
    called_fits = fits[which(scales == 1)]
    if(length(called_fits) == 0)
        stop("scales must include 1 so the called ploidy can be tested")

    verdict = if(called_fits) "accept" else if(any(fits)) "rescale" else "drop"
    best = if(called_fits) 1 else if(any(fits)) scales[which(fits)[1]] else NA

    ## did the allelic data exclude anything at all?
    identifiable = any(!fits)
    ## a smaller ploidy that also fits is worth surfacing, but is not grounds to
    ## overrule the call on its own
    smaller_fits = any(fits & scales < 1)

    if(is.na(best))
        return(list(verdict = "drop",
                    reason = "no tested scale puts the profile on the integer lattice",
                    note = "profile is not integer-consistent at any scale; treat the cell as unreliable",
                    ploidy.called = ploidy, ploidy.min = NA, scale = NA,
                    identifiable = identifiable, doubling.excluded = FALSE,
                    smaller.scale.fits = smaller_fits,
                    n.informative = sum(inf), scores = scores, segments = NULL))

    ntot = transform_bulk2tumour(logr, purity, ploidy * best)
    nmaj = baf * (purity * ntot + 2 * (1 - purity)) - (1 - purity)

    list(verdict = verdict,
         reason = if(best == 1)
             "the called ploidy is consistent with the allelic data"
         else
             paste0("BAF excludes the called ploidy; nearest consistent scale ",
                    "is ", best, "x"),
         note = paste0("whole-genome doubling cannot be excluded from allelic ",
                       "data, so ploidy.min is a lower bound",
                       if(smaller_fits)
                           "; a smaller ploidy also fits - inspect $scores"
                       else ""),
         ploidy.called = ploidy,
         ploidy.min = ploidy * min(scales[fits]),
         scale = best,
         identifiable = identifiable,
         doubling.excluded = FALSE,
         smaller.scale.fits = smaller_fits,
         n.informative = sum(inf),
         scores = scores,
         segments = data.frame(chr = nprof[, 1],
                               startpos = nprof[, "startpos"],
                               endpos = nprof[, "endpos"],
                               logr = logr,
                               BAF = baf,
                               sdBAF = sdbaf,
                               used = keep,
                               informative = inf,
                               ntot = ntot,
                               nMajor = nmaj,
                               nMinor = ntot - nmaj,
                               stringsAsFactors = FALSE))
}
