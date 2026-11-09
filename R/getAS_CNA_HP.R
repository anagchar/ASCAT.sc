## Allele-specific copy number from haplotagged reads (HP tag), per bin.
## Counterpart of the SNP-based getAS_CNA(): the inputs are the per-bin
## hap1/hap2 read counts that getCoverageTrack.10XBAM(hp_tag="HP") adds
## to each cell's track.

## BAF of one segment from per-bin hap1/hap2 counts, robust to phase switch errors.
##
## Each block of `block` consecutive bins contributes
##     P(h1 | n, p) + P(h1 | n, 1 - p)
## i.e. the haplotype labels of a block are treated as exchangeable. A phase
## switch between blocks therefore cannot bias the estimate, but the labels
## (which haplotype is the amplified one) are not used. The result is the
## folded BAF in [0.5, 1].
##   block = 1  : fully switch-robust, least power (sparse single-cell counts)
##   block > 1  : more power per term; a switch inside a block dilutes it
## Returns c(lower 5%, best, upper 95%); NA if the segment has no phased reads.
fitBinom.hap <- function(h1,
                         h2,
                         block=1,
                         maxdepth=1000,
                         values=seq(.5, 1, .001))
{
    ok <- !is.na(h1) & !is.na(h2) & (h1 + h2) > 0
    h1 <- h1[ok]
    h2 <- h2[ok]
    if (length(h1) == 0) return(rep(NA_real_, 3))
    grp <- ceiling(seq_along(h1) / block)
    counts <- as.numeric(tapply(h1, grp, sum))
    depths <- as.numeric(tapply(h1 + h2, grp, sum))
    big <- depths > maxdepth
    counts[big] <- round(counts[big] / depths[big] * maxdepth)
    depths[big] <- maxdepth
    llh <- vapply(values, function(p)
    {
        a <- dbinom(counts, size=depths, prob=p, log=TRUE)
        b <- dbinom(counts, size=depths, prob=1 - p, log=TRUE)
        m <- pmax(a, b)                       # log-sum-exp: no underflow at high depth
        m[is.infinite(m)] <- 0                # both terms impossible (e.g. p = 1, mixed bin): -Inf, not NaN
        sum(m + log(exp(a - m) + exp(b - m)))
    }, numeric(1))
    llh <- llh - max(llh)
    post <- exp(llh) / sum(exp(llh))
    cs <- cumsum(post)
    c(values[which(cs > .05)[1]], values[which.max(llh)], values[which(cs > .95)[1]])
}
