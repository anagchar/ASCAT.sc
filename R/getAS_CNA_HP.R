fitBinom.hap = function(h1, h2, block = 1, maxdepth = 1000, values=seq(.5, 1, .001))
{
    ok = !is.na(h1) & !is.na(h2) & (h1+h2) > 0
    h1 = h1[ok]
    h2 = h2[ok]
    if(length(h1) == 0) return(rep(NA_real_, 3))
    grp = ceiling(seq_along(h1)/block)
    counts = as.numeric(tapply(h1, grp, sum))
    depths = as.numeric(tapply(h1+h2, grp, sum))
    big = depths > maxdepth
    counts[big] = round(counts[big] / depths[big] * maxdepth)
    depths[big] = maxdepth
    llh = vapply(values, function(p)
    {
        a = dbinom(counts, depths, p, log=TRUE)
        b = dbinom(counts, depths, 1-p, log=TRUE)
        m = pmax(a, b)
        m[is.infinite(m)] = 0
        sum(m + log(exp(a-m) + exp(b-m)))
    }, numeric(1))
    llh = llh - max(llh)
    post = exp(llh) / sum(exp(llh))
    cs = cumsum(post)
    c(values[which(cs > .05)[1]], values[which.max(llh)], values[which(cs > .95)[1]])
}