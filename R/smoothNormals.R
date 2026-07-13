smoothNormals <- function(logr, lNormals)
{
    notalreadyinpanel <- !sapply(lNormals,function(x) cor(log2(unlist(lapply(x,function(y) y$records)) +1),logr))>.999
    ## Never drop the sole column of a length-1 (LOO) panel: cell i is already
    ## excluded from the pseudobulk, so it is not a self-duplicate even if it
    ## correlates > .999. Multi-column panels are untouched.
    if(length(notalreadyinpanel)==1) notalreadyinpanel[1] <- TRUE
    normals <- sapply(lNormals,function(x) log2(unlist(lapply(x,function(y) y$records)) +1))
    normals <- as.matrix(normals) ## keep a length-1 panel as a 1-column matrix
    cat(".")
    predicted <- lm(y ~ . - 1, data = data.frame(y = 2^logr-1,
                                                 X = 2^normals[, notalreadyinpanel, drop=FALSE])-1)$fitted.values
    predicted[predicted < 1] <- 1
    print(summary(log2(predicted+1)))
    logr <- logr-log2(predicted+1)
    logr
}
