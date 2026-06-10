## #####################################################################
## User-facing output layer for ASCAT.sc result objects ("ascat.sc").
##
## The run_* entry points return a flat list `res`; every internal
## consumer reads its fields directly. To keep that contract intact we do
## NOT restructure the list. Instead we tag it with class "ascat.sc" and
## add a clean print/summary plus accessor functions on top. The
## accessors are the "profiles" view; the full object stays the
## "detailed" view.
## #####################################################################

.ascatsc_mode <- function(res)
{
    if(!is.null(res$mode)) res$mode else "unknown"
}

## Per-sample purity/ploidy table. Extracted from printResults_all() so it
## is available even when print_results=FALSE. `solutions`/`profiles`
## default to the main fit, but can point at the refitted slots.
.ascatsc_summary_table <- function(res,
                                   solutions=res$allSolutions,
                                   profiles=res$allProfiles)
{
    .mytry <- function(x, retVal=NA)
    {
        out <- try(x, silent=TRUE)
        if(inherits(out, "try-error")) return(retVal)
        out
    }
    getploidy <- function(tt)
    {
        tt <- data.frame(chromosome=as.character(tt[,"chromosome"]),
                         start=as.numeric(tt[,"start"]),
                         end=as.numeric(tt[,"end"]),
                         total_copy_number=as.numeric(tt[,"total_copy_number"]))
        sizes <- (tt$end-tt$start)/1000000
        isna <- is.na(sizes) | is.na(tt$total_copy_number)
        sum(tt$total_copy_number[!isna]*sizes[!isna], na.rm=TRUE)/sum(sizes[!isna], na.rm=TRUE)
    }
    samplenames <- if(!is.null(names(res$allTracks))) names(res$allTracks)
                   else names(res$allProfiles)
    data.frame(samplename=samplenames,
               purity=sapply(solutions, function(x) .mytry(x$purity)),
               ploidy=sapply(solutions, function(x) .mytry(x$ploidy)),
               ploidy.tumour=sapply(profiles, function(x) .mytry(getploidy(x))),
               stringsAsFactors=FALSE)
}

print.ascat.sc <- function(x, ...)
{
    try({
        mode <- .ascatsc_mode(x)
        n <- length(x$allProfiles)
        cat("ASCAT.sc result (", mode, ") - ", n,
            " sample", if(n != 1) "s" else "", "\n", sep="")
        meta <- character(0)
        if(!is.null(x$build))   meta <- c(meta, paste0("build: ", x$build[1]))
        if(!is.null(x$binsize)) meta <- c(meta, paste0("binsize: ",
                                                       format(x$binsize, scientific=FALSE, big.mark=",")))
        if(!is.null(x$isPON))   meta <- c(meta, paste0("PON: ", isTRUE(x$isPON)))
        if(length(meta)) cat("  ", paste(meta, collapse="  |  "), "\n", sep="")
        hasAS   <- "allProfiles_AS" %in% names(x)
        hasASsm <- "allProfiles_AS_smoothed" %in% names(x)
        cat("  allele-specific: ", if(hasAS) "yes" else "no",
            if(hasAS) paste0(" (smoothed: ", if(hasASsm) "yes" else "no", ")") else "",
            "\n", sep="")
        tab <- try(.ascatsc_summary_table(x), silent=TRUE)
        if(!inherits(tab, "try-error") && !is.null(tab) && nrow(tab)) {
            cat("\n")
            print(utils::head(tab, 6), row.names=FALSE)
            if(nrow(tab) > 6)
                cat("  ... ", nrow(tab) - 6, " more (use summary(res))\n", sep="")
        }
        cat("\nUse summary(res), getProfiles(res), getProfilesAS(res), getMetadata(res).\n")
    }, silent=TRUE)
    invisible(x)
}

summary.ascat.sc <- function(object, ...)
{
    .ascatsc_summary_table(object)
}

## --- Accessors: the user-facing "profiles" API ----------------------

getProfiles <- function(res)
{
    res$allProfiles
}

getSolutions <- function(res)
{
    res$allSolutions
}

getSummaryTable <- function(res)
{
    .ascatsc_summary_table(res)
}

getProfilesAS <- function(res)
{
    if(!"allProfiles_AS" %in% names(res)) {
        message("No allele-specific profiles in this object (mode: ",
                .ascatsc_mode(res), ").")
        return(NULL)
    }
    res$allProfiles_AS
}

getProfilesASsmoothed <- function(res)
{
    if(!"allProfiles_AS_smoothed" %in% names(res)) {
        message("No smoothed allele-specific profiles in this object.")
        return(NULL)
    }
    res$allProfiles_AS_smoothed
}

getRefitted <- function(res)
{
    out <- list(auto=res$allProfiles.refitted.auto,
                manual=res$allProfiles.refitted.manual)
    out <- out[!sapply(out, is.null)]
    if(length(out) == 0) {
        message("No refitted profiles in this object.")
        return(NULL)
    }
    out
}

getMetadata <- function(res)
{
    fields <- c("mode", "build", "binsize", "sex", "chr", "purs", "ploidies",
                "maxtumourpsi", "isPON", "multipcf", "segmentation_alpha")
    res[fields[fields %in% names(res)]]
}
