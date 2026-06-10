## #####################################################################
## Unified entry point. Thin dispatcher that infers which pipeline to run
## from the supplied inputs and forwards to the existing run_* function.
## The three pipeline functions remain available for direct use.
## #####################################################################

run_ascat <- function(..., type = c("auto", "sc", "targeted", "methylation"))
{
    type <- match.arg(type)
    args <- list(...)
    if (type == "auto")
    {
        type <- if (!is.null(args$idat_dir) || !is.null(args$rgSet)) "methylation"
                else if (!is.null(args$bed_file))                    "targeted"
                else if (!is.null(args$tumour_bams))                 "sc"
                else stop("run_ascat(): cannot infer the pipeline type. ",
                          "Pass type= explicitly, or provide one of ",
                          "tumour_bams (sc), bed_file (targeted) or idat_dir (methylation).")
    }
    fn <- switch(type,
                 sc          = run_sc_sequencing,
                 targeted    = run_targeted_sequencing,
                 methylation = run_methylation_array)
    message(sprintf("run_ascat: dispatching to '%s' pipeline", type))
    do.call(fn, args)
}
