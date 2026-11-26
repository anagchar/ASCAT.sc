\name{run_sc_sequencing}
\alias{run_sc_sequencing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main function to run on shallow-coverage/single-cell (sc) sequencing data.
}
\description{
This function runs ASCAT.sc on single-cell or shallow-coverage
whole-genome sequencing data.
It bins the genome into small bins and keeps the ones with
high consistency of depth of coverage in diploid samples. It counts
(high mapping quality and properly mapped) reads/pairs of reads within
those bins then applies the ASCAT methodology for smoothing (GC
correction with loess), segmenting (CBS or multipcf for multi-region sampling) and fitting of purity/ploidy and integer
copy-number.
It will also automatically refit the profiles to match the feature space
of high-quality consensus PCAWG calls.
Finally it prints the results and profiles to disk and returns a list with all raw and
processed data.
The output format is meant to keep all important info, such that
subtasks of the function (for example segmentation by providing a
different segmentation_alpha) can be re-run on the same object without
re-running the full pipeline.
}
\usage{
run_sc_sequencing(tumour_bams,
                              allchr=paste0("",c(1:22,"X")),
                              sex=c("female","male","female"),
                              chrstring_bam="",
                              purs = seq(0.1, 1, 0.01),
                              ploidies = seq(1.7,5, 0.01),
                              maxtumourpsi=5,
                              binsize=500000,
                              segmentation_alpha=0.01,
                              predict_refit=TRUE,
                              print_results=TRUE,
                              build=c("hg19","hg38"),
                              MC.CORES=1,
                              outdir="./",
                              probs_filters=.1,
                              sc_filters=FALSE,
                              projectname="project",
                              multipcf=FALSE,
                              hp_tag=NULL,
                              hp_allele_specific=TRUE,
                              hp_min_reads=20,
                              hp_baf_tol=0.05,
                              hp_logr_tol=0.5,
                              hp_max_cn=6,
                              hp_coseg=TRUE,
                              hp_coseg_llr=10,
                              hp_coseg_min_bins=5,
                              hp_coseg_method="circular")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tumour_bams}{
    vector of characters. paths to the bams (sorted and indexed).
  }
  \item{res}{
    result object which was output by the function run_sc_sequencing of
    ASCAT.sc. default is NULL for first run.
  }
  \item{allchr}{
    vector of characters. chromosome names, must be in c(1:22,"X", "Y"),
    can have a "chr" prefix. typically would be "chr" for hg38 and "" for hg19.
  }
  \item{sex}{
    vector of characters. sex in "male" or "female" for each of the
    patients with tumours in tumour_bams.
  }
  \item{chrstring_bam}{
    character string. either "" or "chr" corresponding to presence of string in
    contig names in the bam.
  }
  \item{purs}{
    numeric vector. purity grid values. default suggested for bulk. for
    single cells, try: purs=c(0.5,seq(.99,1,.01)).
  }
  \item{ploidies}{
    numeric vector. ploidy grid values. default from 1.7 to 5 can be
    adjusted depending on tumour type. steps of 0.01 can be increased
    for faster runs.
  }
  \item{maxtumourpsi}{
    numeric value. maximum tumour ploidy above which solutions will be
    masked (distance set to infinity in the grid search). helps reduce
    runtime, as some combinations of psi and purs would lead to
    tumour ploidy >> maxtumourpsi.
  }
  \item{binsize}{
    the size of the bins, should be above 30kb and will be set to a multiple of 30kb, which is
    the current smallest bin size. precomputed 5kb bins are available on demand.
  }
  \item{segmentation_alpha}{
    segmentation parameter for the circular binary
    segmentation. penalty=1/segmentation_alpha when using multipcf=T
  }
  \item{predict_refit}{
    boolean. use xgboost predictor to predict refitted ploidy and refit the
    profiles automatically.
  }
  \item{print_results}{
    boolean. print results to disk to directory "outdir".
  }
  \item{build}{
    character string. either "hg19" or "hg38".
  }
  \item{MC.CORES}{
    numeric value. number of cores to use, relies on parallel's mclapply
  }
  \item{barcodes_10x}{
    vector of character strings. legacy for 10X scDNA-seq bam files:
    barcodes for GEL beads associated with cells/nuclei. typically
    located in the "outs/barcodes.csv" of the 10x cell-ranger pipeline.
  }
  \item{outdir}{
    character string. path to output directory.
  }
  \item{probs_filters}{
      numeric value. if sc_filters is TRUE, proportion of bad cells
      expected. leave at .1, more will be identified and removed anyways.
    }
    \item{sc_filters}{
      boolean. for single-cell data, should bad cells be identified and
      flagged. will generate a pdf with quality checks.
    }
  \item{projectname}{
    character string to use as flag in output file names.
  }
  \item{smooth_sc}{
    boolean. for single-cell data and allele-specific CNA calling
    only. Smooths the inferred integer values associated to the BAF/logR
    of each segment by MLE-like fitting per segment across cells.
  }
  \item{multipcf}{
    boolean. use multi track segmentation when bams come from regions of
    the same tumour (i.e. breakpoints are likely shared).
}
  \item{hp_tag}{
    character string or NULL. for haplotagged 10X-style BAMs: name of the
    haplotype tag (e.g. "HP"). reads tagged HP:i:1 / HP:i:2 are counted
    per cell per bin (columns hap1 / hap2 of each track).
  }
  \item{hp_allele_specific}{
    boolean. with hp_tag set, call allele-specific copy number from the
    hap1 / hap2 counts (getAS_CNA_HP) into res$allProfiles_AS.
  }
  \item{hp_min_reads}{
    numeric value. minimum number of haplotagged reads in a segment (and in
    each piece of a BAF split) to fit it.
  }
  \item{hp_baf_tol}{
    numeric value. a logR-compatible state is accepted if its expected
    major-haplotype fraction is within hp_baf_tol of the segment's fitted
    fraction (or its likelihood is close to the best state's).
  }
  \item{hp_logr_tol}{
    numeric value. a state is logR-compatible if its total copy number is
    within hp_logr_tol of the logR copy number of the segment.
  }
  \item{hp_max_cn}{
    numeric value. largest total copy number considered for the
    allele-specific states.
  }
  \item{hp_coseg}{
    boolean. after fitting the tag error and overdispersion of a cell on
    its logR segments, search each logR segment for BAF change-points and
    call the pieces separately (finds copy-neutral LOH inside a logR-flat
    segment). pieces have seg_source "BAF" in res$allProfiles_AS.
  }
  \item{hp_coseg_llr}{
    numeric value. minimum log-likelihood gain to accept a BAF
    change-point. 10 gave about 0.02 false breakpoints per cell in
    simulations on COLO829BL single cells.
  }
  \item{hp_coseg_min_bins}{
    numeric value. minimum number of bins with haplotagged reads in each
    piece of a BAF split.
  }
  \item{hp_coseg_method}{
    character string. "circular" (default): best middle piece against the
    rest of the segment, as circular binary segmentation; "binary": single
    cuts only (misses changes in the middle of a segment).
  }
}
\details{
  
}
\value{

  
}
\references{

}
\author{
maxime.tarabichi@ulb.be
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
run_methylation_data; run_target_sequencing;
}
\examples{


## Please see vignette for more examples
## Setup for CNA calling from bulk shallow coverage WGS in ctDNA (low
## purity).

res <- run_sc_sequencing(tumour_bams=bams,
                         allchr=paste0("chr",c(1:22,"X")),
                         sex=rep("male",length(bams)),
                         binsize=500000,
                         chrstring_bam="chr",
                         purs = seq(0.01, .2, 0.001), # try to force the purity to be low; decrease step size
                         ploidies = seq(1.7,5, 0.01),
                         maxtumourpsi=5,
                         build="hg38",
                         MC.CORES=1,
                         outdir="./",
                         projectname="ctDNA_test",
                         segmentation_alpha=0.01,
                         predict_refit=TRUE,
                         multipcf=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
