run_sc_sequencing <- function(tumour_bams,
                              res=NULL,
                              allchr=paste0("",c(1:22,"X")),
                              sex=c("female","male","female"),
                              chrstring_bam="",
                              purs = seq(0.1, 1, 0.01),
                              ploidies = seq(1.7,5, 0.01),
                              maxtumourpsi=5,
                              binsize=500000,
                              segmentation_alpha=0.01,
                              pcchromosome = 0.8,
                              predict_refit=TRUE,
                              print_results=TRUE,
                              build=c("hg19","hg38"),
                              MC.CORES=1,
                              barcodes_10x=NULL,
                              normal_bams=NULL,
                              normal_barcodes=NULL,
                              outdir="./",
                              is_pdf=F,
                              probs_filters=.1,
                              path_to_phases=NULL,
                              list_ac_counts_paths=NULL,
                              sc_filters=FALSE,
                              projectname="project",
                              steps=NULL,
                              smooth_sc=FALSE,
                              multipcf=FALSE,
                              normalize=FALSE,
                              lExclude=NULL,
                              svinput=NULL,
                              lSVinput=NULL,
                              sc_exclude_badbins=FALSE,
                              clustering=FALSE,
                              fragments=NULL,
                              rna_counts=NULL,
                              readlen_correction=NULL,
                              readlen_stat=c("median","mean"),
                              readlen_by_cluster=FALSE,
                              readlen_min_reads_per_bin=10)
{
    checkArguments_scs(c(as.list(environment())))
    ## Read-length correction of the coverage on top of GC (see smoothCoverageTrack / getReadLengthTrack):
    ##   readlen_correction: "pseudobulk" = read length per bin pooled across cells, a copy-number-free
    ##                       marker of Tn5 accessibility (default with 10X barcodes); "cell" = each track's
    ##                       own median read length (default for bulk BAMs); "none" = GC only
    ##   readlen_stat:       pooled statistic, "median" (from median_qwidth) or "mean" (nucleotides/records)
    ##   readlen_by_cluster: FALSE = one track pooled over all cells; TRUE = one track per cluster of
    ##                       res$clustering$clusters (clustering=TRUE); or a named vector cell -> group
    ##   readlen_min_reads_per_bin: bins with fewer pooled reads get the genome-wide median read length
    readlen_stat <- match.arg(readlen_stat)
    if(is.null(readlen_correction)) readlen_correction <- if(!is.null(barcodes_10x)) "pseudobulk" else "cell"
    readlen_correction <- match.arg(readlen_correction, c("pseudobulk","cell","none"))
    if(readlen_correction=="pseudobulk" & is.null(barcodes_10x))
    {
        print("## readlen_correction='pseudobulk' needs 10X barcodes - using 'cell' (per-sample read length)")
        readlen_correction <- "cell"
    }

    ## Cell clustering via getCellGroups()
    if (clustering) {
        if (is.null(fragments)) stop("clustering=TRUE requires 'fragments' (path to fragments.tsv.gz)")
        print("## Running getCellGroups() for cell clustering")
        clustering_result <- getCellGroups(fragments = fragments, rna_counts = rna_counts)
        if (is.null(res)) res <- list()
        res$clustering <- list(
            clusters = clustering_result$clusters,
            method = clustering_result$method,
            seurat_obj = clustering_result$seurat_obj
        )
        print(paste0("## Clustering done: ", length(unique(clustering_result$clusters)),
                      " clusters, ", length(clustering_result$clusters),
                      " cells (method: ", clustering_result$method, ")"))
    }

    suppressPackageStartupMessages(require(parallel))
    suppressPackageStartupMessages(require(Rsamtools))
    suppressPackageStartupMessages(require(Biostrings))
    suppressPackageStartupMessages(require(DNAcopy))
    suppressPackageStartupMessages(require(copynumber))
    binsize <- as.numeric(binsize)
    print(paste0("## load Bins for Genome Build: ", build))
    START_WINDOW <- 30000
    data("SBDRYs_precomputed",package="ASCAT.sc")
    if(!as.character(segmentation_alpha)%in%names(SBDRYs))
    {
        nperms <- 10000;
        max.ones <- floor(nperms * segmentation_alpha) + 1
        SBDRY <- DNAcopy::getbdry(eta=0.05, nperm=nperms, max.ones=max.ones)
    } else
    {
        SBDRY <- SBDRYs[[as.character(segmentation_alpha)]]
    }
    if(binsize<30000)
    {
        print("Current minimum bin size is 30000 - resetting to 30000")
    }
    if(is.null(res))
        res <- list()
    if(build=="hg19")
    {
        data("lSe_filtered_30000.hg19",package="ASCAT.sc")
        data("lGCT_filtered_30000.hg19",package="ASCAT.sc")
        allchr. <- gsub("chr","",allchr)
        res$lSe <- lapply(allchr., function(chr) lSe.hg19.filtered[[chr]])
        names(res$lSe) <- allchr
        names(lGCT.hg19.filtered) <- names(res$lSe)
        res$lGCT <- lapply(allchr, function(chr) lGCT.hg19.filtered[[chr]])
    }
    if(build=="hg38")
    {
        data("lSe_filtered_30000.hg38",package="ASCAT.sc")
        data("lGCT_filtered_30000.hg38",package="ASCAT.sc")
        names(lGCT.hg38.filtered) <- names(lSe.hg38.filtered)
        allchr. <- paste0("chr",gsub(chrstring_bam,"",allchr))
        res$lSe <- lapply(allchr., function(chr) lSe.hg38.filtered[[chr]])
        names(res$lSe) <- allchr
        res$lGCT <- lapply(allchr., function(chr) lGCT.hg38.filtered[[chr]])
        names(res$lGCT) <- allchr
        if(chrstring_bam=="")
            names(res$lGCT) <- names(res$lSe) <- gsub("chr","",names(res$lSe))
    }
    if (build == "mm39")
    {
        data("lSe_unfiltered_5000.mm39", package = "ASCAT.sc")
        data("lGCT_unfiltered_5000.mm39", package = "ASCAT.sc")
        names(lGCT)[1:length(allchr)] <- names(lSe)[1:length(allchr)] <- allchr
        res$lSe <- lapply(allchr, function(x) lSe[[x]])
        res$lGCT <- lapply(allchr, function(x) lGCT[[x]])
        names(res$lGCT) <- names(res$lSe) <- allchr
    }
    if (!build %in% c("hg19", "hg38","mm39"))
    {
        START_WINDOW <- median(unlist(lapply(lSe,function(x) x$ends-x$starts+1)))
        res$lSe <- lSe
        res$lGCT <- lGCT
        names(res$lGCT) <- names(res$lSe) <- allchr
    }
    print("## calculate Target Bin size")
    res$nlGCT <- treatGCT(res$lGCT,window=ceiling(binsize/START_WINDOW))
    res$nlSe <- treatlSe(res$lSe,window=ceiling(binsize/START_WINDOW))
    if(is.null(res$isPON)) res$isPON <- FALSE
    if(!any(names(res)=="binsize")) res$binsize <- binsize
    if(!is.null(barcodes_10x))
    {
        if(!any(names(res)=="allTracks"))
        {
            print("## get Tracks from 10X-like bam")
            res$timetoread_tumours <- system.time(res$allTracks <- lapply(tumour_bams,function(bamfile)
            {
                lCTS.tumour <- getTrackForAll.10XBAM(bamfile=bamfile,
                                                     barcodes=barcodes_10x,
                                                     allchr=allchr,
                                                     chrstring="",
                                                     lSe=res$lSe,
                                                     doSeg=FALSE,
                                                     # Percent of chromosome bins that need to have a barcode to be kept
                                                     pcchromosome = pcchromosome,
                                                     mc.cores=MC.CORES)

                # Skip BAM files with no detected barcodes
                if (length(lCTS.tumour) == 0) return(NULL)

                # If multiple bams, rename the tracks to include bamfile (sample)

                if (length(tumour_bams) > 1) {
                    names(lCTS.tumour) <- paste0(basename(bamfile), "_", names(lCTS.tumour))
                }

                # Fix binsizes to match the desired binsize using treatTrack()
                mclapply(lCTS.tumour, function(x) {
                    list(lCTS.tumour = x,
                         nlCTS.tumour = treatTrack(lCTS = x,
                                                   window = ceiling(binsize / START_WINDOW)))
                }, mc.cores = MC.CORES)
            }))
            # Remove NULL entries (BAM files with no detected barcodes)
            res$allTracks <- Filter(Negate(is.null), res$allTracks)
            # Unlist the top layer to flatten list
            res$allTracks <- unlist(res$allTracks, recursive = FALSE)

            # sex[1] because the assumption is when running multiple bamfiles, they are from the same tumour?
            # This ensures that it has the correct length
            sex <- rep(sex[1], length(res$allTracks))
        }
    }
    if(is.null(barcodes_10x))
    {
        if(!is.null(normal_bams[1]) & is.null(res$nlCTS.normal))
        {
            print("## get all tracks from normal bams")
            timetoread_normals <- system.time(res$lCTS.normal <- mclapply(normal_bams,function(bamfile)
            {
                lCTS.normal <- lapply(allchr, function(chr) getCoverageTrack(bamPath=bamfile,
                                                                                                   chr=chr,
                                                                                                   res$lSe[[chr]]$starts,
                                                                                                   res$lSe[[chr]]$ends,
                                                                                                   mapqFilter=30))
                list(lCTS.normal=lCTS.normal,
                     nlCTS.normal=treatTrack(lCTS=lCTS.normal,
                                             window=ceiling(binsize/START_WINDOW)))
            },mc.cores=MC.CORES))
            res$nlCTS.normal <- combineDiploid(lapply(res$lCTS.normal,function(x) x[[2]]))
            res$lNormals <- lapply(res$lCTS.normal,function(x) x$nlCTS.normal)
            res$isPON <- TRUE
        }
        else
        {
            res$isPON <- FALSE
            res$nlCTS.normal <- NULL
            res$lNormals <- NULL
            res$timetoread_normals <- NULL
        }
        if(any(names(res)=="allTracks") & res$binsize!=binsize)
        {
            print("## adjust Tracks for bin size ")
            res$timetoread_tumours <- system.time(res$allTracks <- mclapply(names(res$allTracks),function(bamfile)
            {
                lCTS.tumour <- res$allTracks[[bamfile]]$lCTS.tumour
                list(lCTS.tumour=lCTS.tumour,
                     nlCTS.tumour=treatTrack(lCTS=lCTS.tumour,
                                             window=ceiling(binsize/START_WINDOW)))
            },mc.cores=MC.CORES))
        }
        if(!any(names(res)=="allTracks"))
        {
            print("## get Tracks from Tumour bams")
            res$timetoread_tumours <- system.time(res$allTracks <- mclapply(tumour_bams,function(bamfile)
            {
                lCTS.tumour <- lapply(allchr, function(chr) getCoverageTrack(bamPath=bamfile,
                                                                             chr=chr,
                                                                             res$lSe[[chr]]$starts,
                                                                             res$lSe[[chr]]$ends,
                                                                             mapqFilter=30))
                list(lCTS.tumour=lCTS.tumour,
                     nlCTS.tumour=treatTrack(lCTS=lCTS.tumour,
                                             window=ceiling(binsize/START_WINDOW)))
            },mc.cores=MC.CORES))
            names(res$allTracks) <- basename(tumour_bams)
        }
    }
    ## Create normal reference from diploid barcodes (10X mode)
    if(!is.null(barcodes_10x) && !is.null(normal_barcodes))
    {
        missing_bc <- normal_barcodes[!normal_barcodes %in% names(res$allTracks)]
        if(length(missing_bc) > 0)
            stop(paste0("normal_barcodes not found in allTracks: ",
                        paste(head(missing_bc, 5), collapse=", "),
                        if(length(missing_bc) > 5) paste0(" ... and ", length(missing_bc)-5, " more")))
        print(paste0("## create normal reference from ", length(normal_barcodes), " diploid barcodes"))
        res$lNormals <- lapply(normal_barcodes, function(bc) res$allTracks[[bc]]$nlCTS.tumour)
        res$isPON <- TRUE
        ## Force re-smoothing so the normal reference is applied
        res$allTracks.processed <- NULL
    }
    ## Ensure sex has correct length for 10X barcodes on rerun
    if(!is.null(barcodes_10x) && length(sex) != length(res$allTracks))
        sex <- rep(sex[1], length(res$allTracks))
    if(sc_exclude_badbins)
    {
        print("## exclude Bad bins inferred from normal samples or tumour samples")
        res <- sc_excludeBadBins(res)
    }
    ## Bin-level read length pooled across cells (same bins as res$nlGCT): one track, or one per group
    lRLT_cells <- NULL
    res$nlRLT <- NULL
    readlen_groups <- NULL
    if(readlen_correction=="pseudobulk")
    {
        if(isTRUE(readlen_by_cluster))
        {
            if(is.null(res$clustering$clusters))
                print("## readlen_by_cluster=TRUE but no clustering in res (run with clustering=TRUE) - using one track for all cells")
            else
                readlen_groups <- res$clustering$clusters
        }
        else if(!isFALSE(readlen_by_cluster))
            readlen_groups <- readlen_by_cluster
        if(!is.null(readlen_groups) && !any(names(res$allTracks) %in% names(readlen_groups)))
        {
            print("## readlen_by_cluster: the group names do not match the cell names - using one track for all cells")
            readlen_groups <- NULL
        }
        print(paste0("## build pseudobulk read-length track (", readlen_stat, ")"))
        res$nlRLT <- getReadLengthTrack(lapply(res$allTracks, function(x) x$nlCTS.tumour),
                                        stat=readlen_stat,
                                        groups=readlen_groups,
                                        min_reads_per_bin=readlen_min_reads_per_bin)
        lRLT_cells <- lapply(seq_along(res$allTracks), function(x)
            pickReadLengthTrack(res$nlRLT, if(is.null(names(res$allTracks))) NA else names(res$allTracks)[x]))
    }
    ## Re-smooth when rerunning on a result smoothed with other read-length settings (or before they existed)
    readlen_settings <- list(correction=readlen_correction,
                             stat=if(is.null(res$nlRLT)) NULL else attr(res$nlRLT, "stat"),
                             groups=readlen_groups,
                             min_reads_per_bin=readlen_min_reads_per_bin)
    if(!is.null(res$allTracks.processed) && !identical(res$readlen_settings, readlen_settings))
    {
        print("## read-length correction settings changed - tracks will be re-smoothed")
        res$allTracks.processed <- NULL
    }
    if(multipcf)
    {
        print("## calculate Multipcf - multi-sample mode - do not use if samples from different tumours")
        res$timetoprocessed <- system.time(res$allTracks.processed <- getLSegs.multipcf(allTracks=lapply(res$allTracks, function(x)
        {
            list(lCTS=x$nlCTS.tumour)
        }),
        lCTS=lapply(res$allTracks,function(x) x$nlCTS.tumour),
        lSe=res$nlSe,
        lGCT=res$nlGCT,
        lNormals=res$lNormals,
        allchr=allchr,
        segmentation_alpha=segmentation_alpha,
        normalize=normalize,
        svinput=svinput,
        lRLT=lRLT_cells,
        readlen=readlen_correction,
        MC.CORES=MC.CORES))
    }
    else
    {
        if(!any("allTracks.processed"%in%names(res)) | res$binsize!=binsize)
        {
            print("## smooth Tracks")
            res$timetoprocessed <- system.time(res$allTracks.processed <- mclapply(1:length(res$allTracks), function(x)
            {
                getTrackForAll(bamfile=NULL,
                               window=NULL,
                               lCT=res$allTracks[[x]][[2]],
                               lSe=res$nlSe,
                               lGCT=res$nlGCT,
                               lNormals=res$lNormals,
                               allchr=allchr,
                               sdNormalise=0,
                               SBDRY=SBDRY,
                               svinput=if(is.null(lSVinput)) NULL else lSVinput[[x]],
                               segmentation_alpha=segmentation_alpha,
                               lRLT=if(is.null(lRLT_cells)) NULL else lRLT_cells[[x]],
                               readlen=readlen_correction)
            },mc.cores=MC.CORES))
        }
    }
    if(is.null(barcodes_10x))
        names(res$allTracks.processed) <- names(res$allTracks) <- basename(tumour_bams)
    else
        names(res$allTracks.processed) <- names(res$allTracks)

    if(!is.list(purs)) 
    {
        purs <- lapply(1:length(res$allTracks.processed), function(x) purs)
        ploidies <- lapply(1:length(res$allTracks.processed), function(x) ploidies)
    }
    print("## fit Purity/Ploidy")
    res$timetofit <- system.time(res$allSols <- mclapply(1:length(res$allTracks.processed), function(x)
    {
        sol <- try(searchGrid(res$allTracks.processed[[x]],
                              purs = purs[[x]],
                              ploidies = ploidies[[x]],
                              maxTumourPhi=maxtumourpsi,
                              ismale=if(sex[x]=="male") T else F,
                              isPON=res$isPON),silent=F)
    },mc.cores=MC.CORES))
    print("## get Fitted CNA Profiles")
    res$allProfiles <- mclapply(1:length(res$allTracks.processed), function(x)
    {
        try(getProfile(fitProfile(res$allTracks.processed[[x]],
                                  purity=res$allSols[[x]]$purity,
                                  ploidy=res$allSols[[x]]$ploidy,
                                  ismale=if(sex[x]=="male") T else F),
                       CHRS=allchr),silent=F)
    },mc.cores=MC.CORES)
    names(res$allProfiles) <- names(res$allSols) <- names(res$allTracks)
    print("## compile Results")
    res <- list(allTracks.processed=res$allTracks.processed,
                allTracks=res$allTracks,
                allSolutions=res$allSols,
                allProfiles=res$allProfiles,
                chr=allchr,
                chrstring_bam=chrstring_bam,
                purs=purs,
                ploidies=ploidies,
                maxtumourpsi=maxtumourpsi,
                build=build,
                binsize=binsize,
                sex=sex,
                segmentation_alpha=segmentation_alpha,
                multipcf=multipcf,
                lSe=res$nlSe,
                lGCT=res$nlGCT,
                isPON=res$isPON,
                timetoread_tumours=res$timetoread_tumours,
                timetoprocessed=res$timetoprocessed,
                timetofit=res$timetofit,
                svinput=svinput,
                lSVinput=lSVinput,
                lRLT=res$nlRLT,
                readlen_settings=readlen_settings,
                clustering=res$clustering,
                mode="sc")
    if(predict_refit)
    {
        print("## predict Refit")
        res <- predictRefit_all(res)
    }
    if(print_results)
    {
        print("## print Results")
        res <- printResults_all(res, svinput=svinput, lSVinput=lSVinput, outdir=outdir, projectname=projectname, is_pdf=is_pdf)
    }
    if(sc_filters)
    {
        print("## get Filters")
        res<- getFilters(res,
                         probs=probs_filters,
                         outdir=outdir,
                         projectname=projectname)
    }
    if(!is.null(list_ac_counts_paths) & !is.null(path_to_phases))
    {
        print("## get Allele-specific CNA")
        res <- getAS_CNA(res,
                         path_to_phases=path_to_phases,
                         list_ac_counts_paths=list_ac_counts_paths,
                         purs=purs,
                         ploidies=ploidies,
                         outdir=outdir,
                         projectname=projectname,
                         steps=steps,
                         mc.cores=MC.CORES)
        names(res$allProfiles_AS) <- names(res$allProfiles)
    }
    if(smooth_sc & any(grepl("_AS",names(res))))
    {
        print("## smooth Across Single-Cells/Nulcei")
        res <- getAS_CNA_smoothed(res,
                                  mc.cores=MC.CORES)
    }
    if(smooth_sc & !any(grepl("_AS",names(res))))
        print("Warning: Smoothing is only possible for allele-specific copy numbers")
    res
}
