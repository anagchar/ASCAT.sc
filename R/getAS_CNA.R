getAS_CNA <- function(res,
                      path_to_phases,
                      list_ac_counts_paths,
                      purs,
                      ploidies,
                      distances=NULL,
                      chrstring="chr",
                      projectname="project",
                      outdir="./",
                      mc.cores=1,
                      steps=NULL)
{
    suppressPackageStartupMessages(require(GenomicRanges))
    suppressPackageStartupMessages(require(data.table))
    getNANBfromTot <- function(baf,tot,purity, retNa=TRUE)
    {
        Na <- round(baf*(purity*tot+(1-purity)*2)-1+purity)
        Nb <- round(tot-Na)
        if(retNa) return(Na)
        return(Nb)
    }
    .searchGrid  <-  function (baf,
                              logr,
                              sizes,
                              purs=seq(.90,.99,.01),
                              ploidies=seq(2,6,.01))
    {
        getNANB <- function(baf,logr,purity,ploidy)
        {
            K <- (2^logr*ploidy-2*(1-purity))/purity
            Na <- (1-purity-baf*(2-2*purity)-(purity*baf-purity)*K)/purity
            Nb <- K-Na
            list(Na=Na,Nb=Nb)
        }

        errors <- function(baf, logr, sizes, purity, ploidy)
        {
            nanb <- getNANB(baf,logr,purity,ploidy)
            ssize <- sum(sizes)
            sum(((nanb$Na-round(nanb$Na))^2+(nanb$Nb-round(nanb$Nb))^2)*(sizes/ssize))
        }

        ggprofile <- function(baf, logr, purity, ploidy)
        {
            nanb <- getNANB(baf,logr,purity,ploidy)
        }
        errs <- t(sapply(purs,function(rho)
            sapply(ploidies,function(psi)
                errors(baf=baf,
                       logr=logr,
                       sizes,
                       purity=rho,ploidy=psi))))
        rownames(errs) <- purs
        colnames(errs) <- ploidies
        purs <- as.numeric(rownames(errs))
        ploidies <- as.numeric(colnames(errs))
        mins <- arrayInd(which.min(errs), dim(errs))
        purity <- purs[mins[1]]
        ploidy <- ploidies[mins[2]]
        return(list(purity=purity,
                    ploidy=ploidy,
                    profile=ggprofile(baf,logr,purity,ploidy),
                    errs=errs))
    }

    readPhases <- function(phasing_paths)
    {
        phasing <- lapply(phasing_paths,function(x) as.data.frame(data.table::fread(x)))
        phases <- lapply(phasing,function(x)
        {
            ##x <- x[x[,10]%in%c("0|1","1|0"),]
            x <- x[grep("0\\|1|1\\|0",x[, 10]), ]
            phase <- gsub("(.*)\\|(.*)","\\1",x[,10])
            phases1 <- x[,"REF"]
            phases1[phase=="1"] <- x[phase=="1","ALT"]
            phases2 <- x[,"ALT"]
            phases2[phases2==phases1] <- x[phases2==phases1,"REF"]
            list(chr=x[,1],
                 pos=x[,2],
                 phases1=phases1,
                 phases2=phases2)
        })
        phasesall <- do.call("rbind",lapply(phases,function(x) data.frame(chr=x[[1]],
                                                                          pos=x[[2]],
                                                                          phase1=x[[3]],
                                                                          phase2=x[[4]])))
    }

    getPhasedInfo <- function(ac, phasing)
    {
        ac[,1] <- gsub("chr","",ac[,1])
        phasing[,1] <- gsub("chr","",phasing[,1])
        ac.. <- merge(ac, phasing,by.x=c("#CHR","POS"), by.y=c("chr","pos"))
        nmsPH1 <- ac..[,"phase1"]
        nmsPH2 <- ac..[,"phase2"]
        letters <- c("A","C","G","T")
        letters_index <- sapply(paste0("Count_",letters),function(x) which(colnames(ac..)%in%x))
        names(letters_index) <-letters
        ind1 <- letters_index[nmsPH1]
        ind2 <- letters_index[nmsPH2]
        ind1 <- cbind(seq_along(ind1), ind1)
        ind2 <- cbind(seq_along(ind2), ind2)
        df <- data.frame(chr=ac..[,1],
                         pos=ac..[,2],
                         counts1=as.numeric(as.character(ac..[ind1])),
                         counts2=as.numeric(as.character(ac..[ind2])))
        df <- df[rowSums(df[,3:4])>0,]
        df
    }

    getHet <- function(snp1,snp2)
    {
        gr1 <- GRanges(gsub("chr","",snp1[,1]),
                       IRanges(snp1[,2],snp1[,2]))
        gr2 <- GRanges(gsub("chr","",snp2[,1]),
                       IRanges(snp2[,2],snp2[,2]))
        ovs <- findOverlaps(gr1,gr2)
        snp1[queryHits(ovs),]
    }

    getAC <- function(ac_counts_paths, phases)
    {
        acs <- do.call("rbind",lapply(ac_counts_paths,function(dd)
        {
            gc()
            ok <- getHet(as.data.frame(data.table::fread(dd)),phases)
        }))
    }

    fitBinom.1dist <- function(counts, depths, steps=NULL, maxdepth=1000)
    {
        if(is.null(steps))
            steps <- if(length(counts)%/%3>10) 5 else 3
        nonas <- !is.na(counts) & !is.na(depths)
        counts <- counts[nonas]
        depths <- depths[nonas]
        haploblocks <- cut(1:length(counts),pmax(length(counts)%/%steps,2))
        if(length(counts)<10) haploblocks <- rep(1,length(counts))
        lcounts <- split(counts,haploblocks)
        ldepths <- split(depths,haploblocks)
        lcounts <- lapply(1:length(lcounts), function(x) if(rnorm(1)<0) ldepths[[x]]-lcounts[[x]] else lcounts[[x]])
        counts <- sapply(lcounts,sum)
        depths <- sapply(ldepths,sum)
        counts[depths>maxdepth] <- round(counts[depths>maxdepth]/depths[depths>maxdepth]*maxdepth)
        depths[depths>maxdepth] <- maxdepth
        values <- seq(.5,1,0.001)
        llh <- sapply(values,function(x)
        {
            sum(log(dbinom(counts,size=depths,prob=x, log=F)+dbinom(counts,size=depths,prob=1-x, log=F)))
        })
        llh <- llh-max(llh)
        normalised <- exp(llh)/sum(exp(llh))
        ##baf <- sum(values*normalised)
        baf <- values[which.max(llh)]
        cs <- cumsum(normalised)
        q95 <- values[c(which(cs>.05)[1],which(cs>.95)[1])]
        c(q95[1],baf,q95[2])
    }

    fitBinom.1dist.noswitch <- function(counts, depths)
    {
        nonas <- !is.na(counts) & !is.na(depths)
        counts <- counts[nonas]
        depths <- depths[nonas]
        values <- seq(0,1,0.001)
        llh <- sapply(values,function(x)
        {
            sum(dbinom(counts,size=depths,prob=x, log=T))
        })
        llh <- llh-max(llh)
        normalised <- exp(llh)/sum(exp(llh))
        ##baf <- sum(values*normalised)
        baf <- values[which.max(llh)]
        cs <- cumsum(normalised)
        q95 <- values[c(which(cs>.05)[1],which(cs>.95)[1])]
        c(q95[1],baf,q95[2])
    }
    is_distant_enough <- function(positions, distance=1000)
    {
        ##ord <- order(positions,decreasing=F)
        ##psort <- positions[ord]
        ##keep <- c(T,diff(psort)>distance)
        ##return(keep[order(ord,decreasing=F)])
        spositions <- sort(positions)
        keep <- numeric(0)
        prev_kept <- -Inf
        for (pos in spositions)
        {
            if (pos - prev_kept >= distance)
            {
                keep <- c(keep, pos)
                prev_kept <- pos
            }
        }
        positions%in%keep
    }
    getProfile <- function(df,
                           prof,
                           purity,
                           ploidy,
                           purs,
                           ploidies,
                           steps=NULL,
                           distance=1000)
    {
        nprof <- data.frame(chr=as.character(prof[,"chromosome"]),
                            startpos=as.numeric(prof[,"start"]),
                            endpos=as.numeric(prof[,"end"]),
                            total_copy_number=as.numeric(prof[,"total_copy_number"]),
                            total_copy_number_logr=as.numeric(prof[,"total_copy_number_logr"]),
                            logr=as.numeric(prof[,"logr"]),
                            logr.sd=as.numeric(prof[,"logr.sd"]),
                            fitted=as.numeric(prof[,"total_copy_number"]),
                            q05=as.numeric(NA),
                            BAF=as.numeric(NA),
                            q95=as.numeric(NA),
                            q05_noswitch=as.numeric(NA),
                            BAF_noswitch=as.numeric(NA),
                            q95_noswitch=as.numeric(NA),
                            stringsAsFactors=F)
        nprof[nprof[,1]=="chr23",1] <- "chrX"
        nprof[nprof[,1]=="23",1] <- "X"
        nprof[nprof[,1]=="chr24",1] <- "chrY"
        nprof[nprof[,1]=="24",1] <- "Y"
        grseg <- GRanges(gsub("chr","",nprof[,"chr"]),
                         IRanges(as.numeric(nprof[,"startpos"]),
                                 as.numeric(nprof[,"endpos"])))
        grsnp <- GRanges(gsub("chr","",df[,1]),
                         IRanges(df[,2],
                                 df[,2]))
        ovs <- findOverlaps(grseg,
                            grsnp)
        qH <- queryHits(ovs)
        sH <- subjectHits(ovs)
        df[,3] <- as.numeric(as.character(df[,3]))
        df[,4] <- as.numeric(as.character(df[,4]))
        for(i in unique(qH))
        {
            inds <- 1:nrow(df)%in%sH[qH==i]
            inds[inds] <- inds[inds] & is_distant_enough(df[inds,2], distance=distance)
            if(sum(inds)>1)
            {
                nprof[i,c("q05","BAF","q95")] <- fitBinom.1dist(df[inds, 3],
                                                                rowSums(df[inds, c(3,4)]),
                                                                steps=steps)
                nprof[i,c("q05_noswitch","BAF_noswitch","q95_noswitch")] <- fitBinom.1dist.noswitch(df[inds, 3],
                                                                                                    rowSums(df[inds, c(3,4)]))
            }
        }
        nona <- !is.na(nprof[,"BAF"])
        sG <- .searchGrid(as.numeric(nprof[nona,"BAF"]),
                         as.numeric(nprof[nona,"logr"]),
                         as.numeric(nprof[nona,"endpos"])-as.numeric(nprof[nona,"startpos"]),
                         purs=purs,
                         ploidies=ploidies)
        sG.fixed <- .searchGrid(nprof[nona,"BAF"],
                               nprof[nona,"logr"],
                               nprof[nona,"endpos"]-nprof[nona,"startpos"],
                               purs=purity,
                               ploidies=ploidy)
        exceptNA <- function(vec,nona)
        {
            newvec <- rep(NA,length(nona))
            newvec[nona] <- vec
            newvec
        }
        nprof[,"ntot_free"] <- transform_bulk2tumour(nprof[,"logr"],sG$purity,sG$ploidy)
        nprof[,"ntot_fixed"] <- transform_bulk2tumour(nprof[,"logr"],sG.fixed$purity,sG.fixed$ploidy)
        list(nprof.free=cbind(nprof,
                              nA=getNANBfromTot(baf=nprof[,"BAF"], tot=nprof[,"fitted"], purity=purity),
                              nB=getNANBfromTot(baf=nprof[,"BAF"], tot=nprof[,"fitted"], purity=purity, retNa=FALSE),
                              nA_sc_raw=exceptNA(sG$profile$Nb,nona),
                              nB_sc_raw=exceptNA(sG$profile$Na,nona),
                              nA_sc=exceptNA(round(sG$profile$Nb),nona),
                              nB_sc=exceptNA(round(sG$profile$Na),nona)),
             purity.free=sG$purity,
             ploidy.free=sG$ploidy,
             nprof.fixed=cbind(nprof,
                               nA=getNANBfromTot(baf=nprof[,"BAF"], tot=nprof[,"fitted"], purity=purity),
                               nB=getNANBfromTot(baf=nprof[,"BAF"], tot=nprof[,"fitted"], purity=purity, retNa=FALSE),
                               nA_sc_raw=exceptNA(sG.fixed$profile$Nb,nona),
                               nB_sc_raw=exceptNA(sG.fixed$profile$Na,nona),
                               nA_sc=exceptNA(round(sG.fixed$profile$Nb),nona),
                               nB_sc=exceptNA(round(sG.fixed$profile$Na),nona)),
             purity.fixed=sG.fixed$purity,
             ploidy.fixed=sG.fixed$ploidy)
        }

    getBinBAF <- function(ac.ph, track = NULL, cell_name = NULL,
                     distance  = 1000,
                     min_snps  = 20,
                     min_reads = 30,
                     steps     = NULL) {
    suppressPackageStartupMessages({
        require(GenomicRanges)
        require(data.table)
    })

    # Sort SNPs per chromosome by position
    df <- data.table(
        chr = gsub("chr", "", as.character(ac.ph[, 1])),
        pos = as.integer(ac.ph[, 2]),
        c1  = as.numeric(ac.ph[, 3]),
        c2  = as.numeric(ac.ph[, 4])
    )
    df[, dp := c1 + c2]
    setorder(df, chr, pos)

    bin_list <- list()

    chr_levels <- c(as.character(1:22), "X", "Y")
    chrs_present <- intersect(chr_levels, unique(df$chr))

    for (ch in chrs_present) {
        sub <- df[chr == ch & dp > 0]
        if (nrow(sub) == 0) next

        n <- nrow(sub)

        i <- 1
        while (i <= n) {
            # Start a new bin at position i
            bin_start <- sub$pos[i]
            last_kept_pos <- -Inf
            snp_idx <- integer(0)
            total_reads <- 0L
            n_snps <- 0L

            j <- i
            while (j <= n) {
                p <- sub$pos[j]
                # Distance-thinning within the bin
                if (p - last_kept_pos >= distance) {
                    snp_idx <- c(snp_idx, j)
                    n_snps <- n_snps + 1L
                    total_reads <- total_reads + sub$dp[j]
                    last_kept_pos <- p
                    # Check thresholds
                    if (n_snps >= min_snps && total_reads >= min_reads) {
                        break
                    }
                }
                j <- j + 1
            }

            # End-of-chromosome: bin may not meet thresholds; only keep if it does
            if (n_snps >= min_snps && total_reads >= min_reads) {
                bin_end <- sub$pos[j]
                c1_v <- sub$c1[snp_idx]
                dp_v <- sub$dp[snp_idx]

                fit    <- tryCatch(fitBinom.1dist(c1_v, dp_v, steps = steps),
                                   error = function(e) c(NA, NA, NA))
                fit_ns <- tryCatch(fitBinom.1dist.noswitch(c1_v, dp_v),
                                   error = function(e) c(NA, NA, NA))

                bin_list[[length(bin_list) + 1]] <- data.table(
                    chr           = ch,
                    start         = bin_start,
                    end           = bin_end,
                    nSNP          = n_snps,
                    nSNP_kept     = n_snps,        # same here — built from kept SNPs
                    total_reads   = total_reads,
                    bin_size_mb   = (bin_end - bin_start) / 1e6,
                    q05           = fit[1],
                    BAF           = fit[2],
                    q95           = fit[3],
                    q05_noswitch  = fit_ns[1],
                    BAF_noswitch  = fit_ns[2],
                    q95_noswitch  = fit_ns[3]
                )

                # Move to next SNP after the last kept one
                i <- j + 1
            } else {
                # Couldn't fill a bin at the chromosome end — discard remainder, move on
                break
            }
        }
    }

    if (length(bin_list) == 0) {
        warning("No bins met the thresholds.")
        return(data.table())
    }

    bin_df <- rbindlist(bin_list)
    if (!is.null(cell_name)) bin_df[, cell := cell_name]
    as.data.frame(bin_df)
    }

    getAS_CNA_sample <- function(track,
                                 profile,
                                 ac_counts_paths,
                                 purs,
                                 ploidies,
                                 purity,
                                 ploidy,
                                 phases=NULL,
                                 path_to_phases=NULL,
                                 cell_name=NULL,
                                 steps=NULL,
                                 distance=1000)
    {
        if(is.null(phases))
            phases <- readPhases(path_to_phases)
        ac <- getAC(ac_counts_paths, phases)
        ac.ph <- getPhasedInfo(ac, phases)
        bin_baf <- tryCatch(getBinBAF(ac.ph, track, cell_name=cell_name, distance=distance),
                            error=function(e) { warning("getBinBAF failed: ", conditionMessage(e)); NULL })
        prof <- getProfile(ac.ph,
                           prof=profile,
                           steps=steps,
                           purity=purity,
                           ploidy=ploidy,
                           purs=purs,
                           ploidies=ploidies,
                           distance=distance)
        c(prof, list(bin_baf=bin_baf))
    }


    cell_names <- names(res$allTracks.processed)
    if (is.null(names(list_ac_counts_paths)))
        stop("list_ac_counts_paths must be a named list with cell names matching allTracks.processed. ",
             "Got an unnamed list — cannot guarantee correct cell-to-allele-counts matching.")
    missing_cells <- setdiff(cell_names, names(list_ac_counts_paths))
    if (length(missing_cells) > 0)
        stop("list_ac_counts_paths is missing entries for cells: ",
             paste(missing_cells, collapse=", "))
    extra_cells <- setdiff(names(list_ac_counts_paths), cell_names)
    if (length(extra_cells) > 0)
        warning("list_ac_counts_paths has entries not in allTracks.processed (will be ignored): ",
                paste(extra_cells, collapse=", "))

    if (!is.null(distances)) {
        if (is.null(names(distances)))
            stop("distances must be a named list with cell names matching allTracks.processed. ",
                 "Got an unnamed list — cannot guarantee correct cell-to-distance matching.")
        missing_dist_cells <- setdiff(cell_names, names(distances))
        if (length(missing_dist_cells) > 0)
            stop("distances is missing entries for cells: ",
                 paste(missing_dist_cells, collapse=", "))
    }

    phases <- NULL
    if(length(path_to_phases)==1)
    {
        print("## read-in Phases")
        phases <- readPhases(path_to_phases[[1]])
    }
    print("## derive Allele-specific Profiles")
    res$allProfiles_AS <- parallel::mclapply(setNames(seq_along(cell_names), cell_names), function(x)
    {
        cell <- cell_names[x]
        cat(".")
        getAS_CNA_sample(track=res$allTracks.processed[[cell]],
                         profile=res$allProfiles[[cell]],
                         ac_counts_paths=list_ac_counts_paths[[cell]],
                         phases=phases,
                         purity=if(any(grepl("refitted",names(res)))) res$allSolutions.refitted.auto[[cell]]$purity
                                else res$allSolutions[[cell]]$purity,
                         ploidy=if(any(grepl("refitted",names(res)))) res$allSolutions.refitted.auto[[cell]]$ploidy
                                else res$allSolutions[[cell]]$ploidy,
                         purs=purs[[x]],
                         ploidies=ploidies[[x]],
                         path_to_phases=if(length(path_to_phases)>1) path_to_phases[[x]] else NULL,
                         cell_name=cell,
                         steps=steps,
                         distance=if(!is.null(distances)) distances[[cell]] else 1000)
    },mc.cores=mc.cores)
    print("## write to disk and plot Allele-specific Profiles")
    as_cna_dir <- file.path(outdir, "as_cna_profile")
    dir.create(as_cna_dir, showWarnings=FALSE, recursive=TRUE)
    pdf(paste0(outdir,"/all_as_cna_profiles_",projectname,".pdf"),width=15,height=5)
    tnull <- lapply(1:length(res$allProfiles_AS), function(x)
    {
        try({
            plot_AS_profile(res$allProfiles_AS[[x]]$nprof.fixed)
            title(paste0(names(res$allTracks)[x]," - bam",x) ,cex=.5)
        })
        write.table(res$allProfiles_AS[[x]]$nprof.fixed,
                    sep="\t",col.names=T,row.names=F,quote=F,
                    file=file.path(as_cna_dir, paste0("as_cna_profile_",names(res$allTracks)[x],"_bam",x,".txt")))
    })
    dev.off()
    res
}
