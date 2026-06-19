#' @title Normalize logR by cluster pseudobulk
#' @description For each cluster, computes a per-bin pseudobulk logR (median across
#'   cells in the cluster), then subtracts it from each cell's logR. This removes
#'   shared accessibility-driven biases within each cluster while preserving
#'   cell-specific copy-number deviations. Tracks are re-segmented after normalization.
#' @param res Result object from run_sc_sequencing (must contain allTracks.processed)
#' @param clusters Named character vector: barcode -> cluster_id (from getCellGroups)
#' @param segmentation_alpha Segmentation sensitivity for re-segmentation (default: 0.01)
#' @param MC.CORES Number of cores for parallel re-segmentation (default: 1)
#' @return Updated res with normalized smoothed logR, re-segmented tracks, and
#'   \code{$cluster_pseudobulks} storing the per-cluster median logR per bin.
#' @export
normalizeByCluster <- function(res, clusters,
                                segmentation_alpha = 0.01,
                                MC.CORES = 1) {

  require(DNAcopy)
  require(parallel)

  cell_names <- names(res$allTracks.processed)
  allchr <- names(res$allTracks.processed[[1]]$lCTS)
  cluster_ids <- unique(clusters)

  matched <- intersect(names(clusters), cell_names)
  if (length(matched) == 0) {
    bc_to_track <- sapply(names(clusters), function(bc) {
      idx <- grep(paste0(bc, "$"), cell_names)
      if (length(idx) == 1) cell_names[idx] else NA
    })
    bc_to_track <- bc_to_track[!is.na(bc_to_track)]
    if (length(bc_to_track) == 0)
      stop("No overlap between cluster barcodes and allTracks cell names.\n",
           "  Cluster barcodes: ", paste(head(names(clusters), 3), collapse = ", "), "\n",
           "  Track names: ", paste(head(cell_names, 3), collapse = ", "))
    new_clusters <- clusters[names(bc_to_track)]
    names(new_clusters) <- bc_to_track
    clusters <- new_clusters
    matched <- names(clusters)
  }

  message("Normalizing ", length(matched), " cells across ",
          length(unique(clusters[matched])), " clusters")

  nperms <- 10000
  max.ones <- floor(nperms * segmentation_alpha) + 1
  SBDRY <- DNAcopy::getbdry(eta = 0.05, nperm = nperms, max.ones = max.ones)

  pseudobulks <- list()

  for (cl in cluster_ids) {
    cells_cl <- names(clusters[clusters == cl & names(clusters) %in% matched])
    cell_idx <- match(cells_cl, cell_names)

    if (length(cell_idx) < 2) {
      message("  Cluster ", cl, ": ", length(cell_idx), " cell(s), skipping")
      next
    }

    message("  Cluster ", cl, ": ", length(cell_idx), " cells")

    pb_chr <- list()
    for (chr in allchr) {
      n_bins <- length(res$allTracks.processed[[cell_idx[1]]]$lCTS[[chr]]$smoothed)

      mat <- vapply(cell_idx, function(i) {
        res$allTracks.processed[[i]]$lCTS[[chr]]$smoothed
      }, numeric(n_bins))

      pseudobulk <- apply(mat, 1, median, na.rm = TRUE)
      pb_chr[[chr]] <- pseudobulk

      for (i in cell_idx) {
        res$allTracks.processed[[i]]$lCTS[[chr]]$smoothed <-
          res$allTracks.processed[[i]]$lCTS[[chr]]$smoothed - pseudobulk
      }
    }
    pseudobulks[[as.character(cl)]] <- pb_chr

    # Per-cell centering: shift each cell's genome-wide baseline to zero
    for (i in cell_idx) {
      all_vals <- unlist(lapply(allchr, function(chr) {
        res$allTracks.processed[[i]]$lCTS[[chr]]$smoothed
      }))
      cell_median <- median(all_vals, na.rm = TRUE)
      for (chr in allchr) {
        res$allTracks.processed[[i]]$lCTS[[chr]]$smoothed <-
          res$allTracks.processed[[i]]$lCTS[[chr]]$smoothed - cell_median
      }
    }
  }

  message("Re-segmenting normalized tracks...")
  normalized_idx <- match(matched, cell_names)

  res$allTracks.processed[normalized_idx] <- mclapply(normalized_idx, function(i) {
    track <- res$allTracks.processed[[i]]
    track$lSegs <- lapply(seq_along(allchr), function(j) {
      segmentTrack(track$lCTS[[j]]$smoothed,
                   chr = paste0(j),
                   sd = 0,
                   starts = track$lCTS[[j]]$start,
                   ends = track$lCTS[[j]]$end,
                   SBDRY = SBDRY,
                   ALPHA = segmentation_alpha)
    })
    names(track$lSegs) <- paste0(seq_along(allchr))
    track
  }, mc.cores = MC.CORES)

  res$cluster_pseudobulks <- pseudobulks
  message("Cluster normalization complete")
  res
}
