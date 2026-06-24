#' @title Get cell groups from multimodal single-cell data
#' @param fragments Path to a fragments.tsv.gz file (with .tbi index alongside)
#' @param rna_counts Sparse matrix of RNA counts (genes x cells), or path to 10x directory (optional)
#' @param barcode_map Path to CSV mapping ATAC<->GEX barcodes (default: bundled file in inst/extdata).
#'   Must have columns "ATAC" and "GEX". Used when both modalities are provided.
#' @param resolution Clustering resolution parameter (default: 0.8)
#' @param n_dims Number of LSI/PCA dimensions to use (default: 30)
#' @return A list with: \code{$clusters} (named character vector of barcode -> cluster label) and \code{$seurat_obj} (full Seurat object)
#' @export
getCellGroups <- function(fragments,
                          rna_counts = NULL,
                          barcode_map = NULL,
                          resolution = 0.8,
                          n_dims = 30) {

  # Step 1: Validate inputs
  if (!file.exists(fragments)) stop("Fragment file not found: ", fragments)
  if (!file.exists(paste0(fragments, ".tbi"))) stop("Fragment index (.tbi) not found: ", paste0(fragments, ".tbi"))

  # Validate fragment file format (chr, start, end, barcode, count)
  frag_head <- read.table(fragments, nrows = 5, sep = "\t", comment.char = "#")
  if (ncol(frag_head) < 5) {
    stop("Fragment file must be tab-separated with 5 columns: chr, start, end, barcode, count\n",
         "  Found ", ncol(frag_head), " columns in: ", fragments)
  }
  if (!is.integer(frag_head[[2]]) || !is.integer(frag_head[[3]])) {
    stop("Fragment file columns 2-3 (start, end) must be integers\n",
         "  Expected: chr1\\t10000\\t10449\\tACGTACGT\\t1")
  }

  # Validate rna_counts (if provided)
  if (!is.null(rna_counts)) {
    if (is.character(rna_counts)) {
      if (!dir.exists(rna_counts)) stop("RNA counts directory not found: ", rna_counts)
    } else if (!inherits(rna_counts, c("matrix", "dgCMatrix"))) {
      stop("rna_counts must be a sparse matrix (dgCMatrix), matrix, or path to 10x directory")
    }
  }

  # Validate numeric parameters
  if (!is.numeric(resolution) || resolution <= 0) stop("resolution must be a positive number")
  if (!is.numeric(n_dims) || n_dims < 2) stop("n_dims must be >= 2")

  # Step 2: Load barcode map (if RNA provided)
  if (is.null(barcode_map) && !is.null(rna_counts)) {
    barcode_map <- system.file("extdata", "barcodes_atac_gex.csv", package = "ASCAT.sc")
    if (nchar(barcode_map) == 0)
      stop("Barcode map (ATAC<->GEX) not found in installed ASCAT.sc package. ",
           "Pass it explicitly via the barcode_map argument.")
  }

  # Step 3: Build ATAC peak matrix + LSI

  # Count fragments per cell and filter low-quality cells
  total_counts <- Signac::CountFragments(fragments)
  barcodes <- total_counts[total_counts$frequency_count > 1000, ]$CB

  # Create fragment object with filtered barcodes
  frag_obj <- Signac::CreateFragmentObject(path = fragments, cells = barcodes)

  # Call peaks from fragments — find macs2 or macs3 on PATH
  macs_bin <- Sys.which("macs2")
  if (macs_bin == "") macs_bin <- Sys.which("macs3")
  if (macs_bin == "")
    stop("Neither macs2 nor macs3 found on PATH. Install MACS: https://macs3-project.github.io/MACS/")
  peaks_gr <- Signac::CallPeaks(frag_obj, macs2.path = macs_bin)

  # Build peak x cell count matrix
  peak_matrix <- Signac::FeatureMatrix(fragments = frag_obj, features = peaks_gr, cells = barcodes)

  # Create Seurat object with ChromatinAssay
  chrom_assay <- Signac::CreateChromatinAssay(
    counts = peak_matrix,
    fragments = frag_obj,
    min.cells = 10
  )
  seurat_obj <- Seurat::CreateSeuratObject(counts = chrom_assay, assay = "ATAC")

  # LSI dimensionality reduction
  seurat_obj <- Signac::RunTFIDF(seurat_obj)
  seurat_obj <- Signac::FindTopFeatures(seurat_obj, min.cutoff = "q0")
  seurat_obj <- Signac::RunSVD(seurat_obj)

  # Step 4: Build RNA PCA (if rna_counts provided)
  if (!is.null(rna_counts)) {
    # Load RNA counts
    if (is.character(rna_counts)) {
      rna_matrix <- Seurat::Read10X(data.dir = rna_counts)
    } else {
      rna_matrix <- rna_counts
    }
    
    # Map GEX barcodes to ATAC barcodes
    bc_map <- read.csv(barcode_map)
    atac_cells <- colnames(seurat_obj)
    bc_map <- bc_map[bc_map$ATAC %in% atac_cells, ]
    rna_matrix <- rna_matrix[, colnames(rna_matrix) %in% bc_map$GEX]
    
    # Rename RNA columns to ATAC barcodes
    idx <- match(colnames(rna_matrix), bc_map$GEX)
    colnames(rna_matrix) <- bc_map$ATAC[idx]

    # Subset to cells present in both modalities
    shared_cells <- intersect(colnames(seurat_obj), colnames(rna_matrix))
    seurat_obj <- subset(seurat_obj, cells = shared_cells)
    rna_matrix <- rna_matrix[, shared_cells]
    message("Cells in both ATAC and RNA: ", length(shared_cells))

    # Add RNA assay and run PCA
    seurat_obj[["RNA"]] <- Seurat::CreateAssayObject(counts = rna_matrix)
    Seurat::DefaultAssay(seurat_obj) <- "RNA"
    seurat_obj <- Seurat::NormalizeData(seurat_obj)
    seurat_obj <- Seurat::FindVariableFeatures(seurat_obj)
    seurat_obj <- Seurat::ScaleData(seurat_obj)
    seurat_obj <- Seurat::RunPCA(seurat_obj)
  }

  # Step 5: Clustering (ATAC-only or WNN)
  if (!is.null(rna_counts)) {
    # WNN clustering (ATAC + RNA)
    seurat_obj <- Seurat::FindMultiModalNeighbors(
      seurat_obj,
      reduction.list = list("lsi", "pca"),
      dims.list = list(2:n_dims, 1:n_dims)
    )
    seurat_obj <- Seurat::FindClusters(seurat_obj, graph.name = "wsnn",
                                        algorithm = 3, resolution = resolution)
  } else {
    # ATAC-only clustering
    seurat_obj <- Seurat::FindNeighbors(seurat_obj, reduction = "lsi", dims = 2:n_dims)
    seurat_obj <- Seurat::FindClusters(seurat_obj, algorithm = 3, resolution = resolution)
  }

  # Step 6: Return results
  clusters <- setNames(as.character(seurat_obj$seurat_clusters), colnames(seurat_obj))
  method <- if (!is.null(rna_counts)) "wnn" else "lsi"
  message("Clustering complete: ", length(unique(clusters)), " clusters, ",
          length(clusters), " cells (method: ", method, ")")

  list(
    clusters = clusters,
    method = method,
    seurat_obj = seurat_obj
  )
}
