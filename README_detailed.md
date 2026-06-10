# ASCAT.sc — Detailed Function Reference

**ASCAT.sc** (Allele-Specific Copy number Analysis of Tumours — single cell) is an R package for LogR-based copy number calling from:
- Single-cell DNA sequencing (shallow or deep)
- Shallow-coverage bulk sequencing
- Targeted/exome sequencing (off-target reads)
- Illumina methylation arrays (450K, EPICv1, EPICv2)

It produces per-cell (or per-sample) total and allele-specific copy number profiles, purity, and ploidy estimates.

**Version:** 0.1 | **Author:** Maxime Tarabichi | **License:** GNU-GPL V3

---

## Table of Contents

1. [Pipeline Overview](#pipeline-overview)
2. [Main Entry-Point Functions](#main-entry-point-functions)
3. [Coverage Extraction](#coverage-extraction)
4. [Binning and Genomic Coordinate Functions](#binning-and-genomic-coordinate-functions)
5. [GC Correction and Smoothing](#gc-correction-and-smoothing)
6. [Segmentation](#segmentation)
7. [Purity and Ploidy Estimation](#purity-and-ploidy-estimation)
8. [Copy Number Profile Generation](#copy-number-profile-generation)
9. [Allele-Specific Copy Numbers](#allele-specific-copy-numbers)
10. [Quality Control and Filtering](#quality-control-and-filtering)
11. [Refit and Refinement](#refit-and-refinement)
12. [Methylation Array Functions](#methylation-array-functions)
13. [Targeted Sequencing Functions](#targeted-sequencing-functions)
14. [Visualization and Output](#visualization-and-output)
15. [Utility and Helper Functions](#utility-and-helper-functions)
16. [Key Data Structures](#key-data-structures)
17. [Dependencies](#dependencies)

---

## Pipeline Overview

```
Input (BAM / IDAT / BED)
        │
        ▼
Bin Definition (lSe)  ←  getStartsEnds / getBins_lSe_GC_from_fasta
        │
        ▼
Coverage Extraction  ←  getCoverageTrack / getCoverageTrack.10XBAM
        │
        ▼
GC Correction & Smoothing  ←  smoothCoverageTrack / treatTrack / treatlSe
        │
        ▼
Normal Reference Subtraction  ←  combineDiploid / smoothNormals
        │
        ▼
Segmentation (CBS / PCF / multipcf)  ←  segmentTrack / getLSegs.multipcf
        │
        ▼
Normalisation  ←  normaliseByPloidy
        │
        ▼
Grid Search: Purity × Ploidy  ←  searchGrid / buildDistanceMatrix / geterrors
        │
        ▼
Profile Fitting  ←  fitProfile / getProfile
        │
        ▼
Solution Validation / ML Refit  ←  predictRefit_all / isGoodSolution
        │
        ▼
Allele-Specific CNA (optional)  ←  getAS_CNA / getAS_CNA_smoothed
        │
        ▼
Output & Visualisation  ←  printResults_all / plotSolution / sc_plotHeat
```

---

## Main Entry-Point Functions

> **Public API.** ASCAT.sc exports ~40 user-facing functions (entry points, result accessors, plotting, and the documented low-level building blocks). All other functions are internal but remain callable with `ASCAT.sc:::fn()`.

### `run_ascat()` — `run_ascat.R`

**Recommended single entry point.** Thin dispatcher that auto-detects the pipeline from its inputs and forwards to `run_sc_sequencing()` / `run_targeted_sequencing()` / `run_methylation_array()`:

- `idat_dir` or `rgSet` present → methylation
- `bed_file` present → targeted
- otherwise `tumour_bams` → single-cell / shallow-coverage

Override with `type = "sc" | "targeted" | "methylation"`. Arguments are forwarded by name, so use the same arguments documented for the underlying pipeline.

```r
res <- run_ascat(tumour_bams = bams, sex = "female", outdir = "out/")
```

### `run_sc_sequencing()` — `run_sc_sequencing.R`

The primary pipeline for single-cell and shallow-coverage whole-genome sequencing.

**Arguments:**

| Argument | Type | Description |
|---|---|---|
| `tumour_bams` | character vector | Paths to tumour BAM files, or a single 10X BAM |
| `sex` | character vector | `"male"` or `"female"` per sample |
| `allchr` | character vector | Chromosomes to analyse (e.g. `paste0("chr", 1:22)`) |
| `purs` | numeric vector | Purity values to search (e.g. `seq(0.1, 1, 0.01)`) |
| `ploidies` | numeric vector | Ploidy values to search (e.g. `seq(1.7, 5, 0.01)`) |
| `maxtumourpsi` | numeric | Maximum allowed tumour ploidy |
| `binsize` | integer | Bin size in bp (e.g. `30000`, `500000`) |
| `segmentation_alpha` | numeric | CBS significance threshold (smaller = fewer breakpoints) |
| `barcodes_10x` | character vector | Cell barcodes for 10X BAM mode |
| `normal_bams` | character vector | Normal BAM paths for reference construction |
| `outdir` | character | Output directory |
| `projectname` | character | Prefix for output files |
| `multipcf` | logical | Use joint multi-sample PCF segmentation |
| `smooth_sc` | logical | Smooth copy numbers across cells |
| `normal_barcodes` | character vector | Barcodes of diploid reference cells (10X mode) |
| `svinput` | data.frame | SV breakpoints to guide segmentation |

**What it does, step by step:**

1. **Loads reference bins.** If a precomputed `binsize` is available (30 kb for hg19/hg38), loads the bundled `lSe` and `lGCT` objects. Otherwise calls `getStartsEnds()` per chromosome. Stores in `res$lSe` and `res$lGCT`.

2. **Rebins to requested `binsize`.** Calls `treatlSe()`, `treatGCT()`, and `treatTrack()` to aggregate the fine-resolution reference bins into the user-specified window. Result stored in `res$nlSe`.

3. **Extracts coverage.** For each tumour BAM (or each barcode from a 10X BAM), calls `getTrackForAll()` / `getTrackForAll.10XBAM()`. Coverage is stored per-sample in `res$allTracks`.

4. **Builds a normal reference.** If `normal_bams` or `normal_barcodes` are provided, calls `combineDiploid()` to sum counts across normal samples. The normal reference is used later to subtract systematic biases.

5. **Excludes bad bins.** Optionally calls `sc_excludeBadBins()` to remove bins with extreme variance across cells.

6. **Segments the tracks.** Either jointly with `getLSegs.multipcf()` (when `multipcf=TRUE`) or individually via `segmentTrack()` / `segmentTrack_pcf()`. Segmented tracks per cell are stored in `res$allTracks.processed`.

7. **Grid searches purity and ploidy.** Calls `searchGrid()` per cell over the supplied `purs × ploidies` grid. Stores best solution per cell in `res$allSols`.

8. **Generates profiles.** Calls `fitProfile()` then `getProfile()` per cell. Stores in `res$allProfiles`.

9. **ML-guided refit (optional).** Calls `predictRefit_all()` which uses a pre-trained XGBoost model to select among the top local minima per cell. Updates `res$allProfiles`.

10. **Prints all results.** Calls `printResults_all()` to write TSVs, PNGs, and a summary table.

11. **Allele-specific CNA (optional).** If phase information is provided, calls `getAS_CNA()` and `getAS_CNA_smoothed()`.

**Returns:** An object of class `ascat.sc` (a named list). Printing it shows a clean overview rather than dumping every slot, and a set of accessor functions provides the commonly-needed pieces. All fields remain reachable directly with `$`, so existing code keeps working.

**Accessors (the "profiles" view):**

| Call | Returns |
|---|---|
| `summary(res)` / `getSummaryTable(res)` | Per-sample `samplename` / `purity` / `ploidy` / `ploidy.tumour` table |
| `getProfiles(res)` | List of total copy-number profiles |
| `getSolutions(res)` | List of best purity/ploidy solutions |
| `getProfilesAS(res)` / `getProfilesASsmoothed(res)` | Allele-specific profiles (single-cell only; `NULL` otherwise) |
| `getRefitted(res)` | Refitted profiles (`auto` / `manual`) if present |
| `getMetadata(res)` | Run parameters present on the object |

**Main slots (the "detailed" view):**

| Slot | Description |
|---|---|
| `res$lSe` | Reference bin coordinates (list by chromosome) |
| `res$nlSe` | Rebinned coordinates at analysis resolution |
| `res$lGCT` | GC content per bin (list by chromosome) |
| `res$allTracks` | Raw per-cell coverage tracks |
| `res$allTracks.processed` | Smoothed and segmented tracks per cell |
| `res$allSolutions` | Best purity/ploidy solution per cell |
| `res$allProfiles` | Copy number profile per cell (data.frame) |
| `res$allProfiles_AS` | Allele-specific profiles (if computed) |
| `res$summary` | Table with purity, ploidy per cell (written to disk by `printResults_all()`) |
| `res$mode` | Pipeline that produced the object: `"sc"`, `"targeted"`, or `"methylation"` |
| `res$binsize` | Bin size used |

The same `ascat.sc` object and accessors are returned by `run_targeted_sequencing()` and `run_methylation_array()`.

---

### `run_methylation_array()` — `run_methylation_arrays.R`

End-to-end pipeline for copy number calling from Illumina methylation arrays.

**Arguments:**

| Argument | Type | Description |
|---|---|---|
| `idat_dir` | character | Directory containing IDAT files |
| `rgSet` | RGChannelSet | Pre-loaded minfi object (alternative to `idat_dir`) |
| `sex` | character vector | `"male"` / `"female"` per sample |
| `purs` | numeric vector | Purity search grid |
| `ploidies` | numeric vector | Ploidy search grid |
| `platform` | character | `"450K"`, `"epicv1"`, or `"epicv2"` |
| `conumee` | logical | Use conumee binning instead of custom bins |
| `outdir` | character | Output directory |
| `projectname` | character | Output file prefix |

**What it does, step by step:**

1. Reads IDAT files with `minfi::read.metharray.exp()` or uses a supplied `rgSet`.
2. Extracts total intensity (sum of methylated + unmethylated) per probe.
3. Loads platform-specific bad probe lists.
4. Infers sex from chrX/chrY intensity ratios.
5. Builds a Panel of Normals (PoN) using sex-matched reference samples.
6. Normalises intensity to logR via `meth_getLogR()`.
7. Bins probes into genomic windows using `meth_bin()` or conumee bins.
8. Smooths and winsorises per-sample tracks (`meth_smooth()`, `meth_winsorise_ascat()`).
9. Segments each track.
10. Grid-searches purity and ploidy (`searchGrid()`).
11. Fits and stores profiles (`fitProfile()`, `getProfile()`).
12. Optionally refits with ML (`predictRefit_all()`).
13. Writes output via `printResults_all()`.

**Returns:** Same `res` structure as `run_sc_sequencing()`, with an additional `res$logr` matrix (probes × samples).

---

### `run_targeted_sequencing()` — `run_targeted_sequencing.R`

Copy number analysis using off-target reads from targeted/exome sequencing.

**Arguments:**

| Argument | Type | Description |
|---|---|---|
| `tumour_bams` | character vector | Tumour BAM paths |
| `bed_file` | character | BED file of capture targets (to exclude on-target reads) |
| `normal_bams` | character vector | Normal BAM paths |
| `purs` | numeric vector | Purity search grid |
| `ploidies` | numeric vector | Ploidy search grid |
| `multipcf` | logical | Joint segmentation |
| `outdir` / `projectname` | character | Output settings |

**What it does:** Identical pipeline to `run_sc_sequencing()` but with an extra initial step that maps the BED file to bins and removes all reads overlapping capture targets (`ts_getExcludeFromBedfile()`, `ts_removeOnTargets()`). This ensures only uninformative off-target reads drive the copy number signal, avoiding capture bias.

**Returns:** Same `res` structure.

---

## Coverage Extraction

### `getCoverageTrack()` — `getCoverageTrack.R`

Counts reads per genomic bin from a BAM file.

**Arguments:** `bamPath`, `chr`, `starts`, `ends`, `mapqFilter` (default 20)

**What it does:**
1. Builds a `ScanBamParam` specifying the genomic regions (one chromosome at a time) and a mapping quality filter.
2. Calls `Rsamtools::countBam()` to count reads landing in each bin.
3. Returns a data.frame with columns `records` (read count) and `nucleotides` (total bases) per bin.

**Returns:** `data.frame` with one row per bin: `chr`, `start`, `end`, `records`, `nucleotides`.

---

### `getCoverageTrack.Fix.R` — `getCoverageTrack.Fix.R`

Identical to `getCoverageTrack()` but handles edge cases in BAM files with non-standard chromosome naming. Used internally when the standard function fails.

---

### `getCoverageTrack.10XBAM()` — `getCoverageTrack.10XBAM.R`

Extracts per-barcode coverage from a 10X Genomics BAM file.

**Arguments:** `bamPath`, `chr`, `starts`, `ends`, `barcodes`, `pcchromosome`

**What it does:**
1. Scans the BAM for reads tagged with the `CB` (cell barcode) tag.
2. If `barcodes` is not provided, guesses valid barcodes by requiring each barcode to appear in at least `pcchromosome`% of bins on the chromosome.
3. For each barcode, counts reads and nucleotides per bin.
4. Returns a list indexed by barcode.

**Returns:** Named list of coverage data.frames, one entry per barcode.

---

### `getTrackForAll()` — `getTrackForAll.R`

Complete single-sample processing: coverage extraction + GC correction + segmentation.

**Arguments:** `bamfile`, `window`, `lCT` (pre-loaded counts), `lSe` (bin coords), `lGCT` (GC fractions), `lNormals` (normal reference), `segmentation_alpha`, `svinput`, `doSmooth`, `doSeg`

**What it does:**
1. Calls `getCoverageTrack()` per chromosome if `lCT` not provided.
2. Calls `smoothCoverageTrack()` for GC/normal correction.
3. Calls `segmentTrack()` (or PCF variant) per chromosome.
4. Assembles per-chromosome results into a unified list.

**Returns:** List with slots `lCTS.tumour` (raw coverage), `lSmooth` (GC-corrected logR), `lSegs` (segmentation).

---

### `getTrackForAll.10XBAM()` — `getTrackForAll.10XBAM.R`

Same as `getTrackForAll()` but for 10X barcoded BAMs. Runs `getCoverageTrack.10XBAM()` and then processes each barcode independently in parallel.

**Returns:** Named list of processed tracks per barcode.

---

### `getTrackForAll.bins.R` — `getTrackForAll.bins.R`

Variant of `getTrackForAll()` that operates on pre-defined bins without re-computing GC content. Used when bins and GC are already loaded from package data.

---

## Binning and Genomic Coordinate Functions

### `getStartsEnds()` — `getStartsEnds.R`

Generates genomic bin coordinates for a single chromosome.

**Arguments:** `window` (bin size in bp), `chr`, `lengthChr`, `dna` (reference sequence), `pathWindows`, `excludebed`, `centromeres`

**What it does:**
1. Divides the chromosome length into non-overlapping windows of size `window`.
2. If a BED file of regions to exclude is provided, removes bins overlapping those regions.
3. Calls `removeCentromeres()` to drop centromere bins.
4. Calls `excludeBadBins()` to drop bins with >5% N content in the reference FASTA.
5. Returns the final set of valid start/end coordinates.

**Returns:** List with `$starts` and `$ends` (integer vectors, length = number of valid bins).

---

### `getstartends()` — `getstartends.R`

Lightweight utility: given a total length and a window size, returns a set of start/end indices. Used internally by `treatlSe()` for rebinning.

**Returns:** List with `$starts` and `$ends`.

---

### `getBins_lSe_GC_from_fasta()` — `getBins_lSe_GC_from_fasta.R`

Batch wrapper: calls `getStartsEnds()` and `gcTrack()` for all chromosomes from a FASTA file in one call. Used to precompute the bundled package data objects.

**Returns:** List with `$lSe` (bin coordinates by chromosome) and `$lGCT` (GC fractions by chromosome).

---

### `treatlSe()` — `treatlSe.R`

Aggregates fine-resolution bin coordinates into larger windows.

**Arguments:** `lSe` (list of per-chromosome bin coords), `window` (number of fine bins per large bin)

**What it does:** For each chromosome, takes every `window`-th start and end coordinate, effectively merging `window` fine bins into one large bin. The merged bin spans from the start of the first fine bin to the end of the last.

**Returns:** New `lSe` list with coarser bin coordinates (same structure, fewer rows).

---

### `removeCentromeres()` — `removeCentromeres.R`

Removes bins that overlap centromeric regions.

**Arguments:** `lSe` (bin list for one chromosome), `centromeres` (data.frame with centromere coordinates), `chr`

**Returns:** Filtered `lSe` with centromeric bins removed.

---

### `excludeBadBins()` — `excludeBadBins.R`

Removes bins with high N-base content from a reference FASTA.

**Arguments:** `lSe`, `chr`, `dna` (DNAStringSet), `max.N.freq` (default 0.05)

**What it does:** For each bin, computes the fraction of N bases in the reference. Removes bins exceeding the threshold.

**Returns:** Filtered `lSe`.

---

### `getRefGenome()` — `getRefGenome.R`

Loads a reference genome FASTA into memory.

**Arguments:** `fasta` (file path), `CHRS` (chromosomes to load)

**Returns:** Named list of `DNAStringSet` objects, one per chromosome.

---

### `getBed_from_TS()` — `getBed_from_TS.R`

Reads and parses a BED file of targeted regions. Used in targeted sequencing mode to define the on-target regions that should be excluded.

**Returns:** `data.frame` with `chr`, `start`, `end` columns.

---

## GC Correction and Smoothing

### `gcTrack()` — `gcTrack.R`

Computes GC content per bin from a reference FASTA.

**Arguments:** `chr`, `starts`, `ends`, `dna`, `window` (sub-window within bin for GC estimation), `starts.exclude`, `ends.exclude`

**What it does:**
1. For each bin, extracts a central sub-window of the reference sequence.
2. Counts G and C bases.
3. Divides by sub-window length to get GC fraction.
4. If excluded regions are provided, subtracts their GC contribution.

**Returns:** Numeric vector of GC fractions, one per bin.

---

### `gcTrack.fixed()` — `gcTrack.fixed.R`

Alternative GC computation using a fixed 50 bp window centred on each bin midpoint, regardless of bin size. More stable for very large bins.

**Returns:** Numeric vector of GC fractions.

---

### `smoothCoverageTrack()` — `smoothCoverageTrack.R`

GC-corrects and optionally normal-subtracts a coverage track using LOESS regression.

**Arguments:** `lCT` (coverage data.frames by chromosome), `lSe` (bin coords), `lGCT` (GC fractions), `lNormals` (normal reference counts), `method`

**What it does:**
1. Log2-transforms read counts per bin: `log2(records / nucleotides)`.
2. If a normal reference is provided, log2-transforms the normal counts and subtracts them (tumour-vs-normal logR).
3. Fits a bivariate LOESS model: `logR ~ GC_fraction + read_length` using `myloess()`.
4. Residuals from this fit are the GC-corrected logR values.
5. Returns both the raw and corrected logR per bin.

**Returns:** Named list by chromosome; each element is a `data.frame` with columns `start`, `end`, `logr` (raw), `logr.fit` (LOESS-fitted), `logr.corrected` (residual).

---

### `myloess()` — `myloess.R`

Adaptive LOESS wrapper that automatically selects the `loess()` implementation based on sample size. For large datasets uses an approximate fast version to avoid memory issues.

**Arguments:** `LL` (number of data points), plus all standard `loess()` arguments passed via `...`

**Returns:** LOESS model fit object.

---

### `treatGCT()` — `treatGCT.R`

Aggregates per-fine-bin GC values to larger windows by taking the mean GC within each window.

**Arguments:** `lGCT` (list of GC vectors by chromosome), `window`

**Returns:** Rebinned list of GC vectors.

---

### `treatTrack()` — `treatTrack.R`

Aggregates raw coverage counts to larger bins by summing reads and nucleotides within each window.

**Arguments:** `lCTS` (list of coverage data.frames by chromosome), `window`

**Returns:** Rebinned list of coverage data.frames with summed `records` and `nucleotides`.

---

### `smoothNormals()` — `smoothNormals.R`

Fits a tumour logR against a panel of normal references using linear regression (iteratively excluding outliers). Used in methylation mode to remove batch effects.

**Arguments:** `logr` (sample logR vector), `lNormals` (matrix of normal logR vectors)

**Returns:** Normalised logR vector (residuals from the linear fit).

---

### `combineDiploid()` — `combineDiploid.R`

Merges multiple normal/diploid cell coverage tracks into a single reference by summing their read counts per bin.

**Arguments:** `lNormals` (list of per-cell coverage data.frames)

**Returns:** Single data.frame representing the pooled normal reference.

---

## Segmentation

### `segmentTrack()` — `segmentTrack.R`

Segments a single-chromosome logR track using Circular Binary Segmentation (CBS).

**Arguments:** `covtrack` (logR vector), `chr`, `starts` (bin positions), `sd` (noise estimate), `min.width` (min segment length), `ALPHA` (CBS significance), `SBDRY` (pre-computed segment boundaries)

**What it does:**
1. Creates a `CNA` object (from the `DNAcopy` package).
2. Applies `smooth.CNA()` to reduce outlier influence.
3. Calls `DNAcopy::segment()` with the specified `alpha` threshold and precomputed `SBDRY` boundaries for speed.
4. Returns the CBS segmentation result with segment means and boundaries.

**Returns:** `CNA` segment object with columns `loc.start`, `loc.end`, `num.mark`, `seg.mean`.

---

### `segmentTrack_pcf()` — `segmentTrack_pcf.R`

Segments a logR track using Piecewise Constant Fitting (PCF) from the `copynumber` package.

**Arguments:** `covtrack`, `chr`, `starts`, `sd`, `transform`, `min.width`, `ALPHA`

**Returns:** Formatted list with segment output compatible with the rest of the pipeline.

---

### `segmentCoverageTrack()` — `segmentCoverageTrack.R`

Wrapper that calls either `segmentTrack()` or `segmentTrack_pcf()` depending on the user-selected method. Handles per-chromosome iteration and error catching.

**Returns:** List of segmentation results by chromosome.

---

### `getLSegs.multipcf()` — `getLSegs.multipcf.R`

Joint multi-sample segmentation: segments all cells simultaneously using shared breakpoints.

**Arguments:** `allTracks` (list of raw coverage tracks per cell), `lCTS`, `lSe`, `lGCT`, `lNormals`, `segmentation_alpha`, `normalize`, `svinput`

**What it does:**
1. For each chromosome, stacks the logR tracks of all cells into a matrix.
2. Calls `copynumber::multipcf()`, which finds breakpoints shared across samples, improving sensitivity in single-cell data.
3. Applies original bin coordinates to the output.
4. If SV breakpoints are provided, calls `sc_getBreakpoints_multipcf()` to force breaks at SV positions.
5. Returns one processed track per cell, using the common breakpoints.

**Returns:** Named list of processed track objects, one per cell.

---

### `sc_getBreakpoints_multipcf()` — `sc_getBreakpoints_multipcf.R`

Guides PCF segmentation with known SV breakpoints.

**Arguments:** `tmpdata` (logR matrix), `svinput` (data.frame of SV breakpoints with `chr`, `pos`), `penalties`

**What it does:**
1. Finds bins closest to each SV breakpoint.
2. Splits the logR matrix at those positions.
3. Runs `multipcf()` independently on each fragment.
4. Reassembles the fragments into a single segmentation result.

**Returns:** Segmentation result with forced breaks at SV positions.

---

## Purity and Ploidy Estimation

### `searchGrid()` — `searchGrid.R`

Core optimisation: finds the best purity and ploidy by minimising copy number fit error over a grid.

**Arguments:** `tracksSingle` (processed track for one cell/sample), `purs`, `ploidies`, `forcepurity`, `maxTumourPhi`, `distance`, `ismale`, `isPON`, `gamma`

**What it does:**
1. Normalises the track to zero mean (`normaliseByPloidy()`).
2. Calls `buildDistanceMatrix()` to compute the fit error for every (purity, ploidy) combination.
3. Calls `findLocalMinima()` to identify the top candidate solutions.
4. For each local minimum, validates with `isGoodSolution()` (rejects solutions where >20 Mb have CN ≤ 0).
5. Returns the minimum-error valid solution.

**Returns:** List with:
- `$rho` — best purity
- `$psi` — best ploidy
- `$distance` — error value at best solution
- `$distancematrix` — full error grid (matrix)
- `$localMinima` — all detected local minima
- `$ambiguous` — logical flag (TRUE if multiple minima are close in error)

---

### `buildDistanceMatrix()` — `buildDistanceMatrix.R`

Constructs the purity × ploidy error matrix.

**Arguments:** `meansSeg` (segment means), `weights` (segment sizes), `purs`, `ploidies`, `maxTumourPhi`, `distance`, `gamma`, `sds`

**What it does:**
1. Iterates over all combinations of `purs` and `ploidies`.
2. Skips combinations where tumour ploidy (from `getTumourPhi()`) exceeds `maxTumourPhi`.
3. For each valid combination, calls `geterrors()` to compute the weighted fit error.
4. Fills the corresponding cell in the matrix.

**Returns:** Numeric matrix (rows = purity values, columns = ploidy values) of fit errors.

---

### `geterrors()` — `geterrors.R`

Computes the goodness-of-fit between observed segment logR values and expected integer copy numbers, for a given purity and ploidy.

**Arguments:** `rho` (purity), `phi` (ploidy), `meansSeg` (segment mean logR values), `weights` (segment lengths), `distance` (`"mse"` or `"statistical"`)

**What it does:**
1. Transforms each segment's mean logR to copy number using `transform_bulk2tumour()`.
2. Rounds to the nearest integer.
3. Computes the weighted sum of squared residuals between the observed logR and the expected logR for that integer CN.
4. For the `"statistical"` method, additionally incorporates segment-level uncertainty (standard deviation).

**Returns:** Single numeric value — the weighted error for this (purity, ploidy) combination.

---

### `getTumourPhi()` — `getTumourPhi.R`

Converts observed average ploidy to tumour-cell-specific ploidy, accounting for normal cell contamination.

**Formula:** `tumourPhi = (psi - 2 * (1 - rho)) / rho`

Where `psi` is the average observed ploidy (tumour + normal mixture) and `rho` is purity.

**Returns:** Numeric tumour ploidy.

---

### `transform_bulk2tumour()` — `transform_bulk2tumour.R`

Converts a bulk logR value to an absolute copy number in the tumour cells.

**Formula:** `CN = (psi_t * 2^(x/gamma) - 2 * (1 - rho)) / rho`

Where `psi_t` is tumour ploidy, `x` is the logR value, `gamma` corrects for non-linearity (important for methylation arrays), and `rho` is purity.

**Returns:** Numeric copy number value (not rounded).

---

### `normaliseByPloidy()` — `normaliseByPloidy.R`

Centres the segmented logR signal to a global mean of zero before grid search.

**Arguments:** `tracksSingle`, `ismedian` (use median instead of mean)

**What it does:** Computes the weighted mean (or median) of segment means across all chromosomes, then subtracts it from every segment. This accounts for global ploidy shift so that diploid segments are near logR = 0.

**Returns:** Updated track object with normalised segment means.

---

### `findLocalMinima()` — `findMinima.R`

Identifies local minima in the purity × ploidy error matrix.

**Arguments:** `mat` (error matrix), `N` (number of top minima to return)

**What it does:**
1. For each cell in the matrix, checks whether it is smaller than all 8 neighbours (Moore neighbourhood).
2. Clusters nearby minima that are within a proximity threshold.
3. Returns the top `N` by error value.

**Returns:** List with `$best` (coordinates of the global minimum) and `$all` (data.frame of all local minima sorted by error).

---

### `isGoodSolution()` — `isGoodSolution.R`

Validates that a purity/ploidy solution is biologically plausible.

**Arguments:** `meansSeg` (segment means with associated sizes), `maxSize0` (max allowed size of CN ≤ 0 segments in bp)

**Logic:** Returns `TRUE` if the total genomic length of segments with copy number ≤ 0 is less than 20 Mb. Solutions with large deletions to zero (or below) are rejected as likely fitting artefacts.

**Returns:** `logical` (TRUE = valid solution).

---

## Copy Number Profile Generation

### `fitProfile()` — `fitProfile.R`

Fits a copy number profile for a given purity and ploidy.

**Arguments:** `tracksSingle`, `purity`, `ploidy`, `ismale`, `isPON`, `gamma`, `ismedian`, `gridsearch`

**What it does:**
1. Normalises the track to zero mean.
2. For each segmented region on each chromosome, extracts the mean and standard deviation of the raw logR values within the segment boundaries.
3. Applies `transform_bulk2tumour()` to convert segment mean logR to copy number.
4. Stores per-segment statistics: mean, SD, number of bins, rounded copy number.

**Returns:** Per-chromosome list; each element is a list of segments, where each segment contains `$start`, `$end`, `$mu` (mean logR), `$sd`, `$num.mark`, `$roundmu` (rounded CN).

---

### `getProfile()` — `getProfile.R`

Formats the output of `fitProfile()` into a standardised flat data.frame.

**Arguments:** `profile` (output of `fitProfile()`), `CHRS` (chromosome names)

**What it does:** Iterates over chromosomes and segments, binding all segment information into one data.frame.

**Returns:** `data.frame` with columns:

| Column | Description |
|---|---|
| `chromosome` | Chromosome name |
| `start` | Segment start position (bp) |
| `end` | Segment end position (bp) |
| `num.mark` | Number of bins in segment |
| `total_copy_number` | Rounded integer copy number |
| `total_copy_number_logr` | Non-integer fitted CN |
| `logr` | Observed mean logR of segment |
| `logr.sd` | Standard deviation of logR within segment |

---

### `sc_getCalls()` — `sc_getCalls.R`

Generates discrete copy number calls for each segment.

**Arguments:** `tracksSingle`, `purity`, `ploidy`

**Returns:** List of segments per chromosome with rounded integer CN values.

---

### `sc_getMatrixOfCalls()` — `sc_getMatrixOfCalls.R`

Compiles per-cell copy numbers into a genome-wide matrix.

**Arguments:** `listTracks` (list of processed tracks per cell), `lSolutions` (list of solutions per cell), `lSe` (bin coordinates)

**What it does:** For each cell, generates CN calls per bin by propagating segment-level calls to their constituent bins. Stacks all cells into a matrix.

**Returns:** Numeric matrix with dimensions `(number of bins) × (number of cells)`. Values are integer copy numbers.

---

## Allele-Specific Copy Numbers

### `getAS_CNA()` — `getAS_CNA.R`

Derives allele-specific (phased) copy numbers from SNP data.

**Arguments:** `res`, `path_to_phases` (file with phasing information), `list_ac_counts_paths` (paths to allelic count files), `purs`, `ploidies`, `outdir`, `projectname`

**What it does:**
1. Reads phasing data: for each heterozygous SNP, records which haplotype (0|1 or 1|0) carries the reference allele.
2. Reads allelic count files (VCF-like format with reference and alternative allele counts per SNP per cell).
3. Merges allelic counts with phasing information to assign counts to haplotype A or B.
4. For each copy number segment, fits a binomial model to estimate the BAF (B-allele frequency).
5. From total CN and BAF, estimates the number of copies of each allele (nA = major allele, nB = minor allele).

**Returns:** Updated `res` with `res$allProfiles_AS`, where each profile data.frame has additional columns `nA` (major allele CN) and `nB` (minor allele CN).

---

### `getAS_CNA_smoothed()` — `getAS_CNA_smoothed.R`

Smooths allele-specific copy numbers across multiple cells.

**Arguments:** `res`, `mc.cores`

**What it does:**
1. Across all cells, collects BAF and total CN estimates per segment.
2. Fits a mixture model (EM) where each component corresponds to a valid integer allele combination (e.g. 1+1, 2+0, 2+1, etc.).
3. Assigns each segment to the most probable allele combination.
4. Uses information shared across cells to resolve ambiguous single-cell assignments.

**Returns:** Updated `res` with consensus allele-specific profiles.

---

## Quality Control and Filtering

### `getFilters()` — `getFilters.R`

Identifies low-quality cells to exclude from downstream analysis.

**Arguments:** `res`, `probs` (quantile thresholds), `threshold_extra_noise`, `outdir`

**What it does:**
1. **Read depth filter:** Computes total reads per cell. Cells below the low-depth threshold are flagged.
2. **Noise filter:** For each cell, computes the logR noise as the standard deviation of consecutive bin differences. Cells with noise above a LOESS-fitted upper envelope are flagged.
3. **Ambiguity filter:** Cells where `searchGrid()` set the `ambiguous` flag (multiple competing solutions) are flagged.
4. Combines all filters into a binary pass/fail per cell.

**Returns:** Updated `res` with `res$filters` (named logical vector) and `res$filters_df` (data.frame with per-filter details).

---

### `sc_excludeBadBins()` — `sc_excludeBadBins.R`

Removes systematically unreliable bins from the analysis.

**Arguments:** `res`

**What it does:**
1. Across all cells, computes the variance of each bin's logR.
2. Ranks bins by variance.
3. Removes bins in the top and bottom decile (extreme variance likely indicates mappability artefacts or segmental duplications).

**Returns:** Updated `res$lSe`, `res$nlSe`, `res$lGCT` with bad bins removed.

---

### `identify_diploid_cells()` — `identify_diploid_cells.R`

Identifies cells that are likely normal diploid, for use as a reference.

**Arguments:** `res`, `ploidy_range` (acceptable ploidy window around 2.0), `exclude_ambiguous`, `min_diploid_fraction`

**What it does:**
1. Selects cells with ploidy within `ploidy_range` of 2.0.
2. Optionally excludes cells flagged as `ambiguous` by `searchGrid()`.
3. For each remaining cell, computes the fraction of the genome at CN = 2.
4. Keeps cells where this fraction exceeds `min_diploid_fraction`.

**Returns:** Character vector of cell names passing all filters.

---

### `filterBins()` — `filterBins.R`

Statistical bin filtering based on cross-sample consistency.

**Arguments:** `allTracks`, `logr`, `lSe`, `fraction`, `k`

**What it does:**
1. Computes the median logR per bin across all cells.
2. Computes the variability (IQR or SD) of each bin across cells.
3. Flags bins in the extreme `fraction` tails of the median distribution.
4. Returns indices of bins to remove.

**Returns:** List of bin indices to remove per chromosome.

---

### `get_filters_targeted_sequencing()` — `get_filters_targeted_sequencing.R`

Identifies unreliable bins specific to targeted sequencing off-target data.

**Arguments:** `res`

**What it does:** Similar to `sc_excludeBadBins()` but tuned to the expected off-target coverage patterns, where read depth is much lower and more variable.

**Returns:** Updated `res` with unreliable bins removed.

---

## Refit and Refinement

### `predictRefit_all()` — `predictRefit_all.R`

Machine-learning-guided solution selection using a pre-trained XGBoost model.

**Arguments:** `res`, `ismedian`, `gamma`

**What it does:**
1. Loads a pre-trained XGBoost model bundled with the package (trained on manually reviewed copy number profiles).
2. For each cell, extracts the top N local minima from the purity/ploidy error landscape.
3. Calls `get_fit_features()` to compute features for each candidate solution.
4. Applies the XGBoost model to score each candidate.
5. Selects the candidate with the highest predicted quality score.
6. Re-fits the profile at the selected solution.

**Returns:** Updated `res$allProfiles` and `res$allSols` with ML-selected solutions.

---

### `predictRefit()` — `predictRefit.R`

Legacy single-sample version of `predictRefit_all()`. Kept for backward compatibility.

---

### `get_fit_features()` — `get_fit_features.R`

Extracts numerical features from a copy number profile for ML scoring.

**Arguments:** `fit_cn` (profile data.frame from `getProfile()`)

**Features computed:**
- Proportion of genome at each integer CN state (0, 1, 2, 3, 4, 5, 6, >6, <0)
- Median, 5th, and 95th percentiles of the logR distribution
- Deviation between observed logR and expected logR for the rounded CN
- Number of segments

**Returns:** Single-row `data.frame` of feature values.

---

### `refitProfile()` — `refitProfile.R`

Manual refit using user-specified known copy numbers in two reference segments.

**Arguments:** `track`, `solution`, `chr1`, `ind1`, `n1` (known CN at region 1), `chr2`, `ind2`, `n2` (known CN at region 2), `gridpur`, `gridpl`

**What it does:**
1. From the observed logR at the two reference segments and their known integer CN values, algebraically solves for the implied purity and ploidy.
2. Runs `searchGrid()` in a narrow window around the computed point.

**Returns:** Refined solution list.

---

### `refitProfile_shift()` — `refitProfile_shift.R`

Tries shifting the entire copy number scale by ±1 and checks if the new solution is better.

**Arguments:** `track`, `solution`, `shift` (+1 or -1), `gridpur`, `gridpl`

**What it does:** Adjusts the global CN offset by `shift`, re-fits, and returns the new solution if it has a lower error than the original.

**Returns:** Updated or original solution.

---

### `run_any_refitProfile()` — `run_any_refitProfile.R`

User-facing wrapper to manually refit a specific sample within a `res` object.

**Returns:** Updated `res` with refitted solution for the specified sample.

---

### `run_any_refitProfile_shift()` — `run_any_refitProfile_shift.R`

Wrapper around `refitProfile_shift()` for one named sample.

**Returns:** Updated `res`.

---

### `runPloidyModifier()` — `runPloidyModifier.R`

Launches an interactive Shiny application for manual inspection and correction of purity/ploidy solutions. Allows the user to browse the sunrise plot, click a new solution, and refit.

**Returns:** Does not return; launches a browser-based GUI.

---

## Methylation Array Functions

### `meth_getLogR()` — `meth_getLogR.R`

Derives probe-level logR from total methylation intensity using Panel of Normals normalisation.

**Arguments:** `totalintensity` (matrix: probes × samples), `totalintensityNormalM` (male normal matrix), `totalintensityNormalF` (female normal matrix), `sex`, `annot` (probe annotation), `GAMMA`

**What it does:**
1. For each tumour sample, selects sex-matched normal references.
2. Fits a linear model: `log2(tumour_intensity) ~ log2(normal_intensities)`.
3. Iteratively re-fits excluding probes with high residuals (outlier probes).
4. Residuals are the normalised logR values.
5. Applies platform-specific `GAMMA` correction for chrX and chrY.

**Returns:** Matrix of logR values (probes × samples).

---

### `meth_derive_pon()` — `meth_derive_pon.R`

Builds a Panel of Normals from a directory of normal IDAT files.

**Arguments:** `idat_dir`, `allchr`, `GAMMA`

**What it does:**
1. Reads all IDAT files in the directory.
2. Preprocesses and computes total intensity per probe.
3. Infers sex of each sample.
4. Computes logR for all samples against each other.
5. Identifies probes with consistently high variability across normals (bad probes).

**Returns:** List with `$totalintensity` (normal intensity matrix) and `$bad_loci` (logical vector of unreliable probes).

---

### `meth_smooth()` — `meth_smooth.R`

Bins and smooths probe-level methylation logR into genomic windows.

**Arguments:** `logr` (probe-level vector), `starts`, `ends` (probe positions), `step` (bin size), `is_isolength` (equal bin length vs. equal probe count)

**What it does:**
1. Divides probes into bins by position.
2. For each bin, computes the median logR across probes.
3. Merges bins with fewer than a minimum number of probes into adjacent bins.
4. Filters edge bins at chromosome boundaries.

**Returns:** `data.frame` with `start`, `end`, `logr` (median) columns.

---

### `meth_smoothAllChr()` — `meth_smoothAllChr.R`

Applies `meth_smooth()` across all chromosomes for a single sample.

**Arguments:** `logr`, `annot` (probe annotations with chromosome/position), `step`

**Returns:** Concatenated smoothed track across all chromosomes.

---

### `meth_bin()` — `meth_bin.R`

Assigns methylation probes to predefined genomic bins (e.g., conumee bins) and computes the median logR per bin.

**Arguments:** `logr`, `starts`, `ends`, `chrs`, `grbins` (GRanges of target bins)

**Returns:** List with `$logr` (binned values), `$starts`, `$ends`, `$chrs`.

---

### `meth_winsorise_ascat()` — `meth_winsorise_ascat.R`

Robust winsorisation of methylation logR to reduce the influence of extreme outlier probes.

**Arguments:** `x` (logR vector), `apply_wins` (logical)

**What it does:**
1. Applies a running median filter.
2. Computes MAD (median absolute deviation) per local window.
3. Clips values beyond median ± k×MAD using asymmetric thresholds (allowing for natural copy number gains).

**Returns:** Winsorised logR vector.

---

### `meth_transform_bulk2tumour()` — `meth_transform_bulk2tumour.R`

Methylation-specific wrapper for `transform_bulk2tumour()` that applies the array non-linearity correction factor `GAMMA`.

**Returns:** Numeric copy number values.

---

### `meth_winsorise_ascat()` — (see above)

---

## Targeted Sequencing Functions

### `ts_treatBed()` — `ts_treatBed.R`

Parses a BED file and expands genomic regions to individual base-pair resolution for overlap testing.

**Returns:** Expanded region coordinates.

---

### `ts_getExcludeFromBedfile()` — `ts_getExcludeFromBedfile.R`

Maps BED capture regions to analysis bins, identifying which bins should be excluded (on-target).

**Returns:** List of bin indices to exclude per chromosome.

---

### `ts_getTrackForAll.excludeTargets()` — `ts_getTrackForAll.excludeTargets.R`

Extracts off-target coverage by reading BAM but ignoring reads that overlap capture targets.

**Returns:** Off-target coverage tracks.

---

### `ts_mergeCountsNormal()` — `ts_mergeCountsNormal.R`

Pools coverage counts from multiple normal BAMs into a single reference.

**Returns:** Merged normal coverage data.frame.

---

### `ts_removeOnTargets()` — `ts_removeOnTargets.R`

Filters out bins overlapping capture targets from an existing coverage track.

**Returns:** Filtered coverage data.frame.

---

### `ts_smoothCoverageTrackAll()` — `ts_smoothCoverageTrackAll.R`

Batch application of GC correction and smoothing to all samples in targeted sequencing mode.

**Returns:** List of smoothed tracks per sample.

---

### `ts_smoothTrack()` — `ts_smoothTrack.R`

Single-sample LOESS smoothing for targeted sequencing off-target tracks.

**Returns:** Smoothed logR values.

---

## Visualization and Output

### `printResults_all()` — `printResults_all.R`

Master output function: generates all plots and files for a complete analysis.

**Arguments:** `res`, `outdir`, `projectname`, `is_pdf`, `rainbowChr`

**What it does:**
1. Creates output subdirectories.
2. For each cell/sample, calls `plotSolution()` and saves as PNG or PDF.
3. Writes per-cell CN profiles to TSV via `writeProfile()`.
4. If `predictRefit_all()` was run, also writes refitted profiles.
5. Compiles a summary table (`results_summary.tsv`) with columns: sample name, purity, ploidy, tumour ploidy, ambiguous flag, filter flags.
6. Saves the full `res` object as `<projectname>.RData`.

**Returns:** Updated `res` with `res$results_summary` added.

---

### `plotSolution()` — `plotSolution.R`

Generates a genome-wide copy number plot for one sample.

**Arguments:** `tracksSingle`, `purity`, `ploidy`, `ylim`, `gamma`, `ismale`, `allchr`, `rainbowChr`, `svinput`

**What it produces:**
- X-axis: genomic position (cumulative across chromosomes)
- Y-axis: LogR value
- Grey dots: individual bin logR values
- Coloured horizontal lines: integer CN levels (red = amplification, blue = deletion, black = diploid)
- Thick bars: fitted segment CN values
- Vertical grey lines: chromosome boundaries with labels
- Top annotation: purity, ploidy, tumour ploidy, error
- Optional: SV breakpoints marked as vertical lines

**Returns:** None (produces a plot).

---

### `plotSunrise()` — `plotSunrise.R`

Heatmap of the purity × ploidy error landscape ("sunrise plot").

**Arguments:** `solution` (output of `searchGrid()`), `localMinima`, `N`

**What it produces:**
- Heatmap where colour encodes the fit error at each (purity, ploidy) combination
- Low-error regions appear as "valleys"
- Best solution marked with a cross
- Local minima marked with circles

**Returns:** Local minima information if `N > 0`.

---

### `plot_AS_profile()` — `plot_AS_profile.R`

Plots allele-specific copy numbers across the genome.

**Arguments:** `t` (data.frame with `nA`, `nB`, `chr`, `start`, `end`, `fitted`)

**What it produces:**
- Dark blue line: minor allele (nB) copy number per segment
- Orange line: major allele (nA) copy number per segment
- Segments are plotted at genomic position with chromosome boundaries

**Returns:** None (produces a plot).

---

### `sc_plotHeat()` — `sc_plotHeat.R`

Hierarchical clustering heatmap of the copy number matrix across all cells.

**Arguments:** `mat` (CN matrix: bins × cells), `fundist` (distance function), `funclust` (clustering function), `centromeres`

**What it produces:**
- Rows: cells (clustered by CN similarity)
- Columns: genomic bins
- Colour: red = deletion, white = diploid, blue = amplification
- Chromosome boundaries marked

**Returns:** None (produces a plot).

---

### `sc_plotGenome()` — `sc_plotGenome.R`

Genome-wide stacked area plot showing frequency of gains and losses across all cells.

**Arguments:** `mat` (CN matrix), `scaleY`, `centromeres`

**What it produces:**
- Stacked bars per genomic position
- Blue: proportion of cells with loss (CN < 2)
- Red: proportion of cells with gain (CN > 2)
- Useful for identifying recurrent alterations

**Returns:** Cumulative position breaks (for axis labelling).

---

### `sc_getCols()` — `sc_getCols.R`

Maps copy number values to RGB colour values for heatmap rendering.

**Arguments:** `mmm` (matrix of CN deviations from diploid, range −2 to +2), `channel` (`"R"`, `"G"`, or `"B"`)

**Returns:** Numeric vector of colour intensities (0–255).

---

### `writeProfile()` — `writeProfile.R`

Writes a single sample's copy number profile to a tab-separated file.

**Arguments:** `prof` (data.frame from `getProfile()`), `samplename`, `outdir`

**Output file:** `<outdir>/<samplename>.profile.tsv`

**Returns:** None (writes file).

---

## Utility and Helper Functions

### `checkArguments_scs()` — `checkArguments_scs.R`

Validates all arguments to `run_sc_sequencing()` before the analysis starts.

**Checks performed:**
- All BAM files exist on disk
- `sex` vector matches length of `tumour_bams`
- `outdir` exists and is writable
- `purs` and `ploidies` are valid numeric ranges
- Result object (if supplied) is complete

**Returns:** None (calls `stop()` on first invalid argument).

---

### `checkArguments_meth()` — `checkArguments_meth.R`

Same validation logic for `run_methylation_array()`.

---

### `getlInds()` — `getlInds.R`

For each chromosome, finds which bin indices overlap a set of exclusion regions.

**Returns:** List of integer vectors (bin indices to remove) per chromosome.

---

### `getlGCT_excluded()` — `getlGCT_excluded.R`

Removes GC fraction values for excluded bins.

**Arguments:** `lGCT`, `lInds` (from `getlInds()`)

**Returns:** Filtered GC list.

---

### `getlSe_excluded()` — `getlSe_excluded.R`

Removes bin coordinate entries for excluded bins.

**Arguments:** `lSe`, `lInds`

**Returns:** Filtered bin coordinate list.

---

### `getnlCTS_excluded()` — `getnlCTS_excluded.R`

Removes coverage rows for excluded bins from all samples.

**Arguments:** `nlCTS` (coverage data.frames), `lInds`

**Returns:** Filtered coverage list.

---

### `findBestSolution()` — `findBestSolution.R`

Deprecated predecessor to `searchGrid()`. Kept for backward compatibility. Uses the same grid search logic but without local minima tracking or ML refit capability.

---

## Key Data Structures

### `res$lSe` — Reference Bin Coordinates

A named list where each element corresponds to one chromosome. Each element is a list with:
- `$starts` — integer vector of bin start positions (bp)
- `$ends` — integer vector of bin end positions (bp)

This is the genome-wide definition of the analysis bins. It is loaded from package data at the start of the pipeline and never changes per-sample. All subsequent data structures are indexed consistently with these bins.

---

### `res$nlSe` — Rebinned Coordinates

Same structure as `res$lSe` but at the user-requested `binsize` resolution (after `treatlSe()`). This is what is actually used for analysis when the user's binsize differs from the 30 kb precomputed default.

---

### `res$lGCT` — GC Content

Named list by chromosome; each element is a numeric vector of GC fractions, one per bin in `res$lSe`.

---

### `res$allTracks` — Raw Coverage Tracks

Named list by sample/cell. Each element is itself a list by chromosome, where each chromosome entry is a `data.frame` with:
- `start`, `end` — bin positions
- `records` — read count
- `nucleotides` — total bases counted

---

### `res$allTracks.processed` — Smoothed and Segmented Tracks

Named list by sample/cell. Each element is a list by chromosome, containing:
- `logr` — GC-corrected logR per bin
- `lSegs` — segmentation result (from CBS or PCF)
- `lSegs$seg.mean` — mean logR per segment
- `lSegs$loc.start`, `lSegs$loc.end` — segment boundaries in bp

This is the key input to `searchGrid()` and `fitProfile()`.

---

### `res$allSols` — Solutions

Named list by sample/cell. Each element contains:
- `$rho` — purity (0–1)
- `$psi` — ploidy (average, mixture of tumour and normal)
- `$distance` — fit error
- `$distancematrix` — full purity × ploidy error grid
- `$localMinima` — data.frame of top local minima
- `$ambiguous` — logical

---

### `res$allProfiles` — Copy Number Profiles

Named list by sample/cell. Each element is a `data.frame` (output of `getProfile()`) with columns:
`chromosome`, `start`, `end`, `num.mark`, `total_copy_number`, `total_copy_number_logr`, `logr`, `logr.sd`

This is the primary output of the pipeline.

---

### `res$results_summary` — Summary Table

`data.frame` with one row per sample/cell and columns:
- `sample` — sample/cell name
- `purity` — estimated purity
- `ploidy` — estimated average ploidy
- `tumour_ploidy` — estimated tumour cell ploidy
- `ambiguous` — TRUE if multiple competing solutions exist
- `filter` — TRUE if cell passed QC filters
- `n_segs` — number of copy number segments

---

## Dependencies

### Bioconductor

| Package | Purpose |
|---|---|
| `Rsamtools` | BAM file reading and read counting |
| `Biostrings` | Reference FASTA loading and GC computation |
| `GenomicRanges` | Genomic interval arithmetic |
| `DNAcopy` | CBS segmentation algorithm |
| `minfi` | Methylation array IDAT reading and preprocessing |
| `conumee` | Alternative methylation binning strategy |

### CRAN

| Package | Purpose |
|---|---|
| `copynumber` | PCF and multipcf segmentation |
| `xgboost` | ML model for solution scoring |
| `parallel` | Multi-core parallel processing (`mclapply`) |
| `shiny` + plugins | Interactive ploidy modifier GUI |
| `spatstat.geom` | Spatial geometry utilities |
| `dipsaus` | Miscellaneous utilities |

### Bundled Package Data

| Object | Description |
|---|---|
| `lSe.hg19.filtered` / `lSe.hg38.filtered` | Precomputed 30 kb bin coordinates for hg19/hg38 |
| `lGCT.hg19.filtered` / `lGCT.hg38.filtered` | Precomputed GC fractions for those bins |
| `lSe_unfiltered_5000.mm39` | 5 kb bins for mouse genome mm39 |
| XGBoost model | Pre-trained model for `predictRefit_all()` |
| SBDRY matrices | Precomputed CBS boundary matrices for fast segmentation |

---

*This document covers all functions in ASCAT.sc v0.1. For usage examples see the [project Wiki](https://github.com/VanLoo-lab/ASCAT.sc/wiki).*
