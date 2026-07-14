# scHelper: An Automated Full-Lifecycle Single-Cell RNA-seq Analysis Assistant

## Overview

We introduce **scHelper** to streamline and automate the entire single-cell transcriptomics analysis workflow. Built seamlessly on top of Seurat, `scHelper` minimizes subjective manual tuning by providing data-driven parameter optimization (e.g., MAD-based QC thresholds, geometric chord distance for HVG/PC selection) and robust cellular context refinement for downstream host-microbe or inter-cellular communication analyses.

## How to run

Since `scHelper` is an R package, the entire automated pipeline is executed directly within your R script or RStudio environment. The standard workflow follows a sequential function execution:

```r
shQC() -> shDoublet() -> shHVG() -> shPC() -> shCellCycle() -> shBatch() -> shCluster() -> shMarkers() -> shFilter() -> shStability()

library(Seurat)
library(scHelper)

# Load your raw Seurat object
object <- readRDS("data/raw_seurat_object.rds")

# 1. Preprocessing & Quality Control
object <- shQC(object, mt_pattern = "^MT-", n_mads = 3)

object <- subset(
    object,
    qc_pass == "Pass"
)

object <- shDoublet(
    object,
    run_scDblFinder = TRUE,
    batch_col = "sample_id",
    ncores = 4
)

# 2. Normalization & Dimensionality Reduction
object <- NormalizeData(object)

object <- shHVG(
    object,
    max_features = 5000
)

object <- ScaleData(object)

object <- RunPCA(object)

object <- shPC(
    object,
    max_pc = 50
)

shBatch(
    object,
    batch_col = "sample_id"
)

# 3. Clustering & Refinement
optimal_pc <- object@misc[["scHelper_optimal_pc"]]

object <- FindNeighbors(
    object,
    dims = 1:optimal_pc
)

object <- RunUMAP(
    object,
    dims = 1:optimal_pc
)

object <- shCluster(
    object,
    res_range = seq(0.2, 1.2, 0.2)
)

# 4. Marker Extraction & Target Filtering
top_markers <- shMarkers(
    object,
    top_n = 5
)

object <- shFilter(
    object,
    core_markers = c("ALB", "HNF4A")
)

stability_res <- shStability(
    object,
    n_iter = 10
)
```

## Parameters

Instead of command-line configuration files, scHelper utilizes a suite of prefix-standardized (`sh*`) R functions. Here are the core parameters across key modules:

- **object**: The primary input containing single-cell RNA sequencing data in Seurat format. It stores gene expression information, along with any associated annotations and metadata.

- **mt_pattern**: A regular expression used in `shQC` to identify mitochondrial genes (e.g., `^MT-` for human data, `^mt-` for mouse data).

- **n_mads**: The number of Median Absolute Deviations (MADs) used in `shQC` to automatically define adaptive filtering thresholds.

- **batch_col**: The metadata column specifying sample or batch groups, required for calculating expected doublet rates in `shDoublet` and evaluating batch effects in `shBatch`.

- **max_features**: The maximum number of highly variable genes (HVGs) to evaluate in `shHVG` using the geometric chord distance method.

- **res_range**: A numeric vector specifying the range of resolutions to test for optimal clustering in `shCluster` (e.g., `seq(0.2, 1.2, 0.2)`).

- **core_markers**: A character vector of essential genes used in `shFilter` to precisely isolate and label the target cell populations.

## Outputs

Running the scHelper workflow appends highly structured metadata and records directly into your Seurat object, alongside visual diagnostics:

- **Metadata Columns (`object@meta.data`)**: Adds informative columns such as `qc_pass` (Pass/Fail status), `scDblFinder_class` (Singlet/Doublet), `optimal_clusters` (auto-selected IDs based on Silhouette score), and `distilled_annotation` (Kept/Excluded for downstream pure subsetting).

- **Misc Records (`object@misc`)**: Stores pipeline-derived optimal parameters (e.g., `scHelper_optimal_pc`, `scHelper_best_res`, `scHelper_doublet_rate`) to ensure automated linkage between upstream and downstream functions.

- **Visual Diagnostics**: High-quality, publication-ready `ggplot2` diagnostic plots (e.g., QC violin plots with cutoffs, LISI density curves, Elbow plots, Silhouette tracking lines, Marker Heatmaps) are printed directly to the R plotting window during execution.

## Installation and Requirements

Before running the workflow, you need to install the required software and dependencies:

- R 4.0+

- Required Core R packages:
  - Seurat (V4 or V5)
  - ggplot2
  - dplyr
  - tidyr
  - cluster

## 1. Install scHelper R Package

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

devtools::install_github("JunChou-x/scHelper")
```

## 2. Install Advanced Module Dependencies

To unlock the full potential and extreme speed of scHelper advanced features, please install the following specialized backend engines:

```r
# [Required for shDoublet] Robust doublet detection and multi-core acceleration

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(
    c(
        "scDblFinder",
        "BiocParallel"
    )
)

# [Required for shMarkers] Ultra-fast marker identification via presto

devtools::install_github("immunogenomics/presto")

# [Required for shBatch] Quantitative batch effect evaluation via LISI

devtools::install_github("immunogenomics/lisi")
```
