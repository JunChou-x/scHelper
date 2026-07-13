# scDistiller

`scDistiller` is a lightweight toolkit designed for single-cell RNA-seq analysis. It aims to provide high-quality matrix inputs for downstream advanced analyses (such as cell-cell communication analysis with CellChat/NicheNet) through automated clustering optimization and precise data cleaning.

## Installation

You can install `scDistiller` from GitHub using `devtools`:

```R
devtools::install_github("JunChou-x/scDistiller")
Core Features
OptimizeRes: Automated clustering optimization
Automatically identifies the optimal resolution based on the Silhouette Score and supports downsampling strategies for large-scale datasets.
PlotTopMarkers: Rapid feature extraction
Quickly identifies marker genes using the presto package and generates publication-quality standardized heatmaps.
PurifyCells: General anchor-based filtering engine
Precisely extracts specific cellular niches through core gene markers and whitelist-based filtering strategies while removing background noise.
CheckStability: Clustering stability assessment
Evaluates the reliability of clustering results using bootstrap resampling.
Usage Example
library(scDistiller)

# 1. Automatically identify the optimal clustering resolution
pbmc <- OptimizeRes(pbmc)

# 2. Rapid marker visualization
markers <- PlotTopMarkers(pbmc)

# 3. Purify target cell populations
# Example: retaining cardiomyocytes and endothelial cells
pbmc <- PurifyCells(
    pbmc,
    core_markers = c("Tnnt2", "Myl2"),
    whitelist_markers = c("Pecam1", "Cdh5")
)

# 4. Evaluate clustering stability
stability <- CheckStability(pbmc)
Contribution

This toolkit is developed and maintained by JunChou-x.

If you encounter any bugs or have suggestions, please submit an Issue.
