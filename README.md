## Core Features

### OptimizeRes: Automated clustering optimization

Automatically determines the optimal clustering resolution based on the **Silhouette Score**. It supports downsampling strategies for large-scale single-cell datasets to improve computational efficiency.

### PlotTopMarkers: Fast marker identification and visualization

Rapidly identifies highly enriched marker genes using the **presto** package and generates publication-quality standardized heatmaps for marker visualization.

### PurifyCells: General anchor-based cell purification engine

A flexible filtering framework that uses core marker genes and whitelist-based strategies to accurately extract specific cellular niches while reducing background contamination and unwanted signals.

### CheckStability: Clustering stability assessment

Evaluates clustering robustness through bootstrap resampling, providing quantitative assessment of clustering reliability.

## Usage Example

```r
library(scDistiller)

# 1. Automatically optimize clustering resolution
pbmc <- OptimizeRes(pbmc)

# 2. Rapid marker identification and visualization
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
