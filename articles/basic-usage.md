# Basic Usage on PBMC3k Data

``` r
suppressPackageStartupMessages({
library(Seurat)
library(SeuratData)
library(recall)
})
```

First, we use the `SeuratData` data package to first download and then
load 2700 PBMCs. The loaded `SeuratObject`, `pbmc3k`, is from an old
version of `Seurat`, and so we update the object to v5.

``` r
set.seed(123)

SeuratData::InstallData("pbmc3k")
data("pbmc3k")

pbmc3k <- UpdateSeuratObject(pbmc3k)
```

Now, we use `Seurat` to perform the usual preprocessing steps that are
performed prior to clustering.

``` r
pbmc3k <- NormalizeData(pbmc3k)
pbmc3k <- FindVariableFeatures(pbmc3k)
pbmc3k <- ScaleData(pbmc3k)
pbmc3k <- RunPCA(pbmc3k)
pbmc3k <- FindNeighbors(pbmc3k)
pbmc3k <- RunUMAP(pbmc3k, dims = 1:10)
```

The `recall` algorithm can be run with a single function call as a
drop-in replacement for the `Seurat` function `FindClusters`.

``` r
pbmc3k <- FindClustersRecall(pbmc3k)
```

The `recall` clusters are set to the idents of the `SeuratObject` that
is returned by `FindClustersRecall`

``` r
DimPlot(pbmc3k)
```

Cluster labels from `FindClustersRecall` care stored in the metadata in
the column `pbmc3k@meta.data$recall_clusters`.

``` r
DimPlot(pbmc3k, group.by = "recall_clusters")
```
