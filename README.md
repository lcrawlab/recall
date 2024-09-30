# recall (Calibrated Clustering with Artificial Variables) <img src="man/figures/recall_logo.png" align="right" alt="" width="120"/>

[![R CMD check](https://github.com/lcrawlab/recall/actions/workflows/check-standard.yml/badge.svg)](https://github.com/lcrawlab/recall/actions/workflows/check-standard.yml)
[![Docker Image CI](https://github.com/lcrawlab/recal/actions/workflows/docker-image.yml/badge.svg)](https://github.com/lcrawlab/recall/actions/workflows/docker-image.yml)

## Introduction

Standard single-cell RNA-sequencing (scRNA-seq) pipelines nearly always include unsupervised clustering as a key step in identifying biologically distinct cell types. A follow-up step in these pipelines is to test for differential expression between the identified clusters. When algorithms over-cluster, downstream analyses will produce inflated P-values resulting in increased false discoveries.
Here, we present `recall` (Calibrated Clustering with Artificial Variables): a new method for protecting against over-clustering by controlling for the impact of double-dipping.
Importantly, our approach can be applied to any clustering algorithm (implemented here are the Louvain and Leiden algorithms with plans to implement the K-means and hierarchical clustering algorithms).
`recall` provides state-of-the-art clustering performance and can rapidly analyze large-scale scRNA-seq studies, even on a personal laptop.

## Installation

You can install the lastest development version by using the [devtools](https://CRAN.R-project.org/package=devtools) library. To install this package with devtools, use this command:

```r
devtools::install_github("lcrawlab/recall")
```

Although it is not explicitly a dependency, making sure you have `presto` installed will make `recall` much faster.

```r
devtools::install_github("immunogenomics/presto")
```


## Tutorial

```r
library(Seurat)
library(SeuratData)

library(recall)

set.seed(123)

# load pbmc3k dataset
SeuratData::InstallData("pbmc3k")
data("pbmc3k")

pbmc3k <- UpdateSeuratObject(pbmc3k)

pbmc3k <- NormalizeData(pbmc3k)
pbmc3k <- FindVariableFeatures(pbmc3k)
pbmc3k <- ScaleData(pbmc3k)
pbmc3k <- RunPCA(pbmc3k)
pbmc3k <- FindNeighbors(pbmc3k)
pbmc3k <- RunUMAP(pbmc3k, dims = 1:10)

pbmc_default <- FindClusters(pbmc3k)
pbmc_recall <- FindClustersRecall(pbmc3k)

DimPlot(pbmc_default) + DimPlot(pbmc_recall)
```
## Overview of the Method

The `recall` algorithm consists of three simple steps:

1. First, we generate synthetic null variables, inspired by knockoff variables (Barber and Candès,2015) , where we augment the single-cell data being analyzed with "fake" genes that are known not to contribute to any unique cell type.
2. Second, we perform both preprocessing and clustering on this augmented dataset.
3. Third, we calibrate the number of inferred clusters by using a hypothesis testing strategy with a data-dependent threshold to determine if there is a statistically significant difference between groups. If any pair of groups does not have statistically significant differences then re-clustering occurs.

The synthetic genes act as negative control variables; they go through the same analytic steps as the real data and are presented with the same opportunity to be identified as marker genes.
The `recall` algorithm uses the guiding principle that well-calibrated clusters (i.e., those representing real groups) should have significantly differentially expressed genes after correcting for multiple hypothesis tests, while over-clustered groups will not.
We use this rule to iteratively re-cluster cells until the inferred clusters are well-calibrated and the observed differences in expression between groups are not due to the effects of double-dipping.

## Relevant Citations
`recall` is currently on the bioRxiv, [here](https://www.biorxiv.org/content/10.1101/2024.03.08.584180v1).

A. DenAdel, M. Ramseier, A. Navia, A. Shalek, S. Raghavan, P. Winter, A. Amini, and L. Crawford. A knockoff calibration method to avoid over-clustering in single-cell RNA-sequencing. _bioRxiv_.

## Questions and Feedback
For questions or concerns with `recall`, please contact
[Alan DenAdel](mailto:alan_denadel@brown.edu) or [Lorin Crawford](lcrawford@microsoft.com). Any feedback on the software, manuscript, and tutorials is appreciated.
