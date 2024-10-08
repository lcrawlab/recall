% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/recall.R
\name{FindClustersCountsplit}
\alias{FindClustersCountsplit}
\title{Runs a typical Seurat workflow on a Seurat object (up to
dimensionality reduction and clustering).}
\usage{
FindClustersCountsplit(
  seurat_obj,
  resolution_start = 0.8,
  reduction_percentage = 0.2,
  num_clusters_start = 20,
  dims = 1:10,
  algorithm = "louvain",
  null_method = "ZIP",
  assay = "RNA",
  cores = 1,
  shared_memory_max = 8000 * 1024^2,
  verbose = TRUE
)
}
\arguments{
\item{seurat_obj}{The Seurat object that will be analyzed.}

\item{resolution_start}{The starting resolution to be used for the
clustering algorithm (Louvain and Leiden algorithms).}

\item{reduction_percentage}{The amount that the starting parameter will be
reduced by after each iteration (between 0 and 1).}

\item{num_clusters_start}{The starting number of clusters to be used for the
clustering algorithm (K-means and Hierarchical clustering algorithms).}

\item{dims}{The dimensions to use as input features (i.e. 1:10).}

\item{algorithm}{The clustering algorithm to be used.}

\item{null_method}{The generating distribution for the synthetic null variables (ZIP, NB, ZIP-copula, NB-copula)}

\item{assay}{The assay to generate artificial variables from.}

\item{cores}{The number of cores to compute marker genes in parallel.}

\item{shared_memory_max}{The maximum size for shared global variables.
Increased this variable if you see the following error:
The total size of the X globals that need to be exported for the future expression
('FUN()') is X GiB. This exceeds the maximum allowed size of 500.00 MiB
(option 'future.globals.maxSize'). The X largest globals are ...}

\item{verbose}{Whether or not to show all logging.}
}
\value{
Returns a Seurat object where the idents have been updated with the
clusters determined via the countsplit algorithm.
Latest clustering results will be stored in the object metadata under
countsplit_clusters'. Note that 'countsplit_clusters' will be overwritten ever
time FindClustersCountsplit is run.
}
\description{
Given a Seurat object, returns a new Seurat that has been
normalized, had variable features identified, scaled, had principal
components computed, hadclusters identified, and had tSNE and UMAP
embeddings determined.
}
