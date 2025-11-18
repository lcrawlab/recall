# Runs a typical Seurat workflow on a Seurat object (up to dimensionality reduction and clustering).

Given a Seurat object, returns a new Seurat that has been normalized,
had variable features identified, scaled, had principal components
computed, had clusters identified, and had tSNE and UMAP embeddings
determined.

## Usage

``` r
seurat_workflow(
  seurat_obj,
  num_variable_features,
  resolution_param = 0.8,
  visualization_method = "umap",
  num_dims = 10,
  algorithm = "louvain"
)
```

## Arguments

- seurat_obj:

  A Seurat object that will be analyzed.

- num_variable_features:

  The number of variable features to use in the analysis.

- resolution_param:

  The resolution parameter to use when clustering.

- visualization_method:

  Either "umap" or "tsne".

- num_dims:

  The number of principal components to use.

- algorithm:

  The clustering algorithm to use, either "louvain" or "leiden".

## Value

A Seurat object containing the relevant analysis results.
