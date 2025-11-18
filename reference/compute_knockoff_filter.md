# Returns the genes selected by the knockoff filter

Given two Seurat objects, returns the the genes selected by the knockoff
filter and their W statistics.

## Usage

``` r
compute_knockoff_filter(
  seurat_obj,
  cluster1,
  cluster2,
  q,
  return_all = FALSE,
  num_cores = 1,
  shared_memory_max
)
```

## Arguments

- seurat_obj:

  A Seurat object

- cluster1:

  The Idents of the cluster of interest in seurat_obj1

- cluster2:

  The Idents of the cluster of interest in seurat_obj2

- q:

  The desired rate to control the FDR at

- return_all:

  Determines if the returned object will contain all genes or just the
  selected genes.

- num_cores:

  The number of cores for computing marker genes in parallel.

- shared_memory_max:

  The maximum size for shared global variables.

## Value

todo
