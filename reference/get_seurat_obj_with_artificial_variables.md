# Returns a Seurat object that contains additional (fake) RNA expression counts.

Given a Seurat object, returns a new Seurat object whose RNA expression
counts includes the variable features from the original object and an
equal number of artificial features.

## Usage

``` r
get_seurat_obj_with_artificial_variables(
  seurat_obj,
  assay = "RNA",
  null_method = "ZIP",
  verbose = TRUE,
  cores
)
```

## Arguments

- seurat_obj:

  A Seurat object containing RNA expression counts.

- assay:

  The assay to generate artificial variables from.

- null_method:

  The generating distribution for the synthetic null variables (ZIP, NB,
  ZIP-copula, NB-copula)

- verbose:

  Whether or not to show logging.

- cores:

  The number of cores to use in generating synthetic null variables.

## Value

A Seurat object that contains the original variable features and an
equal number of artificial features.
