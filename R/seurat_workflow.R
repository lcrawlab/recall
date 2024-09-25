#' @title Runs a typical Seurat workflow on a Seurat object (up to
#' dimensionality reduction and clustering).
#'
#' @description Given a Seurat object, returns a new Seurat that has been
#' normalized, had variable features identified,
#' scaled, had principal components computed, had clusters identified, and had
#' tSNE and UMAP embeddings determined.
#'
#' @param seurat_obj A Seurat object that will be analyzed.
#' @param num_variable_features The number of variable features to use in the
#' analysis.
#' @param resolution_param The resolution parameter to use when clustering.
#' @param visualization_method Either "umap" or "tsne".
#' @param num_dims The number of principal components to use.
#' @param algorithm The clustering algorithm to use, either "louvain" or
#' "leiden".
#' @returns A Seurat object containing the relevant analysis results.
#' @export
#' @name seurat_workflow
seurat_workflow <- function(seurat_obj,
                            num_variable_features,
                            resolution_param = 0.8,
                            visualization_method = "umap",
                            num_dims = 10,
                            algorithm = "louvain") {
  seurat_obj <- Seurat::NormalizeData(seurat_obj)

  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj,
                                             selection.method = "vst",
                                             nfeatures = num_variable_features)

  all_genes <- rownames(seurat_obj)

  seurat_obj <- Seurat::ScaleData(seurat_obj)

  seurat_obj <- Seurat::RunPCA(seurat_obj,
                               features = Seurat::VariableFeatures(object = seurat_obj))

  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:num_dims)

  if (algorithm == "louvain") {
    seurat_obj <- Seurat::FindClusters(seurat_obj,
                                       resolution = resolution_param)
  }

  if (algorithm == "leiden") {
    seurat_obj <- Seurat::FindClusters(seurat_obj,
                                       resolution = resolution_param,
                                       algorithm = 4,
                                       method = "igraph")
  }

  if (visualization_method == "umap") {
    seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:num_dims)
  }
  if (visualization_method == "tsne") {
    seurat_obj <- Seurat::RunTSNE(seurat_obj, dims = 1:num_dims)
  }

if (visualization_method == "both") {
    seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:num_dims)
    seurat_obj <- Seurat::RunTSNE(seurat_obj, dims = 1:num_dims)
  }

  return(seurat_obj)
}
