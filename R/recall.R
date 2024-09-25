
#' @title Returns a Seurat object that contains additional (fake) RNA
#' expression counts.
#'
#' @description Given a Seurat object, returns a new Seurat object whose RNA
#' expression counts includes the
#' variable features from the original object and an equal number of artificial
#' features.
#'
#' @param seurat_obj A Seurat object containing RNA expression counts.
#' @param assay The assay to generate artificial variables from.
#' @param null_method The generating distribution for the synthetic null variables (ZIP, NB, ZIP-copula, NB-copula)
#' @param cores The number of cores to use in generating synthetic null variables.
#' @param verbose Whether or not to show logging.
#' @returns A Seurat object that contains the original variable features and an
#' equal number of artificial features.
#' @name get_seurat_obj_with_artificial_variables
get_seurat_obj_with_artificial_variables <- function(seurat_obj, assay = "RNA", null_method = "ZIP", verbose = TRUE, cores) {

  if (verbose) {
    message("Pulling data from Seurat object")
  }

  var_features <- Seurat::VariableFeatures(seurat_obj)
  seurat_obj_data <- as.data.frame(t(as.matrix(Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts")[var_features, ])))

  #if (verbose) {
  #  message("Estimating the distribution of each gene")
  #}
  if (verbose) {
    message("Computing artificial features")
  }

  if (null_method == "ZIP") {
    estimates <- lapply(seurat_obj_data, estimate_zi_poisson)
    sampling_function <- function(x) {
                                      rzipoisson(nrow(seurat_obj_data),
                                                 x$lambda.hat,
                                                 x$pi.hat)
                                    }
    ko <- as.data.frame(lapply(estimates, sampling_function))

  }
  else if (null_method == "NB") {
    estimates <- lapply(seurat_obj_data, estimate_negative_binomial)
    sampling_function <- function(x) {
                                      stats::rnbinom(nrow(seurat_obj_data),
                                                 size = x$size,
                                                 mu = x$mu)
                                    }
    ko <- as.data.frame(lapply(estimates, sampling_function))
  }
  else if (null_method == "ZIP-copula") {
    ko <- estimate_zi_poisson_copula(seurat_obj_data, cores)
  }
  else if (null_method == "NB-copula") {
    ko <- estimate_negative_binomial_copula(seurat_obj_data, cores)
  }
  else if (null_method == "Poisson-copula") {
    ko <- estimate_poisson_copula(seurat_obj_data, cores)
  }
  else if (null_method == "Gaussian-copula") {
    ko <- estimate_gaussian_copula(seurat_obj_data, cores)
  }
  else {
    stop("You selected a null_method that is not supported. Choose from: ZIP, NB, ZIP-copula, NB-copula, Poisson-copula, Gaussian-copula.")
  }


  num_variable_features <- length(var_features)
  colnames(ko) <- paste0(rep("knockoff", num_variable_features), 1:num_variable_features)
  combined_data <- cbind(seurat_obj_data, ko)

  # sparsify augmented data matrix and transpose for use in Seurat
  combined_data <- Matrix::Matrix(t(combined_data), sparse = TRUE)

  new_project_name <- paste0(seurat_obj@project.name, "_with_knockoffs")
  new_seurat_obj <- Seurat::CreateSeuratObject(counts = combined_data, project = new_project_name)

  return(new_seurat_obj)
}






#' @title Returns the genes selected by the knockoff filter
#'
#' @description Given two Seurat objects, returns the  the genes selected by
#' the knockoff filter and their W statistics.
#'
#' @param seurat_obj A Seurat object
#' @param cluster1 The Idents of the cluster of interest in seurat_obj1
#' @param cluster2 The Idents of the cluster of interest in seurat_obj2
#' @param q The desired rate to control the FDR at
#' @param return_all Determines if the returned object will contain all genes
#' or just the selected genes.
#' @param num_cores The number of cores for computing marker genes in parallel.
#' @param shared_memory_max The maximum size for shared global variables.
#' @returns todo
#' @name compute_knockoff_filter
compute_knockoff_filter <- function(seurat_obj,
                                    cluster1,
                                    cluster2,
                                    q,
                                    return_all = FALSE,
                                    num_cores = 1,
                                    shared_memory_max) {
  options(future.globals.maxSize = shared_memory_max)
  # todo note what this is for, figure this out as a parameter or programmatically
  future::plan("multicore", workers = as.numeric(num_cores))

  markers <- Seurat::FindMarkers(seurat_obj,
                         ident.1 = cluster1,
                         ident.2 = cluster2,
                         logfc.threshold = 0,
                         min.pct = 0)


  # FindMarkers orders by p-value, so we can't rely on position to know which genes are which
  knockoff_indices <- grepl("^knockoff", rownames(markers))
  original_indices <- !knockoff_indices

  # subset the markers data.frame into originals and knockoffs
  knockoff_markers <- markers[knockoff_indices, ]
  original_markers <- markers[original_indices, ]

  all_genes <- rownames(seurat_obj)

  # get indices of knockoffs and originals from seurat_obj, should be [FALSE, ..., FALSE, TRUE, ..., TRUE]
  knockoff_indices_sorted <- grepl("^knockoff", all_genes)
  original_indices_sorted <- !knockoff_indices_sorted

  knockoff_names_sorted <- all_genes[knockoff_indices_sorted]
  original_names_sorted <- all_genes[original_indices_sorted]

  # sort markers data.frames by their original orderings
  knockoff_markers_sorted <- knockoff_markers[knockoff_names_sorted, ]
  original_markers_sorted <- original_markers[original_names_sorted, ]

  original_p_values <- original_markers_sorted$p_val
  knockoff_p_values <- knockoff_markers_sorted$p_val

  log_original_p_values <- -log10(original_p_values)
  log_knockoff_p_values <- -log10(knockoff_p_values)

  W <- log_original_p_values - log_knockoff_p_values

  thres <- knockoff::knockoff.threshold(W, fdr = q, offset = 1)


  if (return_all) {
    all_features <- as.data.frame(list("gene" = original_names_sorted, "W" = W))

    ret <-  list("all_features" = all_features, "threshold" = thres)

    return(ret)
  }
  selected_indices <- which(W >= thres) # todo check if this should be > (case where threshold is Inf, but there are still some Inf -log p)
  #selected_indices <- which(W > thres) # todo check if this should be > (case where threshold is Inf, but there are still some Inf -log p)

  selected_genes <- original_names_sorted[selected_indices]
  selected_Ws <- W[selected_indices]

  selected_features <- as.data.frame(list("selected_gene" = selected_genes, "W" = selected_Ws))

  selected_features <- selected_features[order(selected_features$W, decreasing = TRUE), ]

  ret <-  list("selected_features" = selected_features, "threshold" = thres)

  return(ret)
}








#' @title Runs a typical Seurat workflow on a Seurat object (up to
#' dimensionality reduction and clustering).
#'
#' @description Given a Seurat object, returns a new Seurat that has been
#' normalized, had variable features identified, scaled, had principal
#' components computed, hadclusters identified, and had tSNE and UMAP
#' embeddings determined.
#'
#' @param seurat_obj The Seurat object that will be analyzed.
#' @param resolution_start The starting resolution to be used for the
#' clustering algorithm (Louvain and Leiden algorithms).
#' @param num_clusters_start The starting number of clusters to be used for the
#' clustering algorithm (K-means and Hierarchical clustering algorithms).
#' @param reduction_percentage The amount that the starting parameter will be
#' reduced by after each iteration (between 0 and 1).
#' @param dims The dimensions to use as input features (i.e. 1:10).
#' @param algorithm The clustering algorithm to be used.
#' @param null_method The generating distribution for the synthetic null variables (ZIP, NB, ZIP-copula, NB-copula)
#' @param assay The assay to generate artificial variables from.
#' @param cores The number of cores to compute marker genes in parallel.
#' @param shared_memory_max The maximum size for shared global variables.
#' Increased this variable if you see the following error:
#' The total size of the X globals that need to be exported for the future expression 
#' ('FUN()') is X GiB. This exceeds the maximum allowed size of 500.00 MiB 
#' (option 'future.globals.maxSize'). The X largest globals are ... 
#' @param verbose Whether or not to show all logging.
#' @returns Returns a Seurat object where the idents have been updated with the
#' clusters determined via the recall algorithm.
#' Latest clustering results will be stored in the object metadata under
#' recall_clusters'. Note that 'recall_clusters' will be overwritten ever
#' time FindClustersRecall is run.
#' @name FindClustersRecall
#' @export
FindClustersRecall <- function(seurat_obj,
                                 resolution_start = 0.8,
                                 reduction_percentage = 0.2,
                                 num_clusters_start = 20,
                                 dims = 1:10,
                                 algorithm = "louvain", # todo implement all algos
                                 null_method = "ZIP",
                                 assay = "RNA",
                                 cores = 1,
                                 shared_memory_max = 8000 * 1024^2,
                                 verbose = TRUE) {

  # todo check function arguments for validity

  augmented_with_artificial_variables_seurat_obj <- get_seurat_obj_with_artificial_variables(seurat_obj,
                                                       assay=assay,
                                                       null_method=null_method,
                                                       verbose=verbose,
                                                       cores=cores)

  num_variable_features <- 2 * length(Seurat::VariableFeatures(seurat_obj))


  # Pre-process data
  options(future.globals.maxSize = shared_memory_max)
  # todo log number of cores being used
  future::plan("multicore", workers = as.numeric(cores))
  #options(future.globals.maxSize = 8000 * 1024^2)

    if (verbose) {
      message(paste("Number of cores:", cores))
    }

  #plan("multicore", workers = as.numeric(cores))

  augmented_with_artificial_variables_seurat_obj <- Seurat::NormalizeData(augmented_with_artificial_variables_seurat_obj,
                                               verbose = FALSE)
   
  augmented_with_artificial_variables_seurat_obj <- Seurat::FindVariableFeatures(augmented_with_artificial_variables_seurat_obj,
                                                      selection.method = "vst",
                                                      nfeatures = num_variable_features,
                                                      verbose = FALSE)

  augmented_with_artificial_variables_seurat_obj <- Seurat::ScaleData(augmented_with_artificial_variables_seurat_obj, verbose = FALSE)
  augmented_with_artificial_variables_seurat_obj <- Seurat::RunPCA(augmented_with_artificial_variables_seurat_obj,
                                        features = Seurat::VariableFeatures(object = augmented_with_artificial_variables_seurat_obj),
                                        verbose = FALSE)

  augmented_with_artificial_variables_seurat_obj <- Seurat::FindNeighbors(augmented_with_artificial_variables_seurat_obj,
                                               dims = dims,
                                               verbose = FALSE)

  resolution_param <- resolution_start


  first_iteration <- TRUE

  while (TRUE) {
    if (verbose) {
      message("####################################################################")
      message(paste("Finding clusters with", stringr::str_to_title(algorithm), "algorithm"))
      message(paste("Resolution param:", resolution_param))
    }

    if (algorithm == "louvain") {
      augmented_with_artificial_variables_seurat_obj <- Seurat::FindClusters(augmented_with_artificial_variables_seurat_obj,
                                                  resolution = resolution_param,
                                                  verbose = FALSE)
    }

    if (algorithm == "leiden") {
      #plan("sequential") # todo log number of cores being used # this is a weird one because leiden has a forked job hanging
      augmented_with_artificial_variables_seurat_obj <- Seurat::FindClusters(augmented_with_artificial_variables_seurat_obj,
                                                  resolution = resolution_param,
                                                  algorithm = 4,
                                                  method = "igraph",
                                                  verbose = FALSE)
    }

    # Reduce resolution for next iteration of the loop
    resolution_param <- (1 - reduction_percentage) * resolution_param

    k <- length(levels(Seurat::Idents(augmented_with_artificial_variables_seurat_obj)))
    #knock_idents <- 0:(k-1)

    if (verbose) {
      message("Num clusters:")
      message(k)
    }

    knock_idents <- levels(Seurat::Idents(augmented_with_artificial_variables_seurat_obj))

    num_selected_matrix <- matrix(nrow = k, ncol = k)

    found_no_sign_diff <- FALSE

    num_clusters <- length(knock_idents)


    if (verbose) {
      progress_bar_length <- num_clusters * (num_clusters - 1) / 2
      cli::cli_progress_bar("Processing cluster pairs:",
                            total = progress_bar_length,
                            clear = FALSE)
    }

    m <- 0
    for (i in 1:num_clusters) {
      for (j in 1:num_clusters) {
        if (j >= i) {
          next
        }

        m <- m + 1

        if (verbose) {
          cli::cli_progress_update()
        }

        markers_selected <- compute_knockoff_filter(seurat_obj = augmented_with_artificial_variables_seurat_obj,
                                                    cluster1 = knock_idents[i],
                                                    cluster2 = knock_idents[j],
                                                    q = 0.05,
                                                    num_cores = cores,
                                                    shared_memory_max = shared_memory_max)

        num_selected <- nrow(markers_selected$selected_features)

        if (num_selected == 0) {
          found_no_sign_diff <- TRUE
          break
        }

        num_selected_matrix[i, j] <- num_selected
        num_selected_matrix[j, i] <- num_selected

      }
      if (found_no_sign_diff) {
        if (verbose) {
          cli::cli_progress_done()
          message("Found clusters with no significant differences.")
          message("Progressing to next clustering iteration.")
        }
        first_iteration <- FALSE
        break
      }
    }

    if (found_no_sign_diff) {
      next
    }
    break
  }

  if (first_iteration) {
    warning("Only a single iteration occurred. The inferred cluster labels may be underclustered. To prevent this, you may want to re-run recall with a larger starting parameter.")
  }

  seurat_obj@meta.data$recall_clusters <- Seurat::Idents(augmented_with_artificial_variables_seurat_obj)
  Seurat::Idents(seurat_obj) <- seurat_obj@meta.data$recall_clusters

  return(seurat_obj)
}





#' @title Runs a typical Seurat workflow on a Seurat object (up to
#' dimensionality reduction and clustering).
#'
#' @description Given a Seurat object, returns a new Seurat that has been
#' normalized, had variable features identified, scaled, had principal
#' components computed, hadclusters identified, and had tSNE and UMAP
#' embeddings determined.
#'
#' @param seurat_obj The Seurat object that will be analyzed.
#' @param resolution_start The starting resolution to be used for the
#' clustering algorithm (Louvain and Leiden algorithms).
#' @param num_clusters_start The starting number of clusters to be used for the
#' clustering algorithm (K-means and Hierarchical clustering algorithms).
#' @param reduction_percentage The amount that the starting parameter will be
#' reduced by after each iteration (between 0 and 1).
#' @param dims The dimensions to use as input features (i.e. 1:10).
#' @param algorithm The clustering algorithm to be used.
#' @param null_method The generating distribution for the synthetic null variables (ZIP, NB, ZIP-copula, NB-copula)
#' @param assay The assay to generate artificial variables from.
#' @param cores The number of cores to compute marker genes in parallel.
#' @param shared_memory_max The maximum size for shared global variables.
#' Increased this variable if you see the following error:
#' The total size of the X globals that need to be exported for the future expression 
#' ('FUN()') is X GiB. This exceeds the maximum allowed size of 500.00 MiB 
#' (option 'future.globals.maxSize'). The X largest globals are ... 
#' @param verbose Whether or not to show all logging.
#' @returns Returns a Seurat object where the idents have been updated with the
#' clusters determined via the countsplit algorithm.
#' Latest clustering results will be stored in the object metadata under
#' countsplit_clusters'. Note that 'countsplit_clusters' will be overwritten ever
#' time FindClustersCountsplit is run.
#' @name FindClustersCountsplit
#' @export
FindClustersCountsplit <- function(seurat_obj,
                                   resolution_start = 0.8,
                                   reduction_percentage = 0.2,
                                   num_clusters_start = 20,
                                   dims = 1:10,
                                   algorithm = "louvain", # todo implement all algos
                                   null_method = "ZIP",
                                   assay = "RNA",
                                   cores = 1,
                                   shared_memory_max = 8000 * 1024^2,
                                   verbose = TRUE) {

  options(future.globals.maxSize = shared_memory_max)
  # todo log number of cores being used
  future::plan("multicore", workers = as.numeric(cores))

  num_variable_features <- length(Seurat::VariableFeatures(seurat_obj))


  # follow this issue
  #https://github.com/anna-neufeld/countsplit/issues/8

  # todo only do this for variable features
  split <- countsplit::countsplit(Seurat::GetAssayData(seurat_obj, assay))
  Xtrain <- split[[1]]
  Xtest <- split[[2]]

  seurat_obj_train <- Seurat::CreateSeuratObject(counts = Xtrain)
  seurat_obj_test <- Seurat::CreateSeuratObject(counts = Xtest)


  # process training data
  seurat_obj_train <- Seurat::NormalizeData(seurat_obj_train,
                                            verbose = FALSE)

  seurat_obj_train <- Seurat::FindVariableFeatures(seurat_obj_train,
                                                   selection.method = "vst",
                                                   nfeatures = num_variable_features,
                                                   verbose = FALSE)

  seurat_obj_train <- Seurat::ScaleData(seurat_obj_train, verbose = FALSE)
  seurat_obj_train <- Seurat::RunPCA(seurat_obj_train,
                                     features = Seurat::VariableFeatures(object = seurat_obj_train),
                                     verbose = FALSE)

  seurat_obj_train <- Seurat::FindNeighbors(seurat_obj_train,
                                            dims = dims,
                                            verbose = FALSE)

  # process test data (no need for PCA and FindNeighbors since we are just assigning idents based on training idents)
  seurat_obj_test <- Seurat::NormalizeData(seurat_obj_test,
                                           verbose = FALSE)

  seurat_obj_test <- Seurat::FindVariableFeatures(seurat_obj_test,
                                                  selection.method = "vst",
                                                  nfeatures = num_variable_features,
                                                  verbose = FALSE)

  seurat_obj_test <- Seurat::ScaleData(seurat_obj_test, verbose = FALSE)



  resolution_param <- resolution_start

  # set up multicore for FindMarkers
  future::plan("multicore", workers = as.numeric(cores))

  first_iteration <- TRUE

  while (TRUE) {
    if (verbose) {
      message("####################################################################")
      message(paste("Finding clusters with", stringr::str_to_title(algorithm), "algorithm"))
      message(paste("Resolution param:", resolution_param))
    }
    
    if (algorithm == "louvain") {
      seurat_obj_train <- Seurat::FindClusters(seurat_obj_train,
                                                  resolution = resolution_param,
                                                  verbose = FALSE)
    }
    
    if (algorithm == "leiden") {
      #plan("sequential") # todo log number of cores being used # this is a weird one because leiden has a forked job hanging
      seurat_obj_train <- Seurat::FindClusters(seurat_obj_train,
                                                  resolution = resolution_param,
                                                  algorithm = 4,
                                                  method = "igraph",
                                                  verbose = FALSE)
    }
    
    # Reduce resolution for next iteration of the loop
    resolution_param <- (1 - reduction_percentage) * resolution_param
    
    Seurat::Idents(seurat_obj_test) <- Seurat::Idents(seurat_obj_train)
    
    k <- length(levels(Seurat::Idents(seurat_obj_test)))
    #knock_idents <- 0:(k-1)
    
    if (verbose) {
      message("Num clusters:")
      message(k)
    }
    
    countsplit_idents <- levels(Seurat::Idents(seurat_obj_test))
    
    num_selected_matrix <- matrix(nrow = k, ncol = k)
    
    found_no_sign_diff <- FALSE
    
    num_clusters <- length(countsplit_idents)
    
    
    if (verbose) {
      progress_bar_length <- num_clusters * (num_clusters - 1) / 2
      cli::cli_progress_bar("Processing cluster pairs:",
                            total = progress_bar_length,
                            clear = FALSE)
    }
    
    m <- 0
    for (i in 1:num_clusters) {
      for (j in 1:num_clusters) {
        if (j >= i) {
          next
        }
        
        m <- m + 1
        
        if (verbose) {
          cli::cli_progress_update()
        }
        
        markers_selected <- Seurat::FindMarkers(seurat_obj_test,
                                                ident.1 = countsplit_idents[i],
                                                ident.2 = countsplit_idents[j])
    
        num_selected <- sum(markers_selected$p_val_adj < 0.05)
        
        if (num_selected == 0) {
          found_no_sign_diff <- TRUE
          break
        }
        
        num_selected_matrix[i, j] <- num_selected
        num_selected_matrix[j, i] <- num_selected
        
      }
      if (found_no_sign_diff) {
        if (verbose) {
          cli::cli_progress_done()
          message("Found clusters with no significant differences.")
          message("Progressing to next clustering iteration.")
        }
        first_iteration <- FALSE
        break
      }
    }
    
    if (found_no_sign_diff) {
      next
    }
    break
  }

  if (first_iteration) {
    warning("Only a single iteration occurred. The inferred cluster labels may be underclustered. To prevent this, you may want to re-run FindClustersCountsplit with a larger starting parameter.")
  }


  seurat_obj@meta.data$countsplit_clusters <- Seurat::Idents(seurat_obj_test)
  Seurat::Idents(seurat_obj) <- seurat_obj@meta.data$countsplit_clusters

  return(seurat_obj)

}