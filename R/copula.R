simulate_data_scDesign3 <- function(data_matrix, cores, family) {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = data_matrix))
    SummarizedExperiment::colData(sce)$cell_type <- "1" # scDesign3 needs a cell type so we just make it the same for all cells

    simulated_data <- scDesign3::scdesign3(sce,
                                           celltype = "cell_type",
                                           pseudotime = NULL,
                                           spatial = NULL,
                                           other_covariates = NULL,
                                           empirical_quantile = FALSE,
                                           usebam=TRUE, # to speedup marginal inference
                                           mu_formula = "1",
                                           sigma_formula = "1",
                                           corr_formula = "1",
                                           family_use = family, # this is the key parameter
                                           nonzerovar = FALSE,
                                           n_cores = cores,
                                           parallelization = "mcmapply",
                                           important_feature = "all",
                                           nonnegative = FALSE,
                                           copula = "gaussian",
                                           fastmvn = TRUE)

    ko <- simulated_data$new_count

    return(ko)
}


#' @title todo
#'
#' @description Given data, computes todo
#'
#' @param data_matrix The data to estimate parameters from.
#' @param cores The number of CPU cores to use in estimation by scDesign3.
#' @returns todo
#' @name estimate_negative_binomial_copula
estimate_zi_poisson_copula <- function(data_matrix, cores) {
    family <- "zip"
    ko <- simulate_data_scDesign3(data_matrix, cores, family)
    return(ko)
}


#' @title todo
#'
#' @description Given data, computes todo
#'
#' @param data_matrix The data to estimate parameters from.
#' @param cores The number of CPU cores to use in estimation by scDesign3.
#' @returns todo
#' @name estimate_negative_binomial_copula
estimate_negative_binomial_copula <- function(data_matrix, cores) {
    family <- "nb"
    ko <- simulate_data_scDesign3(data_matrix, cores, family)
    return(ko)
}


#' @title todo
#'
#' @description Given data, computes todo
#'
#' @param data_matrix The data to estimate parameters from.
#' @param cores The number of CPU cores to use in estimation by scDesign3.
#' @returns todo
#' @name estimate_negative_binomial_copula
estimate_poisson_copula <- function(data_matrix, cores) {
    family <- "poisson"
    ko <- simulate_data_scDesign3(data_matrix, cores, family)
    return(ko)
}




#' @title todo
#'
#' @description Given data, computes todo
#'
#' @param data_matrix The data to estimate parameters from.
#' @param cores The number of CPU cores to use in estimation by scDesign3.
#' @returns todo
#' @name estimate_negative_binomial_copula
estimate_gaussian_copula <- function(data_matrix, cores) {
    family <- "gaussian"
    ko <- simulate_data_scDesign3(data_matrix, cores, family)
    return(ko)
}



