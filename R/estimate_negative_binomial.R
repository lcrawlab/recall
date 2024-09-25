
#' @title Maximum likelihood estimation for the negative binomial
#' distribution.
#'
#' @description Given data, computes the maximum likelihood estimators
#' for the negative binomial distribution with parameters: size and mu.
#'
#' @param data The data to estimate parameters from.
#' @returns Maximum likelihood estimators size and mu for the negative 
#' binomial distribution
#' @param verbose Whether or not to show all logging.
#' @name estimate_negative_binomial
estimate_negative_binomial <- function(data, verbose=FALSE) {
  
  if (verbose) { message("Attempting MLE method 1") }
  mle1 <- tryCatch(
    {
      nb_fit <- MASS::fitdistr(data, "negative binomial", method = "Nelder-Mead")
      size <- nb_fit$estimate[["size"]]
      mu <- nb_fit$estimate[["mu"]]

      # check if method returned NaN or NA without throwing an error
      if (is.na(mu) || is.na(size)) { stop() }

      return.list <- list("size" = size, "mu" = mu)
      return(return.list)
    },
    error = function(cond) {
      if (verbose) { message("MLE method 1 failed with an error.") }
      NA
    },
    warning = function(cond) {
      if (verbose) { message("MLE method 2 had a warning. Warning message:\n") }
      if (verbose) { message(cond) }
      if (verbose) { message("\n") }
      NA
    }
  )
  
  if (verbose) { message("Attempting MLE method 2") }
  mle2 <- tryCatch(
    {
      nb_fit <- fitdistrplus::fitdist(data, "nbinom", method="mle")
      size <- nb_fit$estimate[["size"]]
      mu <- nb_fit$estimate[["mu"]]

      # check if method returned NaN or NA without throwing an error
      if (is.na(mu) || is.na(size)) { stop() }

      return.list <- list("size" = size, "mu" = mu)
      return(return.list)

    },
    error = function(cond) {
      if (verbose) { message("MLE method 2 failed with an error.") }
      NA
    },
    warning = function(cond) {
      if (verbose) { message("MLE method 2 had a warning. Warning message:") }
      if (verbose) { message(cond) }
      NA
    }
  )
  
  if (verbose) { message("Attempting MME") }
  mme <- tryCatch(
    {
      nb_fit <- fitdistrplus::fitdist(data, "nbinom", method="mme")
      size <- nb_fit$estimate[["size"]]
      mu <- nb_fit$estimate[["mu"]]

      # check if method returned NaN or NA without throwing an error
      if (is.na(mu) || is.na(size)) { stop() }

      return.list <- list("size" = size, "mu" = mu)
      return(return.list)
    },
    error = function(cond) {
      if (verbose) { message("MME failed with an error.") }
      NA
    },
    warning = function(cond) {
      if (verbose) { message("MME method has a warning. Warning message:") }
      if (verbose) { message(cond) }
      NA
    }
  )

  
  if (verbose) { message("Attempting MME with warnings") }
  mme <- tryCatch(
    {
      nb_fit <- fitdistrplus::fitdist(data, "nbinom", method="mme")
      size <- nb_fit$estimate[["size"]]
      mu <- nb_fit$estimate[["mu"]]

      # check if method returned NaN or NA without throwing an error
      if (is.na(mu) || is.na(size)) { stop() }

      return.list <- list("size" = size, "mu" = mu)
      return(return.list)
    },
    error = function(cond) {
      if (verbose) { message("MME failed with an error.") }
      NA
    }
  )

  if (verbose) { message("Attempting MSE") }
  mme <- tryCatch(
    {
      nb_fit <- fitdistrplus::fitdist(data, "nbinom", method="mse")
      size <- nb_fit$estimate[["size"]]
      mu <- nb_fit$estimate[["mu"]]

      # check if method returned NaN or NA without throwing an error
      if (is.na(mu) || is.na(size)) { stop() }

      return.list <- list("size" = size, "mu" = mu)
      return(return.list)
    },
    error = function(cond) {
      if (verbose) { message("MSE failed with an error.") }
      NA
    },
    warning = function(cond) {
      if (verbose) { message("MSE method failed. Warning message:") }
      if (verbose) { message(cond) }
      NA
    }
  )

  if (verbose) { message("Attempting QME") }
  mme <- tryCatch(
    {
      nb_fit <- fitdistrplus::fitdist(data, "nbinom", method="qme")
      size <- nb_fit$estimate[["size"]]
      mu <- nb_fit$estimate[["mu"]]

      # check if method returned NaN or NA without throwing an error
      if (is.na(mu) || is.na(size)) { stop() }

      return.list <- list("size" = size, "mu" = mu)
      return(return.list)
    },
    error = function(cond) {
      if (verbose) { message("QME failed with an error.") }
      NA
    },
    warning = function(cond) {
      if (verbose) { message("QME method failed. Warning message:") }
      if (verbose) { message(cond) }
      NA
    }
  )


  if (verbose) { message("Attempting MGE") }
  mme <- tryCatch(
    {
      nb_fit <- fitdistrplus::fitdist(data, "nbinom", method="mge")
      size <- nb_fit$estimate[["size"]]
      mu <- nb_fit$estimate[["mu"]]

      # check if method returned NaN or NA without throwing an error
      if (is.na(mu) || is.na(size)) { stop() }

      return.list <- list("size" = size, "mu" = mu)
      return(return.list)
    },
    error = function(cond) {
      if (verbose) { message("MGE failed with an error.") }
      NA
    },
    warning = function(cond) {
      if (verbose) { message("MGE method failed. Warning message:") }
      if (verbose) { message(cond) }
      NA
    }
  )
  


  stop("All negative binomial estimation methods failed.")  

}
