
# https://en.wikipedia.org/wiki/Zero-inflated_model#Estimators_of_ZIP_parameters
# https://math.stackexchange.com/questions/2761563/maximum-likelihood-estimation-for-zero-inflated-poisson-distribution
# https://ieeexplore.ieee.org/document/9032203

#' @title Maximum likelihood estimation for the zero-inflated Poisson distribution
#' with Poisson parameter lambda and zero proportion prop.zero.
#'
#' @description Given data, computes the maximum likelihood estimators
#' for the zero-inflated Poisson distribution.
#'
#' @param data The data to estimate parameters from.
#' @returns Maximum likelihood estimators of the zero-inflated Poisson
#' distribution
#' @name estimate_zi_poisson
estimate_zi_poisson <- function(data) {
  num.zeros <- sum(data == 0)
  r0 <- 1 / length(data) * num.zeros
  
  x.bar = mean(data)
  
  gamma <- x.bar / (1 - r0)
  
  lambda.hat <- lamW::lambertW0(-gamma * exp(-gamma)) + gamma
  
  pi.hat <- 1 - x.bar / lambda.hat


  return.list <- list("lambda.hat" = lambda.hat, "pi.hat" = pi.hat)
  return(return.list)
}


#' @title Random data generation for the zero-infalted Poisson distribution
#' with Poisson parameter lambda and zero proportion prop.zero.
#'
#' @description Given the number of samples desired, a Poisson parameter,
#' lambda, and a zero proportion, prop.zero, simulates the number of desired
#' samples from ZIP(lambda, prop.zero).
#'
#' @param n The number of samples to be simulated.
#' @param lambda The Poisson rate parameter.
#' @param prop.zero The proportion of excess zeroes.
#' @returns Simulated data from ZIP(lambda, prop.zero).
#' @name rzipoisson
rzipoisson <- function(n, lambda, prop.zero) {
  data <- c()


  for (i in 1:n) {
    if (stats::runif(1) < prop.zero) {
      data[i] <- 0
    }
    else {
      data[i] <- stats::rpois(1, lambda)
    }
  }
  return(data)
} 

