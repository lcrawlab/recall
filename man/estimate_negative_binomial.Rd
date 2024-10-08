% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_negative_binomial.R
\name{estimate_negative_binomial}
\alias{estimate_negative_binomial}
\title{Maximum likelihood estimation for the negative binomial
distribution.}
\usage{
estimate_negative_binomial(data, verbose = FALSE)
}
\arguments{
\item{data}{The data to estimate parameters from.}

\item{verbose}{Whether or not to show all logging.}
}
\value{
Maximum likelihood estimators size and mu for the negative
binomial distribution
}
\description{
Given data, computes the maximum likelihood estimators
for the negative binomial distribution with parameters: size and mu.
}
