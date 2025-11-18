# Random data generation for the zero-infalted Poisson distribution with Poisson parameter lambda and zero proportion prop.zero.

Given the number of samples desired, a Poisson parameter, lambda, and a
zero proportion, prop.zero, simulates the number of desired samples from
ZIP(lambda, prop.zero).

## Usage

``` r
rzipoisson(n, lambda, prop.zero)
```

## Arguments

- n:

  The number of samples to be simulated.

- lambda:

  The Poisson rate parameter.

- prop.zero:

  The proportion of excess zeroes.

## Value

Simulated data from ZIP(lambda, prop.zero).
