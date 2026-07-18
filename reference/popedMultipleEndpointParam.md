# Populates Multiple Endpoint Parameters for internal solving

This function populates a numeric vector with parameters and unique
times and also populates the internal C++ global index

## Usage

``` r
popedMultipleEndpointParam(p, times, modelSwitch, maxMT, optTime = TRUE)
```

## Arguments

- p:

  A numeric vector of parameters

- times:

  A numeric vector of times

- modelSwitch:

  An integer vector indicating model switches from PopED

- maxMT:

  An integer specifying the maximum number of time points in the mtimes
  model

## Value

A numeric vector containing the parameters followed by unique times, if
the maximum number of times is greater than the input this will append
the maximum observed times in the input. This assumes the first
parameter is the id and is dropped fro the output.

## Details

- This function first uses the input times and model switches to a
  global time indexer.

- It then creates a new numeric vector that combines the input
  parameters and unique times. If the number of times is less than
  `maxMT`, the remaining elements are filled with the maximum time.

## Author

Matthew L. Fidler

## Examples

``` r

# \donttest{

p <- c(1.0, 2.0, 3.0)
times <- c(0.5, 1.5, 2.5)
modelSwitch <- c(1, 2, 3)
maxMT <- 5
popedMultipleEndpointParam(p, times, modelSwitch, maxMT)
#> [1] 2.0 3.0 0.5 1.5 2.5 2.5 2.5

# }
```
