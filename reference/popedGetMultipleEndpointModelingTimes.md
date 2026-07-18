# Get Multiple Endpoint Modeling Times

This function takes a vector of times and a corresponding vector of IDs,
groups the times by their IDs, initializes an internal C++ global
TimeIndexer, that is used to efficiently lookup the final output from
the rxode2 solve and then returns the sorted unique times.

The `popedMultipleEndpointIndexDataFrame()` function can be used to
visualize the internal data structure inside R, but it does not show all
the indexes in the case of time ties for a given ID. Rather it shows one
of the indexs and the total number of indexes in the data.frame

## Usage

``` r
popedGetMultipleEndpointModelingTimes(times, modelSwitch, sorted = FALSE)

popedMultipleEndpointIndexDataFrame(print = FALSE)
```

## Arguments

- times:

  A numeric vector of times.

- modelSwitch:

  An integer vector of model switch indicator corresponding to the times

- sorted:

  A boolean indicating if the returned times should be sorted

- print:

  boolean for `popedMultipleEndpointIndexDataFrame()` when `TRUE` show
  each id/index per time even though it may not reflect in the returned
  data.frame

## Value

A numeric vector of unique times.

## Examples

``` r


# \donttest{

times <- c(1.1, 1.2, 1.3, 2.1, 2.2, 3.1)
modelSwitch <- c(1, 1, 1, 2, 2, 3)
sortedTimes <- popedGetMultipleEndpointModelingTimes(times, modelSwitch, TRUE)
print(sortedTimes)
#> [1] 1.1 1.2 1.3 2.1 2.2 3.1

# now show the output of the data frame representing the model
# switch to endpoint index

popedMultipleEndpointIndexDataFrame()
#>   time MS:1 N:1 MS:3 N:3 MS:5 N:5
#> 1  1.1    1   1   NA  NA   NA  NA
#> 2  1.2    2   1   NA  NA   NA  NA
#> 3  1.3    3   1   NA  NA   NA  NA
#> 4  2.1   NA  NA    4   1   NA  NA
#> 5  2.2   NA  NA    5   1   NA  NA
#> 6  3.1   NA  NA   NA  NA    6   1

# now show a more complex example with overlaps etc.

times <- c(1.1, 1.2, 1.3, 0.5, 2.2, 1.1, 0.75,0.75)
modelSwitch <- c(1, 1, 1, 2, 2, 2, 3, 3)
sortedTimes <- popedGetMultipleEndpointModelingTimes(times, modelSwitch, TRUE)
print(sortedTimes)
#> [1] 0.50 0.75 1.10 1.20 1.30 2.20

popedMultipleEndpointIndexDataFrame(TRUE) # Print to show individual matching
#> modelSwitch: 2 time: 0.500000: 4
#> modelSwitch: 3 time: 0.750000: 7, 8
#> modelSwitch: 2 time: 1.100000: 6
#> modelSwitch: 1 time: 1.100000: 1
#> modelSwitch: 1 time: 1.200000: 2
#> modelSwitch: 1 time: 1.300000: 3
#> modelSwitch: 2 time: 2.200000: 5
#>   time MS:1 N:1 MS:3 N:3 MS:5 N:5
#> 1 0.50   NA  NA    4   1   NA  NA
#> 2 0.75   NA  NA   NA  NA    7   2
#> 3 1.10    1   1    6   1   NA  NA
#> 4 1.20    2   1   NA  NA   NA  NA
#> 5 1.30    3   1   NA  NA   NA  NA
#> 6 2.20   NA  NA    5   1   NA  NA

# }
```
