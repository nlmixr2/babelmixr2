# Determine standardized rxode2 column names from data

Determine standardized rxode2 column names from data

## Usage

``` r
getStandardColNames(data)
```

## Arguments

- data:

  A data.frame as the source for column names

## Value

A named character vector where the names are the standardized names and
the values are either the name of the column from the data or `NA` if
the column is not present in the data.

## Examples

``` r
getStandardColNames(data.frame(ID=1, DV=2, Time=3, CmT=4))
#>     id   time    amt   rate    dur   evid    cmt     ss     ii   addl     dv 
#>   "ID" "Time"     NA     NA     NA     NA  "CmT"     NA     NA     NA   "DV" 
#>    mdv   dvid   cens  limit 
#>     NA     NA     NA     NA 
```
