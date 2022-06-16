
# babelmixr2

<!-- badges: start -->
<!-- badges: end -->

The goal of babelmixr2 is to convert nlmixr2 syntax to other commonly
used tools.

## Installation

If we decide to submit to CRAN, you can install the released version of
babelmixr2 from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("babelmixr2")
```

Otherwise you can always install from github:

```r
remotes::install_github("nlmixr2/babelmixr2")
```

## Example

After installed, if you use the standard interface, you can convert to monolix with

```r
mod <- nlmixr(nlmixrFun, nlmmixrData, est="monolix")
```

