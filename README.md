
# babelmixr

<!-- badges: start -->
[![R-CMD-check](https://github.com/nlmixrdevelopment/babelmixr/workflows/R-CMD-check/badge.svg)](https://github.com/nlmixrdevelopment/babelmixr/actions)
<!-- badges: end -->

The goal of babelmixr is to convert nlmixr syntax to other commonly
used tools.

## Installation

If decide to submit to CRAN, you can install the released version of
babelmixr from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("babelmixr")
```

Otherwise you can always install from github:

```r
remotes::install_github("nlmixrdevelopment/babelmixr")
```

## Example

After installed, if you use the standard interface, you can convert to monolix with

```r
mod <- nlmixr(nlmixrFun, nlmmixrData, est="monolix")
```

