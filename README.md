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

## Monolix Setup

While not required, you can get/install the R 'lixoftConnectors' package in the
'Monolix' installation, as described at the following url
<https://monolix.lixoft.com/monolix-api/lixoftconnectors_installation/>. When
'lixoftConnectors' is available, R can run 'Monolix' directly instead of using a
command line.

## Example

After installed, if you use the standard interface, you can convert to Monolix with

```r
mod <- nlmixr(nlmixrFun, nlmmixrData, est="monolix")
```

or, you can convert to NONMEM with

```r
mod <- nlmixr(nlmixrFun, nlmmixrData, est="nonmem")
```
