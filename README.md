
<!-- README.md is generated from README.Rmd. Please edit that file -->

# babelmixr2

<!-- badges: start -->

[![R-CMD-check](https://github.com/nlmixr2/babelmixr2/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nlmixr2/babelmixr2/actions/workflows/R-CMD-check.yaml)
[![CRAN
version](http://www.r-pkg.org/badges/version/babelmixr2)](https://cran.r-project.org/package=babelmixr2)
[![CRAN total
downloads](https://cranlogs.r-pkg.org/badges/grand-total/babelmixr2)](https://cran.r-project.org/package=babelmixr2)
[![CRAN total
downloads](https://cranlogs.r-pkg.org/badges/babelmixr2)](https://cran.r-project.org/package=babelmixr2)
[![Codecov test
coverage](https://codecov.io/gh/nlmixr2/babelmixr2/branch/main/graph/badge.svg)](https://app.codecov.io/gh/nlmixr2/babelmixr2?branch=main)
[![CodeFactor](https://www.codefactor.io/repository/github/nlmixr2/babelmixr2/badge)](https://www.codefactor.io/repository/github/nlmixr2/babelmixr2)
![r-universe](https://nlmixr2.r-universe.dev/badges/babelmixr2)
<!-- badges: end -->

The goal of babelmixr2 is to convert nlmixr2 syntax to other commonly
used tools.

## Installation

You can install the released version of babelmixr2 from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("babelmixr2")
```

You can install from r-universe by:

``` r
# Download and install babelmixr2 in R
install.packages('babelmixr2',
                 repos = c(
                   nlmixr2 = 'https://nlmixr2.r-universe.dev',
                   CRAN = 'https://cloud.r-project.org'))
```

Otherwise you can always install from GitHub:

## What you can do with `babelmixr2`

Babelmixr2 can help you by:

  - Running your nlmixr2 model in a commercial nonlinear mixed effects
    modeling tool like
    [`NONMEM`](https://nlmixr2.github.io/babelmixr2/articles/running-nonmem.html)
    or `Monolix`

  - Convert your [`NONMEM` model to a nlmixr2
    model](https://nlmixr2.github.io/nonmem2rx/articles/convert-nlmixr2.html)
    (in conjunction with `nonmem2rx`)

  - Calculate scaling factors and automatically add initial conditions
    based on non-compartmental analysis (using `PKNCA`)

## Monolix Setup

While not required, you can get/install the R ‘lixoftConnectors’ package
in the ‘Monolix’ installation, as described at the following url
<https://monolix.lixoft.com/monolix-api/lixoftconnectors_installation/>.
When ‘lixoftConnectors’ is available, R can run ‘Monolix’ directly
instead of using a command line.

## Example

After installed, if you use the standard interface, you can convert to
Monolix with

``` r
mod <- nlmixr(nlmixrFun, nlmmixrData, est="monolix")
```

or, you can convert to NONMEM with

``` r
mod <- nlmixr(nlmixrFun, nlmmixrData, est="nonmem")
```
