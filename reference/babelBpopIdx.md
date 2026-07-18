# Get the bpop_idx by variable name for a poped database created by `babelmixr2`

This may work for other poped databases if the population parameters are
named.

## Usage

``` r
babelBpopIdx(popedInput, var)
```

## Arguments

- popedInput:

  The babelmixr2 created database

- var:

  variable to query

## Value

index of the variable

## Author

Matthew L. Fidler

## Examples

``` r

if (requireNamespace("PopED", quietly=TRUE)) {

f <- function() {
  ini({
    tV <- 72.8
    tKa <- 0.25
    tCl <- 3.75
    tF <- fix(0.9)
    pedCL <- 0.8

    eta.v ~ 0.09
    eta.ka ~ 0.09
    eta.cl ~0.25^2

    prop.sd <- fix(sqrt(0.04))
    add.sd <- fix(sqrt(5e-6))

  })
  model({
    V<-tV*exp(eta.v)
    KA<-tKa*exp(eta.ka) * (pedCL**isPediatric) # add covariate for pediatrics
    CL<-tCl*exp(eta.cl)
    Favail <- tF

    N <-  floor(t/TAU)+1
    y <- (DOSE*Favail/V)*(KA/(KA - CL/V)) *
      (exp(-CL/V * (t - (N - 1) * TAU)) *
         (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) -
         exp(-KA * (t - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))

    y ~ prop(prop.sd) + add(add.sd)
  })
}

e <- et(c( 1,8,10,240,245))

babel.db <- nlmixr2(f, e, "poped",
                    popedControl(m = 2,
                                 groupsize=20,
                                 bUseGrouped_xt=TRUE,
                                 a=list(c(DOSE=20,TAU=24,isPediatric = 0),
                                        c(DOSE=40, TAU=24,isPediatric = 0))))

babelBpopIdx(babel.db, "pedCL")

}
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#>  
#>  
#>  
#>  
#> [1] 4
```
