# Solve poped problem for appropriate times with single/multiple endpoint models

This really should not be called directly (if not setup correctly can
crash R)

## Usage

``` r
.popedSolveIdME(theta, umt, mt, ms, nend, id, totn)

.popedSolveIdME2(theta, umt, mt, ms, nend, id, totn)
```

## Arguments

- theta:

  parameters (includes covariates and modeling times)

- umt:

  unique times sampled

- mt:

  original unsorted time (to match the f/w against)

- ms:

  model switch parameter integer starting with 1 (related to dvid in
  rxode2)

- nend:

  specifies the number of endpoints in this model

- id:

  this is the design identifier

- totn:

  This is the total number of design points tested

## Value

a data frame with \$f and \$w corresponding to the function value and
standard deviation at the sampling point

## Author

Matthew L. Fidler
