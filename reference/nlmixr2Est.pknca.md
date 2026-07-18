# Estimate starting parameters using PKNCA

Estimate starting parameters using PKNCA

## Usage

``` r
# S3 method for class 'pknca'
nlmixr2Est(env, ...)
```

## Arguments

- env:

  Environment for the nlmixr2 estimation routines.

  This needs to have:

  \- rxode2 ui object in \`\$ui\`

  \- data to fit in the estimation routine in \`\$data\`

  \- control for the estimation routine's control options in \`\$ui\`

- ...:

  Other arguments provided to \`nlmixr2Est()\` provided for flexibility
  but not currently used inside nlmixr

## Value

A model with updated starting parameters. In the model a new element
named "nca" will be available which includes the PKNCA results used for
the calculation.

## Details

Parameters are estimated as follows:

- `ka` 4 half-lives to Tmax but not higher than 3: `log(2)/(tmax/4)`

- `vc` Inverse of dose-normalized Cmax

- `cl` Estimated as the median clearance

- `vp,vp2`2- and 4-fold the `vc`, respectively by default, controlled by
  the `vpMult` and `vp2Mult` arguments to `pkncaControl`

- `q,q2` 0.5- and 0.25-fold the `cl`, respectively by default,
  controlled by the `qMult` and `q2Mult` arguments to `pkncaControl`

The bounds for the parameter estimates are set to 10% of the first
percentile and 10 times the 99th percentile. (For ka, the lower bound is
set to the lower of 10% of the first percentile or 0.03 and the upper
bound is not modified from 10 times the 99th percentile.)

Parameter estimation methods may be changed in a future version.
