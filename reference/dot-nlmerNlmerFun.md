# Surrogate gradient function for lme4::nlmer

Computes predictions and the analytical gradient for every observation
using the nlm C machinery (a single multithreaded
[`nlmixr2est::nlmerSolveGrad()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmerSolveGrad.html)
call against the model kept resident by
[`nlmixr2est::.nlmSetupEnv()`](https://nlmixr2.github.io/nlmixr2est/reference/dot-nlmSetupEnv.html)).
Parameters arrive in positional order matching
`.nlmerGlobal$nlmerEnv$paramNames` (`THETA[i]` order).

## Usage

``` r
.nlmerNlmerFun(TIME, ...)
```

## Arguments

- TIME:

  Observation times vector (unused; kept for formula compatibility – the
  solver uses the resident per-subject event data)

- ...:

  Individual-level parameter vectors (one per param in
  `saemParamsToEstimate`)

## Value

Numeric vector of predictions with `"gradient"` attribute

## Details

lme4 passes one value per observation for each nonlinear parameter;
within a subject these are constant (phi = beta + b), so the per-subject
parameter vector is read from each subject's first observation row.
