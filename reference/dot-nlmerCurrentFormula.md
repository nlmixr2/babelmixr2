# Accessor for the nlmer formula currently being fit

[`lme4::nlmer()`](https://rdrr.io/pkg/lme4/man/nlmer.html) captures
`match.call()$formula` and re-evaluates it deep in its own call stack,
where the babelmixr2 fit locals are not in scope, so the formula must be
referenced through an exported, namespace-qualified call rather than a
`:::` reference to the package's internal `.nlmerGlobal` (which CRAN
disallows). This returns the 3-part formula stashed by
`.nlmerFitModel()` for the in-progress fit.

## Usage

``` r
.nlmerCurrentFormula()
```

## Value

The 3-part nlmer formula for the active fit, or `NULL` when no fit is in
progress.
