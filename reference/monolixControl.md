# Monolix Controller for nlmixr2

Monolix Controller for nlmixr2

## Usage

``` r
monolixControl(
  nbSSDoses = 7,
  useLinearization = FALSE,
  stiff = FALSE,
  addProp = c("combined2", "combined1"),
  exploratoryAutoStop = FALSE,
  smoothingAutoStop = FALSE,
  burnInIterations = 5,
  smoothingIterations = 200,
  exploratoryIterations = 250,
  simulatedAnnealingIterations = 250,
  exploratoryInterval = 200,
  exploratoryAlpha = 0,
  omegaTau = 0.95,
  errorModelTau = 0.95,
  variability = c("none", "firstStage", "decreasing"),
  runCommand = getOption("babelmixr2.monolix", ""),
  rxControl = NULL,
  sumProd = FALSE,
  optExpression = TRUE,
  calcTables = TRUE,
  compress = TRUE,
  ci = 0.95,
  sigdigTable = NULL,
  absolutePath = FALSE,
  modelName = NULL,
  muRefCovAlg = TRUE,
  run = TRUE,
  ...
)
```

## Arguments

- nbSSDoses:

  Number of steady state doses (default 7)

- useLinearization:

  Use linearization for log likelihood and fim.

- stiff:

  boolean for using the stiff ODE solver

- addProp:

  Type of additive-plus-proportional error: \`"combined1"\`, where
  standard deviations add: \$\$y = f + (a + b\times f^c) \times
  \varepsilon\$\$; or \`"combined2"\`, where variances add: \$\$y = f +
  \sqrt{a^2 + b^2\times f^{2\times c}} \times \varepsilon\$\$. Here y =
  observed, f = predicted, a = additive sd, b = proportional/power sd, c
  = power exponent (1 in the proportional case).

- exploratoryAutoStop:

  logical to turn on or off exploratory phase auto-stop of SAEM (default
  250)

- smoothingAutoStop:

  Boolean indicating if the smoothing should automatically stop (default
  `FALSE`)

- burnInIterations:

  Number of burn in iterations

- smoothingIterations:

  Number of smoothing iterations

- exploratoryIterations:

  Number of iterations for exploratory phase (default 250)

- simulatedAnnealingIterations:

  Number of simulating annealing iterations

- exploratoryInterval:

  Minimum number of iterations in the exploratory phase (default 200)

- exploratoryAlpha:

  Convergence memory in the exploratory phase (only used when
  `exploratoryAutoStop` is `TRUE`)

- omegaTau:

  Proportional rate on variance for simulated annealing

- errorModelTau:

  Proportional rate on error model for simulated annealing

- variability:

  This describes the methodology for parameters without variability. It
  could be: - Fixed throughout (none) - Variability in the first stage
  (firstStage) - Decreasing until it reaches the fixed value
  (decreasing)

- runCommand:

  is a shell command or function to run monolix; You can specify the
  default by `options("babelmixr2.monolix"="runMonolix")`. If it is
  empty and 'lixoftConnectors' is available, use lixoftConnectors to run
  monolix. See details for function usage.

- rxControl:

  \`rxode2\` ODE solving options during fitting, created with
  \`rxControl()\`

- sumProd:

  Is a boolean indicating if the model should change multiplication to
  high precision multiplication and sums to high precision sums using
  the PreciseSums package. By default this is `FALSE`.

- optExpression:

  Optimize the rxode2 expression to speed up calculation. By default
  this is turned on.

- calcTables:

  This boolean is to determine if the foceiFit will calculate tables. By
  default this is `TRUE`

- compress:

  Should the object have compressed items

- ci:

  Confidence level for some tables. By default this is 0.95 or 95%
  confidence.

- sigdigTable:

  Significant digits in the final output table. If not specified, then
  it matches the significant digits in the \`sigdig\` optimization
  algorithm. If \`sigdig\` is NULL, use 3.

- absolutePath:

  Boolean indicating if the absolute path should be used for the monolix
  runs

- modelName:

  Model name used to generate the NONMEM output. If `NULL` try to infer
  from the model name (could be `x` if not clear). Otherwise use this
  character for outputs.

- muRefCovAlg:

  When \`TRUE\` (default), algebraic expressions that can be
  mu-referenced are internally rewritten as mu-referenced covariates and
  restored after optimization. Mirrors
  `saemControl(muRefCovAlg=)`/`nlmeControl(muRefCovAlg=)`; for
  [`foceiControl()`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.html)
  only takes effect when `muModel != "none"`.

- run:

  Should monolix be run and the results be imported to nlmixr2? (Default
  is TRUE)

- ...:

  Ignored parameters

## Value

A monolix control object

## Details

If `runCommand` is given as a string, it will be called with the
[`system()`](https://rdrr.io/r/base/system.html) command like:

`runCommand mlxtran`.

For example, if `runCommand="'/path/to/monolix/mlxbsub2021' -p "` then
the command line used would look like the following:

`'/path/to/monolix/mlxbsub2021' monolix.mlxtran`

If `runCommand` is given as a function, it will be called as
`FUN(mlxtran, directory, ui)` to run Monolix. This allows you to run
Monolix in any way that you may need, as long as you can write it in R.
babelmixr2 will wait for the function to return before proceeding.

If `runCommand` is `NA`,
[`nlmixr()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.html)
will stop after writing the model files and without starting Monolix.

Note that you can get the translated monolix components from a
parsed/compiled rxode2 ui object with `ui$monolixModel` and `ui$mlxtran`

## Author

Matthew Fidler
