# PKNCA estimation control

PKNCA estimation control

## Usage

``` r
pkncaControl(
  concu = NA_character_,
  doseu = NA_character_,
  timeu = NA_character_,
  volumeu = NA_character_,
  vpMult = 2,
  qMult = 1/2,
  vp2Mult = 4,
  q2Mult = 1/4,
  dvParam = "cp",
  groups = character(),
  sparse = FALSE,
  ncaData = NULL,
  ncaResults = NULL,
  rxControl = rxode2::rxControl()
)
```

## Arguments

- concu, doseu, timeu:

  concentration, dose, and time units from the source data (passed to
  [`PKNCA::pknca_units_table()`](https://humanpred.github.io/pknca/reference/pknca_units_table.html)).

- volumeu:

  compartment volume for the model (if `NULL`, simplified units from
  source data will be used)

- vpMult, qMult, vp2Mult, q2Mult:

  Multipliers for vc and cl to provide initial estimates for vp, q, vp2,
  and q2

- dvParam:

  The parameter name in the model that should be modified for
  concentration unit conversions. It must be assigned on a line by
  itself, separate from the residual error model line.

- groups:

  Grouping columns for NCA summaries by group (required if
  `sparse = TRUE`)

- sparse:

  Are the concentration-time data sparse PK (commonly used in small
  nonclinical species or with terminal or difficult sampling) or dense
  PK (commonly used in clinical studies or larger nonclinical species)?

- ncaData:

  Data to use for calculating NCA parameters. Typical use is when a
  subset of the original data are informative for NCA.

- ncaResults:

  Already computed NCA results (a PKNCAresults object) to bypass
  automatic calculations. At least the following parameters must be
  calculated in the NCA: tmax, cmax.dn, cl.last

- rxControl:

  Control options sent to
  [`rxode2::rxControl()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)

## Value

A list of parameters
