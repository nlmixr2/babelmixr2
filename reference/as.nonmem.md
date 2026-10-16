# Convert an nlmixr2 model to a NONMEM control stream or Monolix project

These write the files that `nlmixr2(model, data, est="nonmem")` or
`nlmixr2(model, data, est="monolix")` would write, without running
NONMEM or Monolix. The file names are the same as the ones from
`nonmemControl(modelName=)` and `monolixControl(modelName=)`:

## Usage

``` r
as.nonmem(
  x,
  data = NULL,
  modelName = NULL,
  timeVarying = character(0),
  dose = NULL,
  doseType = "bolus",
  ss = FALSE,
  cens = FALSE,
  limit = FALSE,
  control = nonmemControl(),
  ...
)

as.monolix(
  x,
  data = NULL,
  modelName = NULL,
  timeVarying = character(0),
  dose = NULL,
  doseType = "bolus",
  ss = FALSE,
  cens = FALSE,
  limit = FALSE,
  control = monolixControl(),
  ...
)
```

## Arguments

- x:

  nlmixr2 model function, `rxode2` ui or fit

- data:

  optional data to translate with. When `NULL` and `x` is a fit (and
  none of the template arguments are given), the fit's data are used;
  otherwise a template dataset is made from the arguments below.

- modelName:

  the base name of the files. When `NULL`, the control's `modelName` is
  used, then the name of the model variable.

- timeVarying:

  character vector of the covariates that change with time (by default
  none do). Time-varying covariates are not mu-referenced (and are
  Monolix regressors).

- dose:

  the compartments receiving doses (names or numbers); by default the
  first compartment of the model (a model without compartments has no
  doses). Numbers are in the translated order, where a `linCmt()`
  model's compartments come before its ODEs

- doseType:

  the type of dose for each `dose` compartment (recycled): `"bolus"`,
  `"infusion"`, `"modelRate"` (a modeled `rate()`) or `"modelDur"` (a
  modeled `dur()`)

- ss:

  logical; does the data have steady state doses (adds `SS` and `II`)

- cens:

  logical; does the data have censored observations (adds `CENS`)

- limit:

  logical; does the data have censoring limits (adds `LIMIT`, and
  implies `cens=TRUE`)

- control:

  the
  [`nonmemControl()`](https://nlmixr2.github.io/babelmixr2/reference/nonmemControl.md)
  or
  [`monolixControl()`](https://nlmixr2.github.io/babelmixr2/reference/monolixControl.md)
  options for the translation; `runCommand` is ignored since the program
  is not run

- ...:

  other arguments (currently ignored)

## Value

invisibly, a named character vector of the written files

## Details

- NONMEM: `<modelName>-nonmem/<modelName>.nmctl` and
  `<modelName>-nonmem/<modelName>.csv`

- Monolix: `<modelName>-monolix.mlxtran`, `<modelName>-monolix.txt` and
  `<modelName>-monolix.csv`

Part of the translation depends on the data: which covariates are
time-varying (these are not mu-referenced), which compartments are dosed
and how (the Monolix administrations, `RATE`), and if the data has
steady state (`SS`/`II`) or censoring (`CENS`/`LIMIT`) columns. When
`data` is given, all of these come from the data. Otherwise a
one-subject template dataset is made from the model and the
`timeVarying`, `dose`, `doseType`, `ss`, `cens` and `limit` arguments.
This template dataset is written as the `.csv` file so you can see the
columns the model expects; replace it with your data before running
NONMEM or Monolix. With the template dataset, covariate expressions
(like `log(WT/70)`) are written in the model instead of as derived data
columns (`muRefCovAlg=FALSE`), so the files only need the model's
covariates.

When files with the same name but a different model already exist, the
new files are numbered (like `<modelName>-001-nonmem`), as in
[`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.html);
files of the same model are left alone.

## Author

Matthew L. Fidler

## Examples

``` r
# \donttest{
one.cmt <- function() {
  ini({
    tka <- 0.45
    tcl <- log(c(0, 2.7, 100))
    tv <- 3.45
    cl.wt <- 0.1
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl + cl.wt * WT)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

withr::with_tempdir({
  print(as.nonmem(one.cmt))
  print(as.monolix(one.cmt, timeVarying = "WT"))
})
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#>  
#>  
#> ℹ translating linCmt() to NONMEM's closed-form solution
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> ℹ writing nonmem files
#> ℹ done
#> ℹ only exported NONMEM control stream/data
#>                          nmctl                            csv 
#> "one.cmt-nonmem/one.cmt.nmctl"   "one.cmt-nonmem/one.cmt.csv" 
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#>  
#>  
#> ℹ translating linCmt() to Monolix's closed-form solution
#> ℹ writing monolix files
#> ℹ done
#> ℹ run Monolix manually and rerun nlmixr()
#>                   mlxtran                     model                       csv 
#> "one.cmt-monolix.mlxtran"     "one.cmt-monolix.txt"     "one.cmt-monolix.csv" 
# }
```
