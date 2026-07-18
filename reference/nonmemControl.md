# NONMEM estimation control

NONMEM estimation control

## Usage

``` r
nonmemControl(
  est = c("focei", "imp", "its", "posthoc"),
  advanOde = c("advan13", "advan8", "advan6"),
  cov = c("r,s", "r", "s", ""),
  maxeval = 1e+05,
  tol = 6,
  atol = 12,
  sstol = 6,
  ssatol = 12,
  sigl = 12,
  sigdig = 3,
  print = 1,
  extension = getOption("babelmixr2.nmModelExtension", ".nmctl"),
  outputExtension = getOption("babelmixr2.nmOutputExtension", ".lst"),
  runCommand = getOption("babelmixr2.nonmem", ""),
  iniSigDig = 5,
  protectZeros = FALSE,
  muRef = TRUE,
  addProp = c("combined2", "combined1"),
  rxControl = NULL,
  sumProd = FALSE,
  optExpression = TRUE,
  calcTables = TRUE,
  compress = TRUE,
  ci = 0.95,
  sigdigTable = NULL,
  readRounding = FALSE,
  readBadOpt = FALSE,
  niter = 100L,
  isample = 1000L,
  iaccept = 0.4,
  iscaleMin = 0.1,
  iscaleMax = 10,
  df = 4,
  seed = 14456,
  mapiter = 1,
  mapinter = 0,
  noabort = TRUE,
  modelName = NULL,
  muRefCovAlg = TRUE,
  run = TRUE,
  ...
)
```

## Arguments

- est:

  NONMEM estimation method

- advanOde:

  The ODE solving method for NONMEM

- cov:

  The NONMEM covariance method

- maxeval:

  NONMEM's maxeval (for non posthoc methods)

- tol:

  NONMEM tolerance for ODE solving advan

- atol:

  NONMEM absolute tolerance for ODE solving

- sstol:

  NONMEM tolerance for steady state ODE solving

- ssatol:

  NONMEM absolute tolerance for steady state ODE solving

- sigl:

  NONMEM sigl estimation option

- sigdig:

  the significant digits for NONMEM

- print:

  The print number for NONMEM

- extension:

  NONMEM file extensions

- outputExtension:

  Extension to use for the NONMEM output listing

- runCommand:

  Command to run NONMEM (typically the path to "nmfe75") or a function.
  See the details for more information.

- iniSigDig:

  How many significant digits are printed in \$THETA and \$OMEGA when
  the estimate is zero. Also controls the zero protection numbers

- protectZeros:

  Add methods to protect divide by zero

- muRef:

  Automatically mu-reference the control stream

- addProp, sumProd, optExpression, calcTables, compress, ci,
  sigdigTable:

  Passed to
  [`nlmixr2est::foceiControl`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.html)

- rxControl:

  Options to pass to
  [`rxode2::rxControl`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)
  for simulations

- readRounding:

  Try to read NONMEM output when NONMEM terminated due to rounding
  errors

- readBadOpt:

  Try to read NONMEM output when NONMEM terminated due to an apparent
  failed optimization

- niter:

  number of iterations in NONMEM estimation methods

- isample:

  Isample argument for NONMEM ITS estimation method

- iaccept:

  Iaccept for NONMEM ITS estimation methods

- iscaleMin:

  parameter for IMP NONMEM method (ISCALE_MIN)

- iscaleMax:

  parameter for IMP NONMEM method (ISCALE_MAX)

- df:

  degrees of freedom for IMP method

- seed:

  is the seed for NONMEM methods

- mapiter:

  the number of map iterations for IMP method

- mapinter:

  is the MAPINTER parameter for the IMP method

- noabort:

  Add the `NOABORT` option for `$EST`

- modelName:

  Model name used to generate the NONMEM output. If `NULL` try to infer
  from the model name (could be `x` if not clear). Otherwise use this
  character for outputs.

- muRefCovAlg:

  This controls if algebraic expressions that can be mu-referenced are
  treated as mu-referenced covariates by:

  1\. Creating a internal data-variable \`nlmixrMuDerCov#\` for each
  algebraic mu-referenced expression

  2\. Change the algebraic expression to \`nlmixrMuDerCov# \*
  mu_cov_theta\`

  3\. Use the internal mu-referenced covariate for saem

  4\. After optimization is completed, replace \`model()\` with old
  \`model()\` expression

  5\. Remove \`nlmixrMuDerCov#\` from nlmix2 output

  In general, these covariates should be more accurate since it changes
  the system to a linear compartment model. Therefore, by default this
  is \`TRUE\`.

- run:

  Should NONMEM be run (and the files imported to nlmixr2); default is
  TRUE, but FALSE will simply create the NONMEM control stream and data
  file.

- ...:

  optional `genRxControl` argument controlling automatic `rxControl`
  generation.

## Value

babelmixr2 control option for generating NONMEM control stream and
reading it back into `babelmixr2`/`nlmixr2`

## Details

If `runCommand` is given as a string, it will be called with the
[`system()`](https://rdrr.io/r/base/system.html) command like:

`runCommand controlFile outputFile`.

For example, if `runCommand="'/path/to/nmfe75'"` then the command line
used would look like the following:

`'/path/to/nmfe75' one.cmt.nmctl one.cmt.lst`

If `runCommand` is given as a function, it will be called as
`FUN(ctl, directory, ui)` to run NONMEM. This allows you to run NONMEM
in any way that you may need, as long as you can write it in R.
babelmixr2 will wait for the function to return before proceeding.

If `runCommand` is `NA`,
[`nlmixr()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.html)
will stop after writing the model files and without starting NONMEM.

## Author

Matthew L. Fidler

## Examples

``` r

nonmemControl()
#> $est
#> [1] "focei"
#> 
#> $cov
#> [1] "r,s"
#> 
#> $advanOde
#> [1] "advan13"
#> 
#> $maxeval
#> [1] 1e+05
#> 
#> $print
#> [1] 1
#> 
#> $noabort
#> [1] TRUE
#> 
#> $iniSigDig
#> [1] 5
#> 
#> $tol
#> [1] 6
#> 
#> $atol
#> [1] 12
#> 
#> $sstol
#> [1] 6
#> 
#> $ssatol
#> [1] 12
#> 
#> $sigl
#> [1] 12
#> 
#> $muRef
#> [1] TRUE
#> 
#> $sigdig
#> [1] 3
#> 
#> $protectZeros
#> [1] FALSE
#> 
#> $runCommand
#> [1] ""
#> 
#> $outputExtension
#> [1] ".lst"
#> 
#> $addProp
#> [1] "combined2"
#> 
#> $rxControl
#> $scale
#> NULL
#> 
#> $method
#> liblsoda 
#>        2 
#> 
#> $atol
#> [1] 1e-12
#> 
#> $rtol
#> [1] 1e-06
#> 
#> $maxsteps
#> [1] 70000
#> 
#> $hmin
#> [1] 0
#> 
#> $hmax
#> [1] NA
#> 
#> $hini
#> [1] 0
#> 
#> $maxordn
#> [1] 12
#> 
#> $maxords
#> [1] 5
#> 
#> $covsInterpolation
#> nocb 
#>    2 
#> 
#> $addCov
#> [1] TRUE
#> 
#> $returnType
#> rxSolve 
#>       0 
#> 
#> $sigma
#> NULL
#> 
#> $sigmaDf
#> NULL
#> 
#> $nCoresRV
#> [1] 1
#> 
#> $sigmaIsChol
#> [1] FALSE
#> 
#> $sigmaSeparation
#> [1] "auto"
#> 
#> $sigmaXform
#> identity 
#>        4 
#> 
#> $nDisplayProgress
#> [1] 10000
#> 
#> $amountUnits
#> [1] NA
#> 
#> $timeUnits
#> [1] "hours"
#> 
#> $addDosing
#> [1] FALSE
#> 
#> $stateTrim
#> [1] Inf
#> 
#> $updateObject
#> [1] FALSE
#> 
#> $omega
#> NULL
#> 
#> $omegaDf
#> NULL
#> 
#> $omegaIsChol
#> [1] FALSE
#> 
#> $omegaSeparation
#> [1] "auto"
#> 
#> $omegaXform
#> variance 
#>        6 
#> 
#> $nSub
#> [1] 1
#> 
#> $thetaMat
#> NULL
#> 
#> $thetaDf
#> NULL
#> 
#> $thetaIsChol
#> [1] FALSE
#> 
#> $nStud
#> [1] 1
#> 
#> $dfSub
#> [1] 0
#> 
#> $dfObs
#> [1] 0
#> 
#> $seed
#> NULL
#> 
#> $nsim
#> NULL
#> 
#> $minSS
#> [1] 10
#> 
#> $maxSS
#> [1] 10000
#> 
#> $strictSS
#> [1] 1
#> 
#> $infSSstep
#> [1] 12
#> 
#> $istateReset
#> [1] TRUE
#> 
#> $subsetNonmem
#> [1] TRUE
#> 
#> $hmaxSd
#> [1] 0
#> 
#> $maxAtolRtolFactor
#> [1] 0.1
#> 
#> $from
#> NULL
#> 
#> $to
#> NULL
#> 
#> $by
#> NULL
#> 
#> $length.out
#> NULL
#> 
#> $iCov
#> NULL
#> 
#> $keep
#> NULL
#> 
#> $keepF
#> character(0)
#> 
#> $drop
#> NULL
#> 
#> $warnDrop
#> [1] TRUE
#> 
#> $omegaLower
#> [1] -Inf
#> 
#> $omegaUpper
#> [1] Inf
#> 
#> $sigmaLower
#> [1] -Inf
#> 
#> $sigmaUpper
#> [1] Inf
#> 
#> $thetaLower
#> [1] -Inf
#> 
#> $thetaUpper
#> [1] Inf
#> 
#> $indLinPhiM
#> [1] 0
#> 
#> $indLinPhiTol
#> [1] 1e-07
#> 
#> $indLinMatExpType
#> expokit 
#>       2 
#> 
#> $indLinMatExpOrder
#> [1] 6
#> 
#> $idFactor
#> [1] TRUE
#> 
#> $mxhnil
#> [1] 0
#> 
#> $hmxi
#> [1] 0
#> 
#> $warnIdSort
#> [1] TRUE
#> 
#> $ssAtol
#> [1] 1e-12
#> 
#> $ssRtol
#> [1] 1e-06
#> 
#> $safeZero
#> [1] 0
#> 
#> $sumType
#> pairwise 
#>        1 
#> 
#> $prodType
#> long double 
#>           1 
#> 
#> $resample
#> NULL
#> 
#> $resampleID
#> [1] TRUE
#> 
#> $maxwhile
#> [1] 100000
#> 
#> $cores
#> [1] 0
#> 
#> $atolSens
#> [1] 1e-08
#> 
#> $rtolSens
#> [1] 1e-06
#> 
#> $ssAtolSens
#> [1] 1e-08
#> 
#> $ssRtolSens
#> [1] 1e-06
#> 
#> $simVariability
#> [1] NA
#> 
#> $nLlikAlloc
#> NULL
#> 
#> $useStdPow
#> [1] 0
#> 
#> $naTimeHandle
#> ignore 
#>      1 
#> 
#> $addlKeepsCov
#> [1] FALSE
#> 
#> $addlDropSs
#> [1] TRUE
#> 
#> $ssAtDoseTime
#> [1] TRUE
#> 
#> $ss2cancelAllPending
#> [1] FALSE
#> 
#> $naInterpolation
#> locf 
#>    1 
#> 
#> $keepInterpolation
#> na 
#>  2 
#> 
#> $safeLog
#> [1] 1
#> 
#> $safePow
#> [1] 1
#> 
#> $ssSolved
#> [1] TRUE
#> 
#> $linCmtSensType
#> auto 
#>  100 
#> 
#> $linCmtSensH
#> [1] 1e-04
#> 
#> $linCmtGillFtol
#> [1] 0
#> 
#> $linCmtGillK
#> [1] 20
#> 
#> $linCmtGillStep
#> [1] 4
#> 
#> $linCmtGillRtol
#> [1] 1.490116e-08
#> 
#> $linCmtShiErr
#> [1] 1.490116e-08
#> 
#> $linCmtShiMax
#> [1] 20
#> 
#> $linCmtScale
#> [1] 0 0 0 0 0 0 0
#> 
#> $linCmtHcmt
#> [1] 1
#> 
#> $linCmtHmeanI
#> geometric 
#>         2 
#> 
#> $linCmtHmeanO
#> geometric 
#>         2 
#> 
#> $linCmtSuspect
#> [1] 1e-06
#> 
#> $linCmtForwardMax
#> [1] 2
#> 
#> $indOwnAlloc
#> [1] -1
#> 
#> $maxExtra
#> [1] 1000
#> 
#> $tolFactor
#> NULL
#> 
#> $serializeFile
#> NULL
#> 
#> $dense
#> [1] FALSE
#> 
#> $cvodeLinSolver
#> dense 
#>     1 
#> 
#> $stiff2
#> [1] 0
#> 
#> $autoSwitchMaxStiff
#> [1] 10
#> 
#> $autoSwitchMaxNonstiff
#> [1] 3
#> 
#> $autoSwitchStiffFirst
#> [1] 0
#> 
#> $autoSwitchNonstifftol
#> [1] 0.9
#> 
#> $autoSwitchStifftol
#> [1] 0.9
#> 
#> $autoSwitchDtfac
#> [1] 2
#> 
#> $autoSwitchSwitchMax
#> [1] 5
#> 
#> $useLinCmt
#> [1] TRUE
#> 
#> $file
#> NULL
#> 
#> $chunkSize
#> NULL
#> 
#> $parallel
#> [1] 0
#> 
#> $.zeros
#> NULL
#> 
#> attr(,"class")
#> [1] "rxControl"
#> 
#> $sumProd
#> [1] FALSE
#> 
#> $optExpression
#> [1] TRUE
#> 
#> $calcTables
#> [1] TRUE
#> 
#> $compress
#> [1] TRUE
#> 
#> $ci
#> [1] 0.95
#> 
#> $sigdigTable
#> NULL
#> 
#> $readRounding
#> [1] FALSE
#> 
#> $readBadOpt
#> [1] FALSE
#> 
#> $genRxControl
#> [1] TRUE
#> 
#> $niter
#> [1] 100
#> 
#> $isample
#> [1] 1000
#> 
#> $iaccept
#> [1] 0.4
#> 
#> $iscaleMin
#> [1] 0.1
#> 
#> $iscaleMax
#> [1] 10
#> 
#> $df
#> [1] 4
#> 
#> $seed
#> [1] 14456
#> 
#> $mapiter
#> [1] 1
#> 
#> $modelName
#> NULL
#> 
#> $muRefCovAlg
#> [1] TRUE
#> 
#> $run
#> [1] TRUE
#> 
#> attr(,"class")
#> [1] "nonmemControl"
```
