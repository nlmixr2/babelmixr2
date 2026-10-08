# babelmixr2 0.1.11.9000

* Monolix imports (`as.nlmixr2()` of a monolix2rx model) no longer fail
  when the project's steady-state dose count (`nbdoses`) is below 6:
  `minSS`/`maxSS` come from monolix2rx's `.getSsLimits()`, raised to
  rxode2's floor of 5 and 7.

* `simplifyUnit()` and `modelUnitConversion()` work without the
  suggested packages testthat and units installed: `simplifyUnit()`
  checked its arguments with checkmate's `expect_*()` (testthat
  expectations) instead of `assert_*()`, and both now ask for units
  with `rxode2::rxReq()`.  Their examples run only when units is
  installed, and the PKNCA and saemix tests skip when those packages
  are not installed (CRAN's noSuggests check).

* Monolix (fourth stress kit run): `EQUATION:` only has assignments
  and conditions, so PK macros there were rejected.  The model lines
  before the ODEs (with the compartment properties the macros use) are
  now written in the `PK:` block before the macros, and the ODEs stay
  in `EQUATION:`.  Compartment properties set after the ODEs are moved
  before them, and a property set in an `if` is set in every branch
  (rxode2's default in the others) instead of from a default assigned
  before the `if`, which Monolix cannot reassign.  A missing `LIMIT`
  is written as missing instead of `-Inf`.

* Monolix projects Monolix could not load (third stress kit run, with
  Monolix's reasons):

  - "Conflicting variable definition": Monolix assigns a variable once,
    so a reassigned variable (like `if (SEX == 1) cl <- cl * 1.2`, or
    the pruned `else if`) gets a new name (`cl_rx1`); an `if`/`else`
    that changes a defined variable is pruned first.
  - "Undefined variable": compartment properties set in the model
    (`f()`, `alag()`, `dur()`, `rate()`) are defined in `EQUATION:`, so
    the PK macros using them are now written in `EQUATION:` after them
    (not in a `PK:` block before).  A modeled rate uses
    `Tk0=rx_tk0_<state>` with `rx_tk0_<state> = amtDose/rate` (macro
    arguments cannot be calculations).
  - `monolixControl(stiff=TRUE)` wrote `odeType = stiff` where Monolix
    expects an input definition; it is now in `EQUATION:`.
  - A `LIMIT` on uncensored observations (M2) is refused: Monolix only
    uses `LIMIT` with censored values.

* Monolix fits use the number of observations Monolix used for the
  objective function, like NONMEM fits; it was off for data with
  `evid=2` records or time not starting at zero.

* Monolix: a mu-referenced parameter that stays inside `exp()` in the
  model (like `cl <- exp(tcl + eta.cl) * (CRCL/100)^cl.crcl`, or
  `f(depot) <- exp(lfdepot)`) was given a log-normal distribution while
  the model still took `exp()` of it, so `exp()` was applied twice.
  It is now a normal parameter on nlmixr2's scale (its initial value,
  estimate and covariance too).
* When lixoftConnectors cannot load or run a Monolix project, the error
  now includes Monolix's reason (its `[ERROR]` lines) instead of
  "see Monolix's [ERROR] above".

* Fixes found by running the NONMEM/Monolix stress kit with NONMEM 7.4:

  - `nonmemControl(est="imp")` and `est="its"` fits were always
    reported as unsuccessful; NONMEM ends them with "OPTIMIZATION WAS
    NOT TESTED FOR CONVERGENCE", which is now read as finished.
  - `nonmemControl(est="posthoc")` fits could not be read (no `#TERM:`
    block and no parameter history).
  - Box-Cox, Yeo-Johnson and logit + Yeo-Johnson residual errors were
    rejected by NM-TRAN (IPRED redefined in nested/`ELSE IF`
    structures); the transformations are now written without nested
    IFs.  They are still transform-both-sides (DV in `CCONTR`, IPRED in
    `$ERROR`).  Predictions at or below `sqrt(DBL_EPSILON)` are floored
    like rxode2 instead of being set to `-1000000000` (which gave
    `lnorm()` fits a huge objective function), and Box-Cox no longer
    divides by zero when lambda is 0.
  - The `IPRED`/`PRED` comparison with NONMEM now undoes the residual
    transformation first (NONMEM's predictions are on the transformed
    scale).
  - Censored (`CENS`, `CENS` + `LIMIT`, `LIMIT` only) models were
    rejected by NM-TRAN ("random variable is defined in a nested IF
    structure").
  - A `$PK` parameter changed later in the model (like
    `if (SEX == 1) cl <- cl * 1.2`) was changed in `$DES`, which NONMEM
    does not allow, and used undefined in `$ERROR`; it is now defined
    as `RXPK_<name>` in `$PK`.
  - Models without random effects are refused (NONMEM treats them as
    single-subject data, where `METHOD=COND` is invalid).
* Monolix fixes found by running the NONMEM/Monolix stress kit:

  - Infusions given by rate (`RATE`) or duration (`TINF`) were written
    to the data but not declared in `[CONTENT]`, so Monolix gave them
    as bolus doses.
  - Endpoints with dotted names (like `concentration.in.plasma`) kept
    the dots in the Monolix project (`OUTPUT`, observation names,
    predictions), which Monolix cannot load.
  - Reading a fit with a fixed parameter (like `tfdepot <- fix(0.8)`)
    failed with "subscript out of bounds"; Monolix's covariance has
    only the estimated parameters.
  - A `monolixControl(runCommand=)` command that finished without
    writing Monolix's output now stops with an error instead of
    waiting for it forever.

* `nonmemControl(cov="")` now skips the covariance step as documented;
  before it gave an error because `match.arg()` cannot match `""`.

* The NONMEM/Monolix stress test (`inst/stress`) is now a kit to run on
  a machine with NONMEM and/or Monolix: `install-kit.R` installs the
  versions to test, `run-stress.R --check` finds NONMEM/Monolix, and
  `run-stress.R --kit` runs every case end to end, checks that rxode2
  reproduces the NONMEM/Monolix predictions and that a second fit reads
  the saved output, and zips the results to send back.  New cases cover
  missing observations, `evid=2`, character IDs, time not starting at
  zero, extra data columns, a single subject, oral plus iv doses, ODE
  infusions (rate, steady state, modeled duration), two endpoints,
  bounded thetas, several covariates, models without random effects
  and the NONMEM/Monolix estimation options.

* NONMEM control streams now use `$ABBR PROTECT` so NM-TRAN replaces
  `LOG`, `EXP`, `SQRT`, division and powers with NONMEM's protected
  functions (NONMEM 7.4 or later).  Turn it off with
  `nonmemControl(protect=FALSE)` or `options(babelmixr2.nmProtect=FALSE)`.
  Since NM-TRAN writes `B**E` as `PEXP(E*PLOG(B))`, integer powers are
  written as products (or powers of `x*x`, which is never negative) and
  powers of a positive number with `DEXP()`.  Functions of numbers
  (like `exp(0)`, `log(2*pi)` or `expit(0)`) and divisions of numbers
  are written as numbers, and model
  variables named like a protected function (like `plog`) are renamed.
  babelmixr2's own zero protection (`protectZeros`) is now only used
  with `protect=FALSE`, like for NONMEM before 7.4 (#62).

* NONMEM control streams now explain the zero-protection code
  babelmixr2 adds: each `RXDZ###` `IF` block is preceded by a comment
  saying what it keeps the variable away from and why, as is the
  `IF (W1 .EQ. 0.0)` residual variance protection (#91).
  The protection for `lfactorial()`/`lgamma1p()` now keeps its argument
  above `-1` (it previously clamped it just below `-1`, making `x+1`
  negative), `log(x)` and `1/x` no longer share one protected
  variable, so `1/x` keeps the sign of a negative `x`, and a variable
  reassigned in the model is protected again instead of reusing the
  protection of its old value.  Zero protection needed by `f()`,
  `alag()`, `rate()` or `dur()` (written in `$PK`) is no longer shared
  with other lines, which could use it before (or without) it being
  calculated.
* `est="pknca"` now works with covariates in the data and with a mix of
  intravascular and extravascular doses (#102).  Doses into a
  compartment that the observations are calculated from are
  intravascular.  With both routes, `ka` is estimated from the
  extravascular doses and `vc` and `cl` from the intravascular doses.
  Intravascular bolus doses have the concentration at the time of
  dosing back-extrapolated (replacing a predose concentration at the
  first dose), and other doses have it imputed (as the predose
  concentration, or zero for the first dose).  When no doses are only
  extravascular, `ka` is not updated.  Multiple-dose data no longer need
  a concentration at each dose time; each dose until the next (with at
  least 2 concentrations) is used, with `vc` from the first dose of each
  route and `cl` from dosing intervals mostly covered by concentrations.
  Doses at the same time are combined.
* When a NONMEM run fails, `est="nonmem"` now says why and where to look
  instead of failing with an unclear error: a run command that was not
  found or wrote no output, a NONMEM license problem, an NM-TRAN error in
  the control stream or the data, NONMEM not starting (for example a
  compiler problem), NONMEM crashing or stopping during estimation, and
  output babelmixr2 cannot read (#46).
* `est="pknca"` now updates the initial estimates of models that are not
  mu-referenced, like `ka <- tka * exp(eta.ka)`, instead of failing with
  "Must have names" (#101).  When there is no `vc`, a central volume
  named `v`, `V`, `Vc`, `VC`, `v1` or `V1` (only one of them) now gets
  the NCA central volume estimate.  Parameters defined as `expit()` of a
  theta are now transformed back correctly.  A message lists the
  parameters that could not be updated.
* NONMEM models can now use nested `if`/`else if`/`else` statements
  (and `ifelse()`): their branches are pruned with rxode2's branch
  pruning before the model is translated to NONMEM (#11).  With the
  default `nonmemControl(prune="auto")` a model whose `if` blocks are
  simple is still written with NONMEM `IF` blocks and only a model that
  needs it is pruned; `prune=TRUE` always prunes and `prune=FALSE` never
  prunes (the error for unsupported `if`/`else` statements suggests
  `prune=TRUE`).

* `monolixControl(prune=)` has the same option for Monolix.  Monolix
  writes `if`/`elseif`/`else` (and nested `if`) statements directly, so
  `prune="auto"` (the default) only prunes a model that uses `ifelse()`,
  which Monolix cannot write; `prune=TRUE` always prunes and
  `prune=FALSE` never prunes.  A logical expression used as a number in
  a Monolix model is written as a 0/1 indicator variable (#11).

* A logical expression used as a number in a NONMEM model (like
  `cl <- tcl * (WT > 70)`) is now written as a 0/1 indicator variable,
  since NONMEM cannot use a logical expression as a number.  A numeric
  `if ()` condition is written as not equal to zero.
* `est="nonmem"` and `est="monolix"` fits with
  `table=tableControl(cwres=TRUE)` no longer fail with "objective
  function 'FOCEi' already present".  The fit keeps both the NONMEM (or
  Monolix) objective, which stays in use, and nlmixr2's FOCEi objective
  (#94).  `as.nlmixr2()` of a `nonmem2rx` or `monolix2rx` model now also
  keeps the imported objective in use with `cwres=TRUE`, and no longer
  fails after an earlier import with `cwres=TRUE` (it started from that
  fit's etas).

* `est="monolix"` now translates compartment properties (`f()`,
  `alag()`, `rate()` and `dur()`) that are expressions, like
  `f(depot) <- exp(lfdepot)`, instead of erroring with "the complex F is
  not supported by babelmixr2" (#115).  The expression is calculated in
  a new variable (like `rx_f_depot`) that the `PK:` macro uses.  A
  property set only inside an `if` now keeps rxode2's default otherwise
  (1 for `f()`, 0 for `alag()`) instead of being applied unconditionally.
* The ACoP 2024 `babelmixr2`/`PopED` abstract
  ([doi:10.70534/XUMG6226](https://doi.org/10.70534/XUMG6226)) is now
  in `citation("babelmixr2")`, the package description and the `PopED`
  article (#155).
* `as.nlmixr2()` of a `nonmem2rx` model now stops with an informative
  error when the model still contains untranslated NONMEM residual
  variables (`eps#` or `err#`) instead of failing inside the estimation
  routine (#95).
* `popedControl(sigdig=)` now sets the ODE solver tolerances itself
  (`atol = rtol = 0.5*10^(-sigdig-2)`, the same values as before) instead
  of using `rxode2::rxControl(sigdig=)`.  rxode2 5.1.5 loosened the
  tolerances that `rxControl(sigdig=)` gives, which made the
  finite-difference FIM, and so the design OFV and RSEs, less accurate
  (#223).
* `est="nonmem"` now fits censored data the way nlmixr2 does (#92).  M3
  (`CENS`), M4 (`CENS` with a finite `LIMIT`) and M2 (`CENS=0` with a
  finite `LIMIT`, including data with a `LIMIT` but no `CENS` column)
  use `F_FLAG` likelihoods in `$ERROR` with `LAPLACIAN` estimation.  A
  missing `LIMIT` is now written as NONMEM's infinity instead of `0`,
  `CENS`/`LIMIT` columns that do not censor anything are dropped, the
  objective function is adjusted so the log-likelihood includes the
  censored observations correctly, and censoring with a transformed
  endpoint (like `lnorm()`) is refused instead of giving the wrong
  likelihood.  The censored observations are left out of the NONMEM
  `PRED` comparison since NONMEM's `PRED` is their likelihood.

* `est="nonmem"` and `est="monolix"` now fit `linCmt()` models.  A pure
  `linCmt()` model uses NONMEM's closed-form solutions (`ADVAN1`-`ADVAN4`,
  `ADVAN11` or `ADVAN12` with `TRANS1` micro-constants) or Monolix's
  `pkmodel()`; a model the closed form cannot represent (other ODEs,
  parameters that change with time, and for Monolix modeled rates or
  durations, doses into more than one compartment or amounts used in
  the model) is translated to ODEs with `rxode2::linToOde()`.
  `nonmemControl(linCmt="ode")` and `monolixControl(linCmt="ode")`
  always use ODEs.  The closed-form solutions need an rxode2 with
  `linCmtMicro()`; with an older rxode2, `linCmt()` models use ODEs.

* `est="poped"` now designs `linCmt()` models (translated to ODEs).
  `est="nlmer"`, `est="fmeMcmc"` and `est="pseudoOptim"` are now tested
  with `linCmt()` models.

* New NONMEM/Monolix stress test (`inst/stress/`), shared by the package
  tests and a command line runner (`run-stress.R`).  The runner can also
  fit every case with NONMEM and/or Monolix end to end on a machine that
  has them, and optionally translate the nlmixr2lib models; see
  `inst/stress/README.md`.

* The NONMEM data dropped the `RATE`, `SS` and `II` items: which items
  were written was decided from settings made only after the data were
  converted, so infusions, steady state doses and modeled rates or
  durations were written as plain bolus doses.  The items are now kept
  when the data use them.

* A steady state dose with a lag time (which rxode2 splits in two) is now
  written as one `SS` dose for NONMEM and Monolix, not as two bolus
  doses.

* NONMEM now writes a modeled duration as `Dn` (it was `DURn`, which
  NONMEM does not know) and a modeled rate as `Rn` (it was dropped).

* NONMEM now uses the NONMEM name of a mu-referenced covariate in `$PK`
  (like `NLMIXRMUDERCOV1`), matching `$INPUT`.

* `probitInv()` (which rxode2 writes with `erf()`) now translates to
  NONMEM (`PHI()`) and Monolix (`normcdf()`).

* `est="monolix"` now refuses models without between-subject
  variability, residual errors other than `add()`, `prop()` and
  `add() + prop()`, and transformations other than `lnorm()` and
  `logitNorm()` up front with a clear error.  Before, a `boxCox()`,
  `yeoJohnson()` or `logitNorm()` residual stopped with "argument must be
  a character string".

* Monolix now separates the arguments of the `empty()` macro and the
  `logitNormal` distribution with commas.

* rxode2's normalized powers (`Rx_pow_di()`, `Rx_pow()`) now translate to
  NONMEM (`**`) and Monolix (`^`).

* NONMEM and Monolix data without doses (like a `$PRED` model with the
  dose as a covariate) no longer stop with "undefined columns selected".

* `est="nonmem"` now refuses residual transformations it cannot write
  (like `probitNorm()`) up front; before, writing the files stopped with
  "can only write character objects".

* Monolix now applies a bioavailability or lag time to the doses of that
  compartment only; a property of a compartment without doses no longer
  stops with "values must be length 1".  A parameter written with
  brackets (like `v <- (tv + eta.v)`) is now a normal parameter, and an
  unsupported parameter transformation (like `log10()`) gives a clear
  error.

* `est="monolix"` now spells a mu-referenced parameter with a `.` in its
  name (like the tainted `rx__cl.wt`, or `ka.x <- exp(tka + eta.ka)`) the
  same way in every section of the Monolix model (#220).  The `EQUATION:`
  block already wrote it as `rx__cl__wt`, but `input=`, `[INDIVIDUAL]`,
  `<PARAMETER>` and the output readers used `rx__cl.wt`, so Monolix saw an
  undeclared variable.

* `est="saemix"` now fits `linCmt()` models.  The prediction was looked
  up in a column named after the endpoint (`rxLinCmt`), which the solved
  model does not output, so saemix stopped with `non-numeric argument to
  function` (#212).

* `est="saemix"` now fits models where a structural theta has no
  between-subject variability (e.g. `v <- exp(tv)`).  Collecting the
  individual etas after the fit failed with `invalid subscript type
  'list'` (#212).

* `est="saemix"` now refuses a model it cannot fit, instead of fitting it
  with a different residual error.  saemix fits one endpoint with an
  `add()`, `prop()`, `add() + prop()` (`combined2`, the only combination
  saemix has) or `lnorm()` residual error, or an `ll()` likelihood.  A
  model with more than one endpoint (previously fit against predictions of
  zero), a `combined1` `add() + prop()` (including
  `saemixControl(addProp="combined1")`), `pow()`, `boxCox()`,
  `yeoJohnson()`, a logit/probit transformation, `lnorm() + prop()` or a
  non-normal residual distribution now stops with an error, as does a
  fixed residual error or between-subject variability, which saemix
  would otherwise estimate anyway.  The checks use the new rxode2
  assertions `assertRxUiTransform()`, `assertRxUiErrType()`,
  `assertRxUiAddProp()`, `assertRxUiNoFixedResiduals()` and
  `assertRxUiNoFixedOmega()`, so this requires rxode2 5.1.8 (#212).

* `est="saemix"` now fits `lnorm()` residual errors with saemix's
  exponential error model; they were previously fit as an additive error
  with a missing starting value (#212).
* `est="monolix"` now accepts normal priors from `ini({})` and writes them
  as Monolix MAP estimation (#207).  A parameter with a prior is estimated
  with `method=MAP`, and its prior is written to a `[POPULATION]` section
  of `<MODEL>`.  Monolix's prior on a typical value has the same
  distribution as the parameter itself, with its `sd` in the Gaussian
  space, so the prior mean is back-transformed like the estimate
  (`exp()`, `expit()`, `probitInv()`) while the prior sd is written as is:
  `prior(tka) ~ dnorm(log(1.5), 0.5)` becomes
  `ka_pop = {distribution=logNormal, typical=1.5, sd=0.5}`, the same
  distribution with no approximation.  Covariate effects get a `normal`
  prior.  Checked with Monolix 2024R1: tight priors pin `ka_pop`,
  `cl_pop`, a covariate effect and a logit-normal parameter at the prior
  mean, and a vague prior leaves the estimate at the MLE.  Priors Monolix
  cannot honour are errors rather than being dropped: priors on omega
  elements or omega blocks, multivariate normal priors, non-normal
  priors, priors on a `probitInv()` parameter with bounds other than
  (0, 1), and priors on residual error parameters -- Monolix accepts a
  MAP prior on `add__sd` but ignores it (every estimate identical to the
  run without it).

* The "PRED absolute difference compared to Monolix PRED" line of a Monolix
  fit's message is now an absolute difference (it printed a relative one).
  The covariance of a fit from Monolix 2020 or later now carries nlmixr2's
  parameter names (`tka`, `cl.wt`) instead of Monolix's (`ka_pop`,
  `beta_cl_lWT`), like fits from older Monolix versions already did.

* Monolix projects with mu-referenced covariates (`cl <- exp(tcl + eta.cl +
  cl.wt * lWT)`) now load in Monolix.  The covariate was missing from the
  `[INDIVIDUAL]` inputs (Monolix: `Undefined variable 'lWT'`), and when
  every covariate was mu-referenced the structural model got a regressor
  line with no name (`= {use=regressor}`, a syntax error).  Reading the
  results of such a fit back failed with `subscript out of bounds`: the
  covariance looked up the covariate effect as `NA_pop` instead of
  `beta_cl_lWT`.  The tests now replay Monolix 2024R1 runs, with and
  without MAP priors.  When
  lixoftConnectors cannot load or run the project, `nlmixr2()` now stops
  with an error instead of waiting forever for output Monolix never
  writes.

* `est="nonmem"` now runs models with `ini({})` priors, translating them
  to NONMEM's `$PRIOR NWPRI` (#205).  Normal priors on population
  parameters (`dnorm()`, `stdNormal()`, the `tcl + tv ~ c(...)` joint
  normal) become `$THETAP`/`$THETAPV`, and `invWishart(nu)` degrees of
  freedom on an omega block become `$OMEGAP`/`$OMEGAPD`, with the block's
  own initial estimate as the prior scale.  NWPRI gives its priors to the
  first THETAs and the first omega blocks, so the parameters with a prior
  have to come first in `ini({})`; otherwise, and for priors NWPRI cannot
  express (`dcauchy()`, a normal prior directly on an omega element, which
  is TNPRI), the model is refused before any file is written instead of
  fitting a different prior.  When the output is read back, the prior
  values NM-TRAN adds as extra THETAs and OMEGAs are dropped, and the
  objective function type says `nwpri` because NONMEM's objective
  includes the prior.

* `$OMEGA BLOCK()` records of 3 or more etas are now written in the order
  NONMEM reads them (row by row down the lower triangle).  They used to be
  written column by column, so NONMEM started from the wrong initial
  omega values.

* The `$PROBLEM` record of a generated NONMEM control stream now carries
  the model name (`$PROBLEM one.cmt translated from babelmixr2`).  It read
  a misspelled getter and was always blank (#209).  Because the control
  stream changes, an existing NONMEM export that has a `.md5` hash file
  will not match and is re-run once in a new numbered directory.
  Moving past a second stale export (`-001-nonmem` also not matching) no
  longer hangs: the export directory kept its cached number and the
  hash check looped forever.
* `est="fmeMcmc"` now uses priors declared in the model's `ini({})` block
  (for example `prior(tka) ~ dnorm(0, 10)`) instead of refusing the model.
  They become the `prior` function `FME::modMCMC()` samples with,
  evaluated with rxode2's shared prior kernel on the natural parameter
  scale, even when `scaleType` makes FME sample a rescaled space.
  Supplying `fmeMcmcControl(prior=)` as well is an error rather than
  silently preferring one of them (#208).

* `nonmemControl(est="its")` now writes `$ESTIMATION METHOD=ITS
  INTERACTION` (iterative two stage).  It wrote `METHOD=IMP`, so NONMEM
  ran importance sampling while the returned fit was labelled with the
  `nonmem its` objective function type (#211).

* A PopED design dataset that gives `cmt` as a compartment *number*
  (`et(amt=180, cmt=1)`) now doses the right compartment.  `et()` keeps
  `cmt` as a character column, so `rxode2::etTrans()` read `"1"` as a
  compartment *name*, found no match and quietly moved the dose to an
  extra compartment; the design built without a warning but every
  prediction was zero and the FIM was degenerate (#201).  This also works
  when the column mixes names and numbers, which is what a multiple
  endpoint design looks like when it names the endpoint on its
  observation records.  A dosing record that still cannot be matched to a
  model compartment is now an error instead of a silently empty design.

* A multiple endpoint PopED design can now name its endpoints with `cmt`
  (`cmt="cp"`, `cmt="eff"`) instead of `dvid`.  The `cmt` fallback was
  already written but unreachable: a dataset without a `dvid` column
  stopped with `attempt to select less than one element in get1index`
  before it was tried.  This applies to the usual design space; a design
  that gives per-`ID` sampling through `popedControl(a=)` still needs
  `dvid`.

* The PopED model translation no longer drops the `if ()` condition that
  guards an adaptive dosing call (`evid_()`, `bolus()`, `infuse()`,
  `infuseDur()`, `reset()`, ...).  The branch pruner used to flatten the
  model unconditionally, so a model like `if (t <= 0) infuseDur(DOSE,
  TINF, cmt=1)` pushed a dose at *every* design point instead of once
  (#131).  The pruner's capture protocol is now used and the guarded call
  is restored after the branches are flattened.

* Added two PopED examples showing how to make the dosing regimen itself
  optimizable (#131):

  - `inst/poped/ex.10.PKPD.HCV.dose-and-tinf.babelmixr2.R` keeps the dose
    record and makes the amount and infusion duration design (`a`)
    variables via `f(depot) <- DOSE` (with `amt=1`) and
    `dur(depot) <- TINF` (with `rate=-2`).

  - `inst/poped/ex.10.PKPD.HCV.adaptive-dosing.babelmixr2.R` drops the
    dose records entirely and pushes the regimen from inside the model
    with `infuseDur()`, which makes the dosing *interval* a design
    variable as well.  This one needs rxode2 > 5.1.7 (rxode2#1214).

  Both are optimized with `poped_optim(..., opt_a=TRUE)` and agree on the
  reference design (OFV 88.27).

* The NONMEM/Monolix fit cache is now written with `saveRDS()` as
  `<model>.rds` / `nlmixr.rds` instead of `qs2`, so `qs2` moved from
  `Imports` to `Suggests`.  Existing run directories keep working: a
  `.qs2` cache is read once (when `qs2` is installed) and rewritten as
  the `.rds`, and if it cannot be read the fit is rebuilt from the run
  output as it would be for any missing cache.
* Each estimation method now carries `type` and `description` attributes so it
  appears in the category-grouped method list nlmixr2est prints for an
  unsupported `est=` (or a bare `nlmixr2()` call): `nonmem`, `monolix`, `pknca`,
  `fmeMcmc` and `pseudoOptim` under "External", `saemix` under "Stochastic EM",
  `nlmer` under "Integral approximation", and `poped` under "Optimal Design".

* The mu-referenced covariate algorithm (`muRefCovAlg`) is now applied
  through the `nlmixr2est` preprocessing/post-final-object hooks instead
  of explicit `nlmixr2est::.uiApplyMu2()`/`.uiFinalizeMu2()` calls in the
  `saemix`, `nonmem`, `monolix`, and `nlmer` estimation methods (#184).
  The `nonmem` and `monolix` methods gained the `mu` method attribute so
  the hooks fire for them.

* The `nlmer` estimation method now prints its iterations during the
  `lme4::nlmer` optimization and records a parameter history, both driven
  by the shared `nlmixr2est` nlm machinery (not lme4).  Each recorded
  `nlmerSolveGrad()` evaluation logs the population parameter estimate
  (per-subject mean of the `phi` columns) into the resident nlm scale;
  the accumulated history is recovered via `nlmixr2est::nlmGetParHist()`
  and stored on the fit as `parHistData`.  No objective column is shown
  (lme4 owns the deviance).  Iteration printing defaults on
  (`nlmerControl(print = 1L)`).  Requires `nlmixr2est (>= 6.2.0)`.

* The `pseudoOptimControl()` and `fmeMcmcControl()` functions now
  accept either the legacy scalar `print` / `printNcol` / `useColor`
  arguments or a pre-built `nlmixr2est::iterPrintControl()` object via
  `print`.  Internally the control list stores a single
  `iterPrintControl` sub-list (matching the upstream `nlmixr2est`
  unification in `nlmixr2est` PR #651), so iteration output from these
  estimators uses the same shared C++ formatter as every other
  `nlmixr2est` estimator.  Requires `nlmixr2est (>= 6.0.1)`.

* The `iterPrintControl` unification now also covers `nlmerControl()`
  and `saemixControl()`.  `nlmerControl()` gains the standard `print` /
  `printNcol` / `useColor` arguments (or a pre-built
  `nlmixr2est::iterPrintControl()` object) and feeds the resulting
  `iterPrintControl` sub-list to the nlm C solving engine instead of a
  hard-coded `print = 0L`.  `saemixControl()` absorbs its legacy
  `print` (logical), `printNcol` and `useColor` arguments into the same
  `iterPrintControl` sub-list; a nonzero `every` enables the `saemix`
  progress output.

* Fix NONMEM export silently dropping the absorption lag (#190).  A
  `lag(depot)`/`alag(depot)` assignment computed the lag parameter in `$PK`
  but never emitted the corresponding `ALAG<n>=` statement, so NONMEM fit the
  model without any lag.  The lag value is now assigned to `ALAG<n>` in `$PK`.

* NONMEM export now announces when a model variable is renamed because it
  collides with a NONMEM reserved name (e.g. a variable named `alag` becomes
  `RXR1`).  The rename was previously silent (#190).

* Added `nlmer` estimation method: fits nlmixr2 models via `lme4::nlmer` using
  analytical gradients from rxode2 sensitivity equations. Supports
  mu-referenced and non-mu-referenced random-effects models. Access via
  `nlmixr(model, data, est = "nlmer")`. The underlying lme4 fit is stored as
  `fit$nlmer`.

* Fix integer type safety in C++ source: loop variables and size variables now
  use `R_xlen_t` (signed) or `size_t` (unsigned) instead of `int`/`unsigned
  int` where appropriate, preventing potential integer overflow and segfaults
  for vectors with more than 2^31 elements.  The specific crash: in
  `getDvid()`, `int j = cmtDvid.size()` when `cmtDvid.size()` ≥ 2^31 wraps to
  `INT_MIN`, the subsequent decrement jumps to `INT_MAX`, and
  `cmtDvid[INT_MAX]` accesses memory far out of bounds.

* Add bounds check in `popedSolveIdME()` and `popedSolveIdME2()` to verify
  that `modelSwitch` values are within the allocated matrix column dimensions
  (`nend`), in addition to the existing check against the number of unique IDs
  in the global time indexer.

* Remove `qs` since it will be archived and replace with `qs2`.

* Added `saemix` estimation method

# babelmixr2 0.1.10

* Bug fix for the new version of `units` (#179)

# babelmixr2 0.1.9

* Added estimation method `fmeMcmc` which runs `FME::modMCMC()`.  It
  is also compatible with the `coda` package; you can convert with
  `as.mcmc(fit)` and then run coda tools like
  `coda::raftery.diag(coda::as.mcmc(fit2))`.

* Added estimation method `pseudoOptim` which runs
  `FME::pseudoOptim()`. This estimation method requires all parameters
  to be bound.

* Added bug fix for rstudio completion

# babelmixr2 0.1.8

* Maintenance fix for upcoming nlmixr2est and rxode2

# babelmixr2 0.1.7

* Maintenance fix for upcoming PKNCA

# babelmixr2 0.1.6

* Use new nlmixr2est covariate selection enforcement for babelmixr2

* Fix a bug where the NONMEM export isn't working well (#839)

* Check loaded `rxode2` information and compare to what the loaded
  model information should be. This allows better checking of which
  model is loaded and even more robust stability.  It requires
  `rxode2` > `3.0.2`.

# babelmixr2 0.1.5

* Fix bug where `PopED` could error with certain `dvid` values

* Fix bug where if/else clauses in the model could cause the model to
  not predict the values correctly.

* Fix bug so that `shrinkage()` calculation works

* Fix bug so that you can mix 2 different `PopED` data bases in an
  analysis without crashing R.  While this didn't occur with every
  database clash, it more frequently occurred when you interleaved
  `PopED` code between two different `PopED` databases, like in issue
  #131.

* Added a new function `babelBpopIdx(poped.db, "par")` which will get
  the poped index for a model generated from `babelmixr2`, which is
  useful when calculating the power (as in example 11).

# babelmixr2 0.1.4

* Added experimental `PopED` integration

* Removed dependence on `rxode2parse`

* Imported `monolix2rx` from the `monolix2rx` package

* Also allow conversion of a model imported from monolix to a
  `nlmixr2` fit.

# babelmixr2 0.1.3

* Changed default NONMEM rounding protection to FALSE

* Added a `run` option to the `monolixControl()` and `nonemControl()`
  in case you only want to export the modeling files and not run the
  models.

# babelmixr2 0.1.2

* Handle algebraic `mu` expressions

* PKNCA controller now contains `rxControl` since it is used for some
  translation options

* This revision will load the pruned ui model to query the compartment
  properties (i.e. bioavailability, lag time, etc) when writing out the
  NONMEM model.  It should fix issues where the PK block does not
  define some of the variables and will have a larger calculated
  variable that can be used in the model instead.

* When `nonmem2rx` has a different `lst` file, as long as
  `nonmem2rx::nminfo(file)` works, then a successful conversion to a
  `nlmixr2` fit object will occur.

* Fix to save parameter history into `$parHistData` to accommodate
  changes in `focei`'s output (`$parHist` is now derived).

* Changed the solving options to match the new steady state options in
  `rxode2` and how NONMEM implements them.  Also changed the iwres
  model to account for the `rxerr.` instead of the `err.` which was
  updated in `rxode2` as well.


# babelmixr2 0.1.1

* Add new method `as.nlmixr2` to convert `nonmem2rx` methods to `nlmixr` fits

* Dropped `pmxTools` in favor of `nonmem2rx` to conserve some of the
  methods

# babelmixr2 0.1.0

* Babelmixr has support for "monolix", "nonmem", and "pknca" methods
  on release.

* Added a `NEWS.md` file to track changes to the package.
