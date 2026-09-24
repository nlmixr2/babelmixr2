# babelmixr2 0.1.11.9000

* `babelmixr2` can now export a mu-referenced `nlmixr2` model to PharmML 0.9
  with `as.pharmml()`.  PharmML is a standard XML description of a
  pharmacometric model -- parameter model, structural model, observation
  model, trial design and estimation step in one self-describing document --
  so this makes an `nlmixr2` model archivable and exchangeable without the
  receiving tool needing to understand `rxode2` syntax.  Unlike the NONMEM and
  Monolix backends this is a writer rather than an estimation method, so there
  is no `est="pharmml"`.

  Supporting functions: `rxToPharmml()` translates a single expression (the
  counterpart of `rxToNonmem()` and `rxToMonolix()`), `pharmmlControl()`
  carries the options, and `pharmmlValidate()` checks a document against the
  schema.  The schemas are vendored in the package because the upstream host
  no longer serves them, so validation is entirely offline and never reaches
  the network.

  Solved (`linCmt()`) models are exported as PharmML PK macros rather than
  being expanded to ODEs, which preserves the structure the model was written
  in.  This is the one thing the PharmML writer can do that the NONMEM and
  Monolix writers cannot -- both of those refuse `linCmt()` outright.

  Categorical covariates are detected from the dataset: a character or factor
  column becomes a PharmML `Categorical` covariate with a category per level,
  and the column mapping records how the exported numeric codes correspond to
  those levels.

  `lnorm()` residual errors are written as a `log` transformation of both
  sides of the observation model.

  Models that PharmML cannot express -- non-normal residuals, power residual
  error, `boxCox()`, `yeoJohnson()` and logit/probit residual
  transformations, inter-occasion variability, mixture models, Michaelis-Menten or
  transit absorption through `linCmt()` -- raise an error naming the construct
  rather than emitting a document that looks plausible but is wrong.  Every
  document is validated against the schema before it is returned.

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
