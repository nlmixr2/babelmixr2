# babelmixr2 0.1.11.9000

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
