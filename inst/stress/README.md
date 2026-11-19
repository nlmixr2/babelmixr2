# babelmixr2 NONMEM/Monolix stress test

The stress test checks the edge cases of the nlmixr2 -> NONMEM and
nlmixr2 -> Monolix translations. It has two modes:

- **translate**: write the NONMEM control stream/data and the Monolix
  project/model/data for every case, and check them. No NONMEM or
  Monolix is needed. This mode runs in the package tests.
- **run**: fit every case with NONMEM and/or Monolix end to end, so you
  can see how the translations behave in the real programs. Use this
  mode on a machine that has NONMEM and/or Monolix.

## Quick start: the kit for a NONMEM/Monolix machine

Everything runs from the R session that is set up for NONMEM/Monolix
(for example RStudio); no `Rscript` is needed.

1. In a fresh session (Session > Restart R), install the babelmixr2
   version to test and the development versions of the nlmixr2
   packages it goes with (into the session's `.libPaths()[1]`):

   ```r
   source("https://raw.githubusercontent.com/nlmixr2/babelmixr2/main/inst/stress/install-kit.R")
   installKit()                     # babelmixr2 main
   installKit(ref = "my-branch")    # or a branch, tag or commit
   ```

   Then restart R.

2. Load the kit and check that NONMEM and Monolix are found:

   ```r
   library(babelmixr2)
   source(system.file("stress", "stress.R", package = "babelmixr2"))
   stressCheck()
   stressCheck(nonmem = "/opt/nm75/run/nmfe75")   # if NONMEM is not found
   ```

   NONMEM is found from `options(babelmixr2.nonmem=)`, an `nmfe7*` on
   the `PATH`, or the usual install directories (like
   `/opt/NONMEM/nm75/run/nmfe75` or `C:/nm75/run/nmfe75.bat`); otherwise
   give it with `nonmem=`. Monolix is found through `lixoftConnectors`
   (or give its run command with `monolix=`).

3. Run the kit:

   ```r
   res <- stressKit()                                   # everything that is found
   res <- stressKit(nonmem = "/opt/nm75/run/nmfe75")    # NONMEM not found
   res <- stressKit(engines = "monolix")                # only one engine
   res <- stressKit(cases = "linCmt 1-cmt oral$|rerun") # a few cases first
   res[res$failed, ]                                    # what failed
   ```

   `stressKit()` translates every case (and a sample of nlmixr2lib
   models), fits every case with each engine that was found, compares
   the fits with nlmixr2, and zips the output
   (`babelmixr2-stress-<date>-<time>.zip` in the working directory;
   `attr(res, "zip")` has its path).

4. Send the zip file back (or attach it to an issue at
   <https://github.com/nlmixr2/babelmixr2/issues>).

The full kit is one NONMEM or Monolix fit per case (about 70 fits per
engine), so it takes a while. `stressList()` lists the cases.

`stressKit()` arguments: `nonmem=`, `monolix=`, `engines=`,
`modes=` (`"translate"` and/or `"run"`), `cases=` (a regular
expression), `nlmixr2lib=` (`"none"`, `"sample"`, `"all"`), `out=`
(output directory), `reference=`, `predTol=` (%), `bundle=`.

## Cases

The cases are in `stress.R` (in this directory). They cover:

- `linCmt()` models (1, 2 and 3 compartments; oral, bolus and infusion
  dosing; every parameterization; lag time, bioavailability, and
  modeled rate or duration; covariates; `linCmt()` together with other
  ODEs; parameters that change with time; amounts used in the model).
  NONMEM uses its closed-form `ADVAN1-4`/`ADVAN11-12` with `TRANS1`, and
  Monolix uses `pkmodel()`. When the closed form cannot represent the
  model, both use ODEs instead.
- residual errors (additive, proportional, combined1/2, pow, lognormal,
  logit, Box-Cox, Yeo-Johnson, t distribution, fixed)
- censoring (`CENS`, `CENS` + `LIMIT`, `LIMIT` only)
- model code (`if`/`else`, `else if`, nested `if`/`else` and `ifelse()`
  (pruned, #11), names reserved by NONMEM, long
  dotted names, fixed and block random effects, parameters without
  random effects, initial conditions, time-varying covariates,
  `probitInv()`, between-occasion variability, models that are not
  mu-referenced)
- data (steady state, additional doses, infusions given by rate or by
  duration, modeled rate/duration, reset-and-dose events, missing
  observations, `evid=2` records, character IDs, time not starting at
  zero, extra unused columns, a single subject, oral and iv doses in
  one subject, ODE infusions and steady state infusions)
- models (two endpoints with different residual errors, bounded
  thetas, several covariates, no random effects)
- estimation options (NONMEM `est="imp"`, `"its"`, `"posthoc"`,
  `cov=""`, `advanOde="advan6"`, `linCmt="ode"`; Monolix
  `useLinearization=TRUE`, `linCmt="ode"`, `stiff=TRUE`,
  `variability="decreasing"`)
- a second fit of the same model, which should read the saved
  NONMEM/Monolix output instead of running again

In run mode every fit is also checked:

- the objective function and the estimates are finite
- rxode2 reproduces the NONMEM/Monolix individual predictions: the
  median relative difference of `IPRED` must be at most `--pred-tol`
  percent (default 5)
- with `--reference`, the largest relative difference from the
  nlmixr2 estimates is reported (not a failure: the methods differ)

Each case says, for each engine, whether it should translate or be
refused with a documented error. A case can also list text that must
appear in the control stream or model file (for example `ADVAN4 TRANS1`
or `pkmodel(`).

Optionally the test also translates the models in
[nlmixr2lib](https://nlmixr2.github.io/nlmixr2lib/) (a sample, or all
of them). Those models must either translate or be refused with a
documented error.

## Requirements

On the machine that runs the stress test:

- R with babelmixr2, nlmixr2est, rxode2 and nlmixr2data installed (plus
  nlmixr2lib for the nlmixr2lib models); `installKit()` installs them
  (see the quick start).
- For run mode with NONMEM: a NONMEM installation and its run
  command (like `nmfe75`), either on the `PATH` or given with its full
  path.
- For run mode with Monolix: Monolix plus the `lixoftConnectors` R
  package that comes with it (see Lixoft's documentation), or a Monolix
  command-line run command.

The closed-form `linCmt()` translations need an rxode2 that has
`rxode2::linCmtMicro()`. With an older rxode2, every `linCmt()` model is
translated to ODEs, so the closed-form checks (`ADVAN2 TRANS1`,
`pkmodel(`) fail.

## Running it with Rscript

Where `Rscript` works with the right library paths, `run-stress.R`
does the same from a shell (`--kit`, `--check`, `--list` and the
options below match the `stressKit()` arguments).

Find the runner script:

```sh
STRESS=$(Rscript -e 'cat(system.file("stress", "run-stress.R", package="babelmixr2"))')
```

Translate only (no NONMEM/Monolix needed):

```sh
Rscript "$STRESS" --mode=translate
```

Run NONMEM end to end:

```sh
Rscript "$STRESS" --mode=run --engine=nonmem --nonmem=nmfe75
# or with the full path
Rscript "$STRESS" --mode=run --engine=nonmem --nonmem=/opt/nm75/run/nmfe75
```

Run Monolix end to end (uses `lixoftConnectors` when it is installed):

```sh
Rscript "$STRESS" --mode=run --engine=monolix
```

Both engines, and compare the estimates with nlmixr2 fits of the same
models (`focei` for NONMEM, `saem` for Monolix):

```sh
Rscript "$STRESS" --mode=run --nonmem=nmfe75 --reference
```

Other options:

```sh
Rscript "$STRESS" --check                         # versions; is NONMEM/Monolix found?
Rscript "$STRESS" --list                          # list the cases
Rscript "$STRESS" --pred-tol=1                    # stricter IPRED check (%)
Rscript "$STRESS" --bundle                        # zip the output directory
Rscript "$STRESS" --cases='linCmt'                # only some cases (regex)
Rscript "$STRESS" --nlmixr2lib=sample             # add nlmixr2lib models
Rscript "$STRESS" --nlmixr2lib=all                # every nlmixr2lib model (slow)
Rscript "$STRESS" --out=my-stress-results         # output directory
```

From an R session the same is `stressKit()` (see the quick start), for
example `stressKit(engines = "nonmem", modes = "run", nonmem = "nmfe75")`.

The run can take a while: each case is a full NONMEM or Monolix fit.
Use `--cases=` to start with a few cases.

## Output

The output directory (by default `babelmixr2-stress-<date>-<time>`) has:

- `results.csv`: one row per case and engine, with these columns:
  - `status`:
    - `ok`: translated (and fit in run mode)
    - `refused`: refused with the expected error
    - `error`: an unexpected error
    - `problem`: the files were written but a check failed
    - `not refused`: expected a refusal but it translated
    - `skipped`: the case is translation only
    - `upstream`: a known problem in another nlmixr2 package (rxode2 or
      nlmixr2est), not in babelmixr2 (not counted as a failure)
  - `message`/`problems`: what went wrong
  - `seconds`: how long it took
  - `objf`: the objective function (run mode)
  - `ipredRelDiff`/`predRelDiff`: median relative difference (%)
    between the rxode2 and the NONMEM/Monolix `IPRED`/`PRED`
  - `rerunSeconds`: how long the second fit took (the rerun case)
  - `maxRelDiffTheta`: the largest relative difference from the nlmixr2
    estimates (with `--reference`)
- `summary.md`: a summary table and the failures.
- `sessionInfo.txt`: the R session (package versions).
- `nonmem/<case>/` and `monolix/<case>/`: the control streams, data,
  and NONMEM/Monolix output for each case, to look at a failure in
  detail (under `translate/` and `run/` with `--kit`). `fit.log` in
  each case's directory has what NONMEM/Monolix printed (like Monolix's
  `[ERROR]` lines).

The simulated data are rounded to 8 significant digits, so every
machine writes the same data; the saved NONMEM/Monolix output in a
returned zip can then be read again (without running NONMEM/Monolix)
on another machine.

The script exits with status 1 when any case fails, so it can also be
used in a CI job on a machine with NONMEM or Monolix.

## Package tests

The translate mode runs in the babelmixr2 tests
(`tests/testthat/test-nonmem-monolix-stress.R` and
`tests/testthat/test-nonmem-monolix-nlmixr2lib.R`). These are part of
slow-test batch 3 (`BABELMIXR2_TEST_BATCH=3`). Set
`BABELMIXR2_STRESS_ALL=true` to test every nlmixr2lib model instead of a
sample.

Please report failures from run mode (with `summary.md` and the case's
output directory) at <https://github.com/nlmixr2/babelmixr2/issues>.
