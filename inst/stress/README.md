# babelmixr2 NONMEM/Monolix stress test

The stress test checks the edge cases of the nlmixr2 -> NONMEM and
nlmixr2 -> Monolix translations. It has two modes:

- **translate**: write the NONMEM control stream/data and the Monolix
  project/model/data for every case, and check them. No NONMEM or
  Monolix is needed. This mode runs in the package tests.
- **run**: fit every case with NONMEM and/or Monolix end to end, so you
  can see how the translations behave in the real programs. Use this
  mode on a machine that has NONMEM and/or Monolix.

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
- model code (`if`/`else`, `else if`, names reserved by NONMEM, long
  dotted names, fixed and block random effects, parameters without
  random effects, initial conditions, time-varying covariates,
  `probitInv()`, between-occasion variability, models that are not
  mu-referenced)
- data (steady state, additional doses, infusions given by rate or by
  duration, modeled rate/duration, reset-and-dose events)

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
  nlmixr2lib for `--nlmixr2lib=`). Install the babelmixr2 version you
  want to test, for example from GitHub:

  ```r
  remotes::install_github("nlmixr2/babelmixr2")
  ```

- For `--mode=run` with NONMEM: a NONMEM installation and its run
  command (like `nmfe75`), either on the `PATH` or given with its full
  path.
- For `--mode=run` with Monolix: Monolix plus the `lixoftConnectors` R
  package that comes with it (see Lixoft's documentation), or a Monolix
  command-line run command.

The closed-form `linCmt()` translations need an rxode2 that has
`rxode2::linCmtMicro()`. With an older rxode2, every `linCmt()` model is
translated to ODEs, so the closed-form checks (`ADVAN2 TRANS1`,
`pkmodel(`) fail.

## Running it

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
Rscript "$STRESS" --list                          # list the cases
Rscript "$STRESS" --cases='linCmt'                # only some cases (regex)
Rscript "$STRESS" --nlmixr2lib=sample             # add nlmixr2lib models
Rscript "$STRESS" --nlmixr2lib=all                # every nlmixr2lib model (slow)
Rscript "$STRESS" --out=my-stress-results         # output directory
```

From R, the same can be done with:

```r
source(system.file("stress", "stress.R", package="babelmixr2"))
res <- stressRun(stressCases(), engines="nonmem", mode="run",
                 dir="stress-out", runCommand=list(nonmem="nmfe75"),
                 progress=TRUE)
res[stressFailed(res), ]
```

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
  - `maxRelDiffTheta`: the largest relative difference from the nlmixr2
    estimates (with `--reference`)
- `summary.md`: a summary table and the failures.
- `nonmem/<case>/` and `monolix/<case>/`: the control streams, data,
  and NONMEM/Monolix output for each case, to look at a failure in
  detail.

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
