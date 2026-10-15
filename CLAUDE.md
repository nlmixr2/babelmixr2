# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Overview

**babelmixr2** lets users fit an `nlmixr2`/`rxode2` model with other
estimation and design software through the usual
[`nlmixr2()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.html)
interface: NONMEM, Monolix, PKNCA, PopED, saemix,
[`lme4::nlmer`](https://rdrr.io/pkg/lme4/man/nlmer.html), FME MCMC, and
a pseudo-optimizer. It also imports finished NONMEM/Monolix runs (via
`nonmem2rx`/`monolix2rx`) as nlmixr2 fit objects with
[`as.nlmixr2()`](https://nlmixr2.github.io/babelmixr2/reference/as.nlmixr2.md).
The package is mostly R with a small C/C++ layer (PopED solving and data
conversion). Version and dependencies are in `DESCRIPTION` (testthat
edition 3, roxygen with markdown).

## Build and Development Commands

``` r

devtools::install()   # or from the shell: R CMD INSTALL .
devtools::document()  # regenerate man/ and NAMESPACE from roxygen
```

`configure`/`configure.win` run `inst/tools/workaround.R`, which writes
`src/Makevars` (or `src/Makevars.win`) from `src/Makevars.in` and strips
the `RcppArmadillo.h` include from `src/RcppExports.cpp`. Do not edit
the generated `src/Makevars*` by hand; edit `src/Makevars.in`.

### Tests

``` sh
NOT_CRAN=true Rscript -e "devtools::test()"
NOT_CRAN=true Rscript -e "devtools::test(filter='nonmem')"   # file name without test- and .R
```

`tests/testthat.R` partitions the suite:

- Locally and on CRAN everything runs (estimation files skip themselves
  on CRAN with `skip_on_cran()`).
- On CI (`CI=true`, not CRAN), R-CMD-check runs only the fast subset:
  all files except those in `.slowBatches`.
- `BABELMIXR2_TEST_BATCH=<n>` runs only slow batch `n` (1 = `saemix`; 2
  = `nonmem-read`, `monolix-read`, `poped`, `as-nlmixr2`; 3 =
  `nonmem-monolix-stress`, `nonmem-monolix-nlmixr2lib`). The weekly
  `.github/workflows/slow-tests.yaml` matrix must match the number of
  batches; update both together.
- `BABELMIXR2_STRESS_ALL=true` makes the nlmixr2lib sweep use every
  model.
- On CI/CRAN tests run serially, and on CRAN rxode2 threads are capped
  at 2.

### R CMD check

``` r

devtools::check()   # or R CMD build . && R CMD check babelmixr2_*.tar.gz
```

CI (`.github/workflows/R-CMD-check.yaml`) checks macOS, Windows and
Ubuntu (release, devel, oldrel-1). Other workflows:
`test-coverage.yaml`, `pkgdown.yaml`, `rhub.yaml`, `slow-tests.yaml`.

## Architecture

### Estimation methods are nlmixr2est S3 plug-ins

Each method is wired into nlmixr2est through S3 methods on the `est`
class: `nlmixr2Est.<est>`, `getValidNlmixrCtl.<est>`,
`nmObjGetControl.<est>`, `nmObjHandleControlObject.<est>Control`, and
`rxUiGet.<name>` accessors on the model UI. Methods for `fmeMcmc` and
`pseudoOptim` (and
[`coda::as.mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html)) are
registered at load time in `R/zzz.R` with
[`rxode2::.s3register()`](https://nlmixr2.github.io/rxode2/reference/dot-s3register.html);
the rest are in `NAMESPACE`.

| Method | Main files |
|----|----|
| NONMEM | `R/nonmem*.R` (`nonmemControl.R`, `nonmemNlmixr2est.R`, `nonmemRxUiGet*.R` build each `$` record, `nonmemReadData.R`, `nonmemFailure.R`) |
| Monolix | `R/monolix*.R` (`monolixRxUiGetMlxtran*.R` build each mlxtran section, `monolixReadData.R`, `monolixPriors.R`) |
| PKNCA | `R/pknca.R` |
| PopED | `R/poped.R`, `src/poped.cpp` |
| saemix | `R/saemix.R`, `R/saemixControl.R` |
| nlmer | `R/nlmer.R`, `R/nlmerControl.R`, `R/nlmerRxUiGet.R` |
| FME MCMC | `R/fmeMcmc.R` |
| pseudoOptim | `R/pseudoOptim.R` |

Shared pieces:

- `R/convert.R`
  ([`bblDatToNonmem()`](https://nlmixr2.github.io/babelmixr2/reference/bblDatToMonolix.md),
  [`bblDatToMonolix()`](https://nlmixr2.github.io/babelmixr2/reference/bblDatToMonolix.md))
  and `src/transform.cpp` convert nlmixr2 data/DV for the external
  tools.
- `R/nmGetDistribution*Lines.R` write residual-error/distribution lines.
- `R/pruneIf.R` prunes `if`/`else` blocks for tools that cannot
  reassign.
- `R/linCmt.R` handles `linCmt()` models; `R/units.R` handles units.
- `R/fitCache.R` caches fits; `R/sigdigControl.R` maps `sigdig` to
  tolerances.
- `R/as.nlmixr2*.R` convert `nonmem2rx`/`monolix2rx` imports to fits.
- `inst/*.txt` are NONMEM code templates (`$ERROR`, `CCONTR`) read with
  [`system.file()`](https://rdrr.io/r/base/system.file.html) for the
  error/transform variants.
- `inst/poped/` holds PopED example scripts; `inst/stress/` is the
  NONMEM/Monolix translation stress kit (see `inst/stress/README.md`).

### External programs

NONMEM and Monolix are found via options such as
`options(babelmixr2.nonmem=)` and `options(babelmixr2.monolix=)` (see
[`nonmemControl()`](https://nlmixr2.github.io/babelmixr2/reference/nonmemControl.md)/[`monolixControl()`](https://nlmixr2.github.io/babelmixr2/reference/monolixControl.md)).
Monolix can also run through the optional `lixoftConnectors` package.
Tests that need these programs must skip when they are absent;
translation-only tests must not need them.

### Compiled code and the rxode2 pointer table

- `src/init.c` is the hand-written `.Call` registration table with
  hardcoded argument counts. After adding or changing an
  `// [[Rcpp::export]]` function, run `Rcpp::compileAttributes(".")` AND
  update the prototype and arity in `src/init.c`.
- babelmixr2 links to rxode2 through rxode2’s function-pointer table:
  `R/zzz.R` records `rxode2.api` at build time and `.iniRxode2Ptr()`
  passes
  [`rxode2::.rxode2ptrs()`](https://nlmixr2.github.io/rxode2/reference/dot-rxode2ptrs.html)
  to `_babelmixr2_iniRxodePtrs`. `src/Makevars.in` passes the table
  length as `-DRXAPI`. A newer rxode2 API requires reinstalling
  babelmixr2.

### Generated files (do not edit manually)

- `R/RcppExports.R`, `src/RcppExports.cpp` –
  [`Rcpp::compileAttributes()`](https://rdrr.io/pkg/Rcpp/man/compileAttributes.html)
- `man/*.Rd`, `NAMESPACE` –
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
- `README.md` – knit from `README.Rmd`

## Code Style

- Exported functions: `camelCase` (e.g. `nonmemControl`,
  `bblDatToNonmem`).
- Internal functions and local variables: `.camelCase` with a leading
  dot.
- S3 methods follow `generic.class` (e.g. `rxUiGet.monolixOmega`).
- Formatting uses air (`air.toml`, line width 120, indent 2). Model DSL
  calls (`ini`, `model`, `rxode2`, `quote`, `lotri`, …) are skipped so
  `d/dt(...)` and `x <- 1; label("...")` are left as written.
- Do not use `:::` in package code (CRAN flags it). Tests run in the
  package namespace, so internals can be called by bare name.
- Use `checkmate::assert_*()` (not `expect_*()`) for argument checks,
  and
  [`rxode2::rxReq()`](https://nlmixr2.github.io/rxode2/reference/rxReq.html)
  for suggested packages; package code must work without Suggests
  installed (CRAN noSuggests check).
- Use American English. ASCII only everywhere (CRAN): `--`, `->`,
  straight quotes, `...`.
- `NEWS.md` has one `# babelmixr2 X.Y.Z` heading per version followed by
  bullet entries describing the user-visible change.

## Testing Conventions

- Test files are `tests/testthat/test-*.R`; reference NONMEM/Monolix
  runs are stored as zip files and directories next to them.
- Guard optional packages with `skip_if_not_installed()` and slow or
  estimation-heavy tests with `skip_on_cran()`.
- A new slow test file belongs in a `.slowBatches` entry in
  `tests/testthat.R`, otherwise it runs in every push/PR CI check.

## Release / CRAN Cautions

- `cran-comments.md`, `CRAN-SUBMISSION`, `benchmarks/`, `pkgdown/`,
  `vignettes/articles/` and similar are excluded via `.Rbuildignore`;
  keep new non-package files there too.
- Vignettes live only in `vignettes/articles/` (pkgdown articles, not
  built on CRAN).
- Examples and CRAN tests must not need NONMEM, Monolix or unlisted
  packages.

## CRAN tests: never assert on timing

- Never put a wall-clock timing expectation in a test that can run on
  CRAN: no `expect_lt()`/`expect_true()` on
  [`Sys.time()`](https://rdrr.io/r/base/Sys.time.html) differences,
  [`system.time()`](https://rdrr.io/r/base/system.time.html),
  [`proc.time()`](https://rdrr.io/r/base/proc.time.html),
  `bench`/`microbenchmark` results, or any other elapsed-time threshold.
  CRAN check machines are shared and heavily loaded, so elapsed time
  there is random and a timing assertion fails at random. Example:
  nlmixr2plot `test-model-diagram.R` asserted a `modelGraph()` call took
  \< 30 s; on CRAN it took 48.4 s and failed.
- Test the result, not the speed. If a speed regression guard is really
  needed, put it only in a test that starts with `skip_on_cran()` (or
  the package’s own CRAN skip helper), and never in examples or
  vignettes.
- Tests that take more than a few seconds should also skip on CRAN.
