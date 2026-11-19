# babelmixr2 0.1.12

* Fixes R-devel's `R CMD check` deprecation of `structure(.Label=)`
  (an empty Monolix administration table now uses `factor(levels=)`).
* Package code and examples work without Suggests installed (noSuggests
  check): argument checks use `checkmate::assert_*()` and `units` is
  requested with `rxode2::rxReq()`.
* Many NONMEM/Monolix translation fixes and improvements to the
  `saemix`, `nlmer`, `fmeMcmc` and `pseudoOptim` estimation methods;
  see NEWS.md.

## R CMD check results

0 errors | 0 warnings | 1 note (local R 4.6.1, Ubuntu 24.04, `--as-cran`)

* Suggests not in mainstream repositories: `lixoftConnectors`.  It ships
  with Monolix, is optional, and is only used when available (as in
  previous releases).

## Reverse dependencies

No exported functions were removed (only new exports were added).
nlmixr2, nlmixr2save and ubiquity pass `R CMD check` with this version.
