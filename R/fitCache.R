#' Load the cached nlmixr2 fit from a run directory
#'
#' The cache is written with [saveRDS()].  Directories produced by
#' babelmixr2 versions that cached with `qs2` hold a `.qs2` file
#' instead; that is read once (when `qs2` is available) and rewritten as
#' the `.rds` so later runs no longer need `qs2`.
#'
#' @param rds path of the `.rds` cache file
#'
#' @return the cached fit, or `NULL` when there is no usable cache and
#'   the fit has to be rebuilt from the run output
#'
#' @noRd
#' @author Matthew L. Fidler
.babelmixr2LoadFitCache <- function(rds) {
  if (file.exists(rds)) {
    return(readRDS(rds))
  }
  .legacy <- sub("\\.rds$", ".qs2", rds)
  if (file.exists(.legacy) && requireNamespace("qs2", quietly=TRUE)) {
    .ret <- try(getExportedValue("qs2", "qs_read")(.legacy), silent=TRUE)
    if (!inherits(.ret, "try-error")) {
      .babelmixr2SaveFitCache(.ret, rds)
      return(.ret)
    }
  }
  NULL
}

#' Save the nlmixr2 fit cache for a run directory
#'
#' @param fit nlmixr2 fit object to cache
#'
#' @param rds path of the `.rds` cache file
#'
#' @return `fit`, invisibly
#'
#' @noRd
#' @author Matthew L. Fidler
.babelmixr2SaveFitCache <- function(fit, rds) {
  saveRDS(fit, rds)
  invisible(fit)
}
