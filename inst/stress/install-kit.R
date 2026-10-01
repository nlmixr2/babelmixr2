# Install the babelmixr2 version to stress test, with the development
# versions of the nlmixr2 packages it is tested with.
#
# From an R session (like RStudio), source this file and call
# installKit() in a fresh session (Session > Restart R) so no nlmixr2
# package is loaded while it is replaced; see the quick start in
# README.md.  The packages go into the session's library
# (.libPaths()[1]).  From a shell: Rscript install-kit.R [--ref=REF]
# [--cran].
#
# lixoftConnectors is not installed here: it comes with Monolix (see
# README.md).

#' Install babelmixr2 and the nlmixr2 packages to stress test
#'
#' @param ref babelmixr2 git branch, tag or commit
#' @param cran use the CRAN versions of the other nlmixr2 packages
#'   instead of their GitHub versions
#' @param lib library to install into
#' @return nothing
#' @noRd
installKit <- function(ref = "main", cran = FALSE, lib = .libPaths()[1]) {
  if (identical(unname(getOption("repos")["CRAN"]), "@CRAN@")) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
  if (!requireNamespace("remotes", quietly = TRUE)) {
    utils::install.packages("remotes", lib = lib)
  }
  # dependency order: each package is built against the ones before it
  .gh <- c(
    "nlmixr2/lotri",
    "nlmixr2/rxode2ll",
    "nlmixr2/rxode2lincmt",
    "nlmixr2/rxode2",
    "nlmixr2/nlmixr2est",
    "nlmixr2/nlmixr2extra",
    "nlmixr2/nlmixr2plot",
    "nlmixr2/nonmem2rx",
    "nlmixr2/monolix2rx"
  )
  utils::install.packages(c("withr", "nlmixr2data", "nlmixr2lib"), lib = lib)
  if (!cran) {
    for (.p in .gh) {
      remotes::install_github(
        .p,
        upgrade = "never",
        dependencies = TRUE,
        lib = lib
      )
    }
  }
  remotes::install_github(
    paste0("nlmixr2/babelmixr2@", ref),
    upgrade = "never",
    dependencies = TRUE,
    force = TRUE,
    lib = lib
  )
  message(
    "\ninstalled into ",
    lib,
    "; restart R, then check NONMEM/Monolix with:"
  )
  message(
    "  source(system.file(\"stress\", \"stress.R\", package = \"babelmixr2\"))"
  )
  message("  stressCheck()")
  invisible()
}

if (!interactive() && sys.nframe() == 0L) {
  .args <- commandArgs(trailingOnly = TRUE)
  .ref <- sub("^--ref=", "", grep("^--ref=", .args, value = TRUE)[1])
  installKit(
    ref = if (is.na(.ref)) "main" else .ref,
    cran = "--cran" %in% .args
  )
}
