#!/usr/bin/env Rscript
# Install the babelmixr2 version to stress test, with the development
# versions of the nlmixr2 packages it is tested with.
#
# Usage:
#   Rscript install-kit.R [--ref=main] [--cran]
#
#   --ref=REF   babelmixr2 git branch, tag or commit (default: main)
#   --cran      use the CRAN versions of the other nlmixr2 packages
#               instead of their GitHub versions
#
# lixoftConnectors is not installed here: it comes with Monolix (see
# README.md).

.args <- commandArgs(trailingOnly = TRUE)
.ref <- sub("^--ref=", "", grep("^--ref=", .args, value = TRUE)[1])
if (is.na(.ref)) {
  .ref <- "main"
}
.cran <- "--cran" %in% .args

.repos <- getOption("repos")
if (is.null(.repos) || identical(unname(.repos["CRAN"]), "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
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
install.packages(c("withr", "nlmixr2data", "nlmixr2lib"))
if (!.cran) {
  for (.p in .gh) {
    remotes::install_github(.p, upgrade = "never", dependencies = TRUE)
  }
}
remotes::install_github(
  paste0("nlmixr2/babelmixr2@", .ref),
  upgrade = "never",
  dependencies = TRUE,
  force = TRUE
)

message("\ninstalled; now check that NONMEM/Monolix are found with:")
message(
  "  Rscript \"",
  system.file("stress", "run-stress.R", package = "babelmixr2"),
  "\" --check"
)
