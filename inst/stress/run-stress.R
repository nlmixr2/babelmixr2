#!/usr/bin/env Rscript
# babelmixr2 NONMEM/Monolix stress test runner
#
# Translates (and optionally fits) every stress case with NONMEM and/or
# Monolix and writes a report.  See README.md in this directory.  The
# same can be run from an R session (like RStudio) with stressCheck(),
# stressList() and stressKit() after sourcing stress.R.
#
# Usage:
#   Rscript run-stress.R [options]
#
# Options:
#   --kit                     everything for a NONMEM/Monolix machine:
#                             translate every case (and a sample of
#                             nlmixr2lib), fit every case with the
#                             engines that are found (--reference too),
#                             and zip the output (--bundle)
#   --check                   report the versions and whether NONMEM and
#                             Monolix are found, then exit
#   --engine=nonmem,monolix   engines to test (default: both; with --kit,
#                             the engines that are found)
#   --mode=translate|run      translate: only write the model/data files
#                             run: fit with NONMEM/Monolix end to end
#                             (default: translate)
#   --cases=REGEX             only the cases whose name matches REGEX
#   --nlmixr2lib=none|sample|all
#                             also translate models from nlmixr2lib
#                             (default: none; always translation only)
#   --out=DIR                 output directory (default:
#                             babelmixr2-stress-YYYYMMDD-HHMMSS)
#   --nonmem=COMMAND          command that runs NONMEM, like nmfe75
#                             (default: option babelmixr2.nonmem, then
#                             an nmfe7* found on the PATH or in the
#                             usual install directories)
#   --monolix=COMMAND         command that runs Monolix (default: option
#                             babelmixr2.monolix, or lixoftConnectors
#                             when it is installed)
#   --reference               in run mode, also fit each model with
#                             nlmixr2 (focei for NONMEM, saem for
#                             Monolix) and compare the estimates
#   --pred-tol=PERCENT        in run mode, the largest median relative
#                             difference between the rxode2 and the
#                             NONMEM/Monolix IPRED (default: 5)
#   --bundle                  zip the output directory to send back
#   --list                    list the cases and exit

.args <- commandArgs(trailingOnly = TRUE)
.opt <- function(name, default = NULL) {
  .w <- grep(paste0("^--", name, "(=|$)"), .args, value = TRUE)
  if (length(.w) == 0L) {
    return(default)
  }
  if (!grepl("=", .w[1])) {
    return(TRUE)
  }
  sub(paste0("^--", name, "="), "", .w[1])
}

# loading babelmixr2 registers the "nonmem" and "monolix" estimation
# methods with nlmixr2est
suppressPackageStartupMessages(loadNamespace("babelmixr2"))

.file <- system.file("stress", "stress.R", package = "babelmixr2")
if (.file == "") {
  # running from a source checkout
  .file <- file.path(
    dirname(sub(
      "^--file=",
      "",
      grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
    )),
    "stress.R"
  )
}
source(.file)

.split <- function(x) if (is.null(x)) NULL else strsplit(x, ",")[[1]]
.nonmem <- .opt("nonmem")
.monolix <- .opt("monolix")

if (isTRUE(.opt("check", FALSE))) {
  stressCheck(nonmem = .nonmem, monolix = .monolix)
  quit(save = "no", status = 0)
}
if (isTRUE(.opt("list", FALSE))) {
  stressList(cases = .opt("cases"), nlmixr2lib = .opt("nlmixr2lib", "none"))
  quit(save = "no", status = 0)
}

.out <- .opt(
  "out",
  paste0("babelmixr2-stress-", format(Sys.time(), "%Y%m%d-%H%M%S"))
)
.predTol <- as.numeric(.opt("pred-tol", 5))
.res <- if (isTRUE(.opt("kit", FALSE))) {
  stressKit(
    nonmem = .nonmem,
    monolix = .monolix,
    engines = .split(.opt("engine")),
    cases = .opt("cases"),
    nlmixr2lib = .opt("nlmixr2lib"),
    out = .out,
    predTol = .predTol
  )
} else {
  stressKit(
    nonmem = .nonmem,
    monolix = .monolix,
    engines = .split(.opt("engine", "nonmem,monolix")),
    modes = .opt("mode", "translate"),
    cases = .opt("cases"),
    nlmixr2lib = .opt("nlmixr2lib", "none"),
    out = .out,
    reference = isTRUE(.opt("reference", FALSE)),
    predTol = .predTol,
    bundle = isTRUE(.opt("bundle", FALSE))
  )
}
quit(save = "no", status = ifelse(any(.res$failed), 1L, 0L))
