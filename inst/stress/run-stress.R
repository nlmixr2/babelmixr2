#!/usr/bin/env Rscript
# babelmixr2 NONMEM/Monolix stress test runner
#
# Translates (and optionally fits) every stress case with NONMEM and/or
# Monolix and writes a report.  See README.md in this directory.
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

.nonmemCmd <- .opt("nonmem")
if (is.null(.nonmemCmd)) {
  .nonmemCmd <- stressFindNonmem()
}
.monolixCmd <- .opt("monolix")
.monolix <- if (is.null(.monolixCmd)) {
  stressMonolixStatus()
} else {
  paste0("command ", .monolixCmd)
}

if (isTRUE(.opt("check", FALSE))) {
  message(paste(stressVersions(), collapse = "\n"))
  message(
    "- NONMEM: ",
    ifelse(
      nzchar(.nonmemCmd),
      .nonmemCmd,
      "not found (use --nonmem=/path/to/nmfe75)"
    )
  )
  message(
    "- Monolix: ",
    ifelse(
      nzchar(.monolix),
      .monolix,
      "not found (install lixoftConnectors or use --monolix=)"
    )
  )
  if (!("linCmtMicro" %in% getNamespaceExports("rxode2"))) {
    message(
      "- rxode2 has no linCmtMicro(): the closed-form linCmt() checks will fail"
    )
  }
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------
# options
# ---------------------------------------------------------------------

.kit <- isTRUE(.opt("kit", FALSE))
if (.kit) {
  .engines <- .opt("engine")
  if (is.null(.engines)) {
    .engines <- c(
      if (nzchar(.nonmemCmd)) "nonmem",
      if (nzchar(.monolix)) "monolix"
    )
    if (length(.engines) == 0L) {
      stop("--kit found neither NONMEM nor Monolix; see --check", call. = FALSE)
    }
  } else {
    .engines <- strsplit(.engines, ",")[[1]]
  }
  .modes <- c("translate", "run")
  .lib <- .opt(
    "nlmixr2lib",
    if (requireNamespace("nlmixr2lib", quietly = TRUE)) "sample" else "none"
  )
  .reference <- TRUE
  .bundle <- TRUE
} else {
  .engines <- strsplit(.opt("engine", "nonmem,monolix"), ",")[[1]]
  .modes <- match.arg(.opt("mode", "translate"), c("translate", "run"))
  .lib <- .opt("nlmixr2lib", "none")
  .reference <- isTRUE(.opt("reference", FALSE))
  .bundle <- isTRUE(.opt("bundle", FALSE))
}
.engines <- match.arg(.engines, c("nonmem", "monolix"), several.ok = TRUE)
.lib <- match.arg(.lib, c("none", "sample", "all"))
.predTol <- as.numeric(.opt("pred-tol", 5))
.out <- .opt(
  "out",
  paste0("babelmixr2-stress-", format(Sys.time(), "%Y%m%d-%H%M%S"))
)

.cases <- stressCases()
.libCases <- list()
if (.lib != "none") {
  if (!requireNamespace("nlmixr2lib", quietly = TRUE)) {
    stop("--nlmixr2lib needs the nlmixr2lib package", call. = FALSE)
  }
  .libCases <- stressLibCases(all = (.lib == "all"))
}
.re <- .opt("cases")
if (!is.null(.re)) {
  .cases <- .cases[grepl(.re, names(.cases))]
  .libCases <- .libCases[grepl(.re, names(.libCases))]
}
if (isTRUE(.opt("list", FALSE))) {
  for (.c in c(.cases, .libCases)) {
    .e <- function(engine) {
      if (!(engine %in% .c$engines)) {
        return("-")
      }
      ifelse(
        .c$expect[[engine]] %in% c("ok", "any"),
        .c$expect[[engine]],
        "refuse"
      )
    }
    cat(sprintf(
      "%-50s nonmem=%-8s monolix=%-8s %s\n",
      .c$name,
      .e("nonmem"),
      .e("monolix"),
      .c$description
    ))
  }
  quit(save = "no", status = 0)
}

.runCommand <- list(
  nonmem = if (nzchar(.nonmemCmd)) .nonmemCmd,
  monolix = .monolixCmd
)
if ("run" %in% .modes) {
  if ("nonmem" %in% .engines && !nzchar(.nonmemCmd)) {
    stop(
      "run mode needs --nonmem=COMMAND (like --nonmem=nmfe75) ",
      "or options(babelmixr2.nonmem=)",
      call. = FALSE
    )
  }
  if ("monolix" %in% .engines && !nzchar(.monolix)) {
    stop(
      "run mode needs lixoftConnectors, --monolix=COMMAND ",
      "or options(babelmixr2.monolix=)",
      call. = FALSE
    )
  }
}

# ---------------------------------------------------------------------
# running
# ---------------------------------------------------------------------

dir.create(.out, showWarnings = FALSE, recursive = TRUE)
.out <- normalizePath(.out)
writeLines(
  utils::capture.output(utils::sessionInfo()),
  file.path(.out, "sessionInfo.txt")
)
message(paste(stressVersions(), collapse = "\n"))
if ("run" %in% .modes) {
  if ("nonmem" %in% .engines) {
    message("- NONMEM: ", .nonmemCmd)
  }
  if ("monolix" %in% .engines) message("- Monolix: ", .monolix)
}
message(
  length(.cases) + length(.libCases),
  " cases; engines: ",
  paste(.engines, collapse = ", "),
  "; mode: ",
  paste(.modes, collapse = ", "),
  "; output: ",
  .out
)

.res <- list()
for (.mode in .modes) {
  # the nlmixr2lib models are translation only
  .c <- if (.mode == "translate") c(.cases, .libCases) else .cases
  .dir <- if (length(.modes) > 1L) file.path(.out, .mode) else .out
  .res[[.mode]] <- stressRun(
    .c,
    engines = .engines,
    mode = .mode,
    dir = .dir,
    runCommand = .runCommand,
    reference = .reference,
    predTol = .predTol,
    progress = TRUE
  )
}
.res <- do.call(rbind, .res)
rownames(.res) <- NULL
.res$failed <- stressFailed(.res)
utils::write.csv(.res, file.path(.out, "results.csv"), row.names = FALSE)

.md <- stressSummary(.res, .modes, .engines, .nonmemCmd, .monolix)
writeLines(.md, file.path(.out, "summary.md"))
message("\n", paste(.md, collapse = "\n"))
message("\nresults: ", file.path(.out, "results.csv"))

if (.bundle) {
  stressBundle(.out)
}
quit(save = "no", status = ifelse(any(.res$failed), 1L, 0L))
