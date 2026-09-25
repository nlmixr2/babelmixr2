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
#   --engine=nonmem,monolix   engines to test (default: both)
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
#                             (default: option babelmixr2.nonmem)
#   --monolix=COMMAND         command that runs Monolix (default: option
#                             babelmixr2.monolix, or lixoftConnectors
#                             when it is installed)
#   --reference               in run mode, also fit each model with
#                             nlmixr2 (focei for NONMEM, saem for
#                             Monolix) and compare the estimates
#   --list                    list the cases and exit

.args <- commandArgs(trailingOnly=TRUE)
.opt <- function(name, default=NULL) {
  .w <- grep(paste0("^--", name, "(=|$)"), .args, value=TRUE)
  if (length(.w) == 0L) return(default)
  if (!grepl("=", .w[1])) return(TRUE)
  sub(paste0("^--", name, "="), "", .w[1])
}

suppressPackageStartupMessages({
  library(babelmixr2)
  library(nlmixr2est)
  library(rxode2)
})

.file <- system.file("stress", "stress.R", package="babelmixr2")
if (.file == "") {
  # running from a source checkout
  .file <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])),
                     "stress.R")
}
source(.file)

.engines <- strsplit(.opt("engine", "nonmem,monolix"), ",")[[1]]
.engines <- match.arg(.engines, c("nonmem", "monolix"), several.ok=TRUE)
.mode <- match.arg(.opt("mode", "translate"), c("translate", "run"))
.lib <- match.arg(.opt("nlmixr2lib", "none"), c("none", "sample", "all"))
.out <- .opt("out", paste0("babelmixr2-stress-", format(Sys.time(), "%Y%m%d-%H%M%S")))
.reference <- isTRUE(.opt("reference", FALSE))

.cases <- stressCases()
if (.lib != "none") {
  if (!requireNamespace("nlmixr2lib", quietly=TRUE)) {
    stop("--nlmixr2lib needs the nlmixr2lib package", call.=FALSE)
  }
  .cases <- c(.cases, stressLibCases(all=(.lib == "all")))
}
.re <- .opt("cases")
if (!is.null(.re)) {
  .cases <- .cases[grepl(.re, names(.cases))]
}
if (isTRUE(.opt("list", FALSE))) {
  for (.c in .cases) {
    cat(sprintf("%-50s nonmem=%-8s monolix=%-8s %s\n", .c$name,
                ifelse(.c$expect$nonmem == "ok", "ok", "refuse"),
                ifelse(.c$expect$monolix == "ok", "ok", "refuse"),
                .c$description))
  }
  quit(save="no", status=0)
}

.runCommand <- list(nonmem=.opt("nonmem"), monolix=.opt("monolix"))
if (.mode == "run") {
  if ("nonmem" %in% .engines && is.null(.runCommand$nonmem) &&
        identical(getOption("babelmixr2.nonmem", ""), "")) {
    stop("run mode needs --nonmem=COMMAND (like --nonmem=nmfe75) or options(babelmixr2.nonmem=)",
         call.=FALSE)
  }
  if ("monolix" %in% .engines && is.null(.runCommand$monolix) &&
        identical(getOption("babelmixr2.monolix", ""), "") &&
        !requireNamespace("lixoftConnectors", quietly=TRUE)) {
    stop("run mode needs lixoftConnectors, --monolix=COMMAND or options(babelmixr2.monolix=)",
         call.=FALSE)
  }
}

dir.create(.out, showWarnings=FALSE, recursive=TRUE)
message("babelmixr2 ", utils::packageVersion("babelmixr2"),
        ", rxode2 ", utils::packageVersion("rxode2"),
        ", nlmixr2est ", utils::packageVersion("nlmixr2est"))
message(length(.cases), " cases; engines: ", paste(.engines, collapse=", "),
        "; mode: ", .mode, "; output: ", normalizePath(.out))

.res <- stressRun(.cases, engines=.engines, mode=.mode, dir=.out,
                  runCommand=.runCommand, reference=.reference, progress=TRUE)
.res$failed <- stressFailed(.res)
utils::write.csv(.res, file.path(.out, "results.csv"), row.names=FALSE)

# markdown summary
.tab <- table(.res$status, .res$engine)
.md <- c("# babelmixr2 NONMEM/Monolix stress test", "",
         paste0("- date: ", format(Sys.time())),
         paste0("- mode: ", .mode),
         paste0("- babelmixr2 ", utils::packageVersion("babelmixr2"),
                ", rxode2 ", utils::packageVersion("rxode2"),
                ", nlmixr2est ", utils::packageVersion("nlmixr2est")),
         paste0("- R ", getRversion(), " on ", R.version$platform),
         "", "## Summary", "",
         paste0("| status | ", paste(colnames(.tab), collapse=" | "), " |"),
         paste0("|---|", paste(rep("---", ncol(.tab)), collapse="|"), "|"),
         vapply(rownames(.tab), function(r) {
           paste0("| ", r, " | ", paste(.tab[r, ], collapse=" | "), " |")
         }, character(1)),
         "", "## Failures", "")
.fail <- .res[.res$failed, ]
if (nrow(.fail) == 0L) {
  .md <- c(.md, "None.")
} else {
  .md <- c(.md, "| case | engine | status | message |", "|---|---|---|---|",
           sprintf("| %s | %s | %s | %s |", .fail$case, .fail$engine, .fail$status,
                   gsub("\\|", "/", substr(paste(.fail$message, .fail$problems), 1, 300))))
}
if (.mode == "run") {
  .ok <- .res[.res$status == "ok", ]
  .md <- c(.md, "", "## Fits", "",
           "| case | engine | seconds | objective | max rel. diff vs nlmixr2 |",
           "|---|---|---|---|---|",
           sprintf("| %s | %s | %.1f | %.3f | %s |", .ok$case, .ok$engine, .ok$seconds,
                   .ok$objf, ifelse(is.na(.ok$maxRelDiffTheta), "",
                                    sprintf("%.3f", .ok$maxRelDiffTheta))))
}
writeLines(.md, file.path(.out, "summary.md"))
message("\n", paste(.md, collapse="\n"))
message("\nresults: ", normalizePath(file.path(.out, "results.csv")))
quit(save="no", status=ifelse(any(.res$failed), 1L, 0L))
