#' Stop explaining why a NONMEM run failed
#'
#' @param msg The message lines to show the user
#' @return Nothing, stops with the message
#' @author Matthew L. Fidler
#' @noRd
.nonmemFailureStop <- function(msg) {
  stop(paste(msg, collapse="\n"), call.=FALSE)
}

#' Read the lines NONMEM and NM-TRAN wrote about a run
#'
#' @param exportPath The NONMEM run directory
#' @param lst The NONMEM output file name
#' @return The lines of the output file followed by NM-TRAN's message
#'   file (`FMSG`), or `NULL` when neither exists
#' @author Matthew L. Fidler
#' @noRd
.nonmemFailureLines <- function(exportPath, lst) {
  .files <- file.path(exportPath, c(lst, "FMSG"))
  .files <- .files[file.exists(.files)]
  if (length(.files) == 0L) return(NULL)
  unlist(lapply(.files, function(f) {
    suppressWarnings(readLines(f, warn=FALSE))
  }), use.names=FALSE)
}

#' Pick the lines around the first match of a pattern
#'
#' @param lines The lines to search
#' @param pattern The regular expression to find
#' @param before The number of lines to keep before the match
#' @param after The number of lines to keep after the match
#' @return The non-empty lines around the first match, or
#'   `character(0)` without a match
#' @author Matthew L. Fidler
#' @noRd
.nonmemFailureContext <- function(lines, pattern, before=0L, after=4L) {
  .w <- grep(pattern, lines, ignore.case=TRUE)
  if (length(.w) == 0L) return(character(0))
  .w <- .w[1]
  .ret <- lines[seq(max(1L, .w - before), min(length(lines), .w + after))]
  .ret <- trimws(.ret, which="right")
  .ret[.ret != ""]
}

#' Classify what went wrong in a NONMEM run from its output
#'
#' @param lines The output lines (from `.nonmemFailureLines()`)
#' @return A list with the `cause` of the failure and the output
#'   `lines` that show it, or `NULL` when no failure is recognized
#' @author Matthew L. Fidler
#' @noRd
.nonmemClassifyFailure <- function(lines) {
  if (length(lines) == 0L) return(NULL)
  # a license about to expire is only a warning, so it is not a failure
  .lic <- grepl("licen[cs]e", lines, ignore.case=TRUE) &
    grepl("expired|not valid|invalid|not found|cannot find|could not find|missing|no valid|unable to|failed",
          lines, ignore.case=TRUE)
  if (any(.lic)) {
    return(list(cause="license",
                lines=.nonmemFailureContext(lines[which(.lic)[1]:length(lines)],
                                            ".", after=2L)))
  }
  if (any(grepl("(DATA ERROR)", lines, fixed=TRUE))) {
    return(list(cause="data",
                lines=.nonmemFailureContext(lines, "\\(DATA ERROR\\)")))
  }
  if (any(grepl("AN ERROR WAS FOUND IN THE CONTROL STATEMENTS", lines, fixed=TRUE))) {
    return(list(cause="controlStream",
                lines=.nonmemFailureContext(lines,
                                            "AN ERROR WAS FOUND IN THE CONTROL STATEMENTS",
                                            after=6L)))
  }
  if (any(grepl("PROGRAM TERMINATED", lines, fixed=TRUE))) {
    return(list(cause="crash",
                lines=.nonmemFailureContext(lines, "PROGRAM TERMINATED", after=6L)))
  }
  NULL
}

#' Show the last lines of NONMEM's output
#'
#' @param lines The output lines
#' @param n The number of non-empty lines to keep
#' @return The last `n` non-empty lines, indented
#' @author Matthew L. Fidler
#' @noRd
.nonmemFailureTail <- function(lines, n=10L) {
  .l <- trimws(lines, which="right")
  .l <- .l[.l != ""]
  paste0("  ", utils::tail(.l, n))
}

#' Explain a NONMEM failure and stop
#'
#' @param ui The rxode2 ui being run
#' @param status The exit status of the NONMEM run command, or `NULL`
#'   when it is not known (as with a user function)
#' @param readError The error from reading NONMEM's output, or `NULL`
#' @return Nothing, stops explaining why NONMEM failed when a cause is
#'   found (always when `readError` is given), otherwise returns `NULL`
#'   invisibly
#' @author Matthew L. Fidler
#' @noRd
.nonmemCheckRun <- function(ui, status=NULL, readError=NULL) {
  .exportPath <- ui$nonmemExportPath
  .lst <- ui$nonmemNmlst
  .lstFile <- file.path(.exportPath, .lst)
  .ctlFile <- file.path(.exportPath, ui$nonmemNmctl)
  .dataFile <- file.path(.exportPath, ui$nonmemCsv)
  .cmd <- rxode2::rxGetControl(ui, "runCommand", "")
  .cmdMsg <- if (is.character(.cmd)) {
    paste0("  run command: '", paste(.cmd, ui$nonmemNmctl, .lst), "' in '", .exportPath, "'")
  } else {
    paste0("  run command: the function given in nonmemControl(runCommand=) in '", .exportPath, "'")
  }
  .statusMsg <- if (is.null(status) || identical(as.integer(status), 0L)) {
    NULL
  } else {
    paste0("  exit status: ", status,
           if (identical(as.integer(status), 127L)) " (the command was not found)" else "")
  }
  .lines <- .nonmemFailureLines(.exportPath, .lst)
  if (!file.exists(.lstFile)) {
    .nonmemFailureStop(
      c(paste0("NONMEM did not create its output file '", .lstFile, "'"),
        .cmdMsg, .statusMsg,
        .nonmemFailureTail(.lines),
        "likely causes:",
        "  - nonmemControl(runCommand=) is not the right command or path to NONMEM (for example 'nmfe75'); check it runs from a terminal",
        "  - the command cannot run NONMEM on this system (ask your IT support)",
        "  - the NONMEM license is missing or expired (see the NONMEM messages printed above)",
        paste0("  - a runCommand function did not write '", .lst, "' in the run directory")))
  }
  .fail <- .nonmemClassifyFailure(.lines)
  if (!is.null(.fail)) {
    .msg <- switch(
      .fail$cause,
      license=c("NONMEM could not run because of a license problem:",
                paste0("  ", .fail$lines),
                "check NONMEM's license file (typically 'nonmem.lic' in the NONMEM 'license' directory)"),
      data=c("NM-TRAN found an error in the NONMEM data:",
             paste0("  ", .fail$lines),
             paste0("  data: '", .dataFile, "'"),
             paste0("  control stream: '", .ctlFile, "'"),
             "this is likely a problem in how babelmixr2 wrote the data; please report it at https://github.com/nlmixr2/babelmixr2/issues"),
      controlStream=c("NM-TRAN found an error in the NONMEM control stream:",
                      paste0("  ", .fail$lines),
                      paste0("  control stream: '", .ctlFile, "'"),
                      "this is likely a problem in how babelmixr2 translated the model; please report it at https://github.com/nlmixr2/babelmixr2/issues"),
      crash=c("NONMEM stopped during the run:",
              paste0("  ", .fail$lines),
              if (file.exists(file.path(.exportPath, "PRDERR"))) {
                paste0("  solving errors: '", file.path(.exportPath, "PRDERR"), "'")
              },
              paste0("  output: '", .lstFile, "'"),
              "changing the initial estimates or the model may help"))
    .nonmemFailureStop(.msg)
  }
  .started <- any(grepl("NONLINEAR MIXED EFFECTS MODEL PROGRAM", .lines, fixed=TRUE))
  .finished <- any(grepl("#TERM:|MINIMIZATION SUCCESSFUL|MINIMIZATION TERMINATED|OPTIMIZATION WAS COMPLETED|OPTIMIZATION WAS NOT COMPLETED|STOCHASTIC PORTION WAS|EXPECTATION ONLY PROCESS",
                           .lines))
  if (!.started) {
    .nonmemFailureStop(
      c(paste0("NONMEM did not start estimation (see '", .lstFile, "')"),
        .cmdMsg, .statusMsg,
        .nonmemFailureTail(.lines),
        "likely causes: the NONMEM license, the Fortran compiler, or NONMEM's installation; see the NONMEM messages printed above"))
  }
  if (!.finished) {
    .nonmemFailureStop(
      c(paste0("NONMEM started but did not finish (it may have crashed, run out of memory or been stopped); see '", .lstFile, "'"),
        .statusMsg,
        "  last output:",
        .nonmemFailureTail(.lines)))
  }
  if (!is.null(readError)) {
    .nonmemFailureStop(
      c(paste0("babelmixr2 could not read NONMEM's output '", .lstFile, "':"),
        paste0("  ", conditionMessage(readError)),
        .nonmemFailureTail(.lines),
        "if NONMEM finished, please report this at https://github.com/nlmixr2/babelmixr2/issues"))
  }
  invisible(NULL)
}
