#' Stop explaining why a NONMEM run failed
#'
#' @param msg The message lines to show the user
#' @return Nothing, stops with the message
#' @author Matthew L. Fidler
#' @noRd
.nonmemFailureStop <- function(msg) {
  stop(paste(msg, collapse = "\n"), call. = FALSE)
}

#' Read a NONMEM file's lines as valid UTF-8
#'
#' @param file The file to read
#' @return The lines, with bytes that are not UTF-8 (like a Latin-1
#'   path) written as `<xx>`
#' @author Matthew L. Fidler
#' @noRd
.nonmemFailureReadLines <- function(file) {
  iconv(
    suppressWarnings(readLines(file, warn = FALSE)),
    "UTF-8",
    "UTF-8",
    sub = "byte"
  )
}

#' Read the lines NONMEM and NM-TRAN wrote about a run
#'
#' NONMEM's output starts with a copy of the control stream; those
#' lines are dropped so the model's own code is never mistaken for a
#' NONMEM message.  Control stream lines NM-TRAN shows with an error
#' are kept.
#'
#' @param exportPath The NONMEM run directory
#' @param lst The NONMEM output file name
#' @param ctl The NONMEM control stream file name
#' @return The lines of the output file followed by NM-TRAN's message
#'   file (`FMSG`), without the lines of the control stream, or `NULL`
#'   when neither exists
#' @author Matthew L. Fidler
#' @noRd
.nonmemFailureLines <- function(exportPath, lst, ctl) {
  .files <- file.path(exportPath, c(lst, "FMSG"))
  .files <- .files[file.exists(.files)]
  if (length(.files) == 0L) {
    return(NULL)
  }
  .ret <- lapply(.files, .nonmemFailureReadLines)
  .ctlFile <- file.path(exportPath, ctl)
  if (file.exists(file.path(exportPath, lst)) && file.exists(.ctlFile)) {
    .ctl <- trimws(.nonmemFailureReadLines(.ctlFile))
    .ctl <- .ctl[.ctl != ""]
    .lst <- .ret[[1]]
    # the copy is before NM-TRAN's messages and NONMEM's banner; later
    # copies of control stream lines are NM-TRAN showing an error
    .end <- grep(
      "NM-TRAN MESSAGES|WARNINGS AND ERRORS|NONLINEAR MIXED EFFECTS MODEL PROGRAM",
      .lst
    )
    .end <- if (length(.end) == 0L) length(.lst) else .end[1] - 1L
    .echo <- seq_along(.lst) <= .end & trimws(.lst) %in% .ctl
    .ret[[1]] <- .lst[!.echo]
  }
  unlist(.ret, use.names = FALSE)
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
.nonmemFailureContext <- function(lines, pattern, before = 0L, after = 4L) {
  .w <- grep(pattern, lines, ignore.case = TRUE)
  if (length(.w) == 0L) {
    return(character(0))
  }
  .w <- .w[1]
  .ret <- lines[seq(max(1L, .w - before), min(length(lines), .w + after))]
  .ret <- trimws(.ret, which = "right")
  .ret[.ret != ""]
}

#' Drop the model name from NONMEM's output lines
#'
#' The model name is in NONMEM's problem title and file names, so it
#' is removed (as a whole name) before the output is searched for
#' messages.
#'
#' @param lines The output lines
#' @param modelName The model name, or `NULL`
#' @return The lines without the model name
#' @author Matthew L. Fidler
#' @noRd
.nonmemDropModelName <- function(lines, modelName) {
  if (is.null(modelName) || modelName == "" || length(lines) == 0L) {
    return(lines)
  }
  .name <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", modelName)
  gsub(
    paste0("(?<![[:alnum:]_.])", .name, "(?![[:alnum:]_])"),
    "",
    lines,
    perl = TRUE
  )
}

#' Classify what went wrong in a NONMEM run from its output
#'
#' @param lines The output lines (from `.nonmemFailureLines()`)
#' @param modelName The model name, which is never read as a license
#'   message (it is in NONMEM's problem title and file names)
#' @return A list with the `cause` of the failure and the output
#'   `lines` that show it, or `NULL` when no failure is recognized
#' @author Matthew L. Fidler
#' @noRd
.nonmemClassifyFailure <- function(lines, modelName = NULL) {
  if (length(lines) == 0L) {
    return(NULL)
  }
  if (any(grepl("(DATA ERROR)", lines, fixed = TRUE))) {
    return(list(
      cause = "data",
      lines = .nonmemFailureContext(lines, "\\(DATA ERROR\\)")
    ))
  }
  if (any(grepl("DATA FILE DOES NOT EXIST", lines, fixed = TRUE))) {
    return(list(
      cause = "data",
      lines = .nonmemFailureContext(
        lines,
        "DATA FILE DOES NOT EXIST",
        before = 3L,
        after = 0L
      )
    ))
  }
  if (
    any(grepl(
      "AN ERROR WAS FOUND IN THE CONTROL STATEMENTS",
      lines,
      fixed = TRUE
    ))
  ) {
    return(list(
      cause = "controlStream",
      lines = .nonmemFailureContext(
        lines,
        "AN ERROR WAS FOUND IN THE CONTROL STATEMENTS",
        after = 6L
      )
    ))
  }
  .term <- grep("PROGRAM TERMINATED", lines, fixed = TRUE)
  .tere <- grep("#TERE:", lines, fixed = TRUE)
  if (length(.tere) > 0L) {
    # after the estimation's termination block, NONMEM stopping is a
    # failure of the covariance or table step; estimation finished
    .term <- .term[.term < .tere[1]]
  }
  if (length(.term) > 0L) {
    .start <- grep("NONLINEAR MIXED EFFECTS MODEL PROGRAM", lines, fixed = TRUE)
    # before NONMEM starts, it is NM-TRAN that stopped
    .cause <- if (length(.start) > 0L && .start[1] < .term[1]) {
      "crash"
    } else {
      "nmtran"
    }
    return(list(
      cause = .cause,
      lines = .nonmemFailureContext(
        lines,
        "PROGRAM TERMINATED",
        before = 2L,
        after = 6L
      )
    ))
  }
  # checked last, so a license warning never hides another failure;
  # the registration line and warnings (like a license about to
  # expire) are not failures
  # directories in file paths (like a compiler's) are not license
  # messages; a license file (.lic) is
  lines <- gsub(
    "[^[:space:]]*[/\\\\]",
    "",
    .nonmemDropModelName(lines, modelName)
  )
  .lic <- grepl("licen[cs]e|[.]lic\\b", lines, ignore.case = TRUE) &
    grepl(
      "expired|not valid|invalid|not found|cannot find|could not find|missing|no valid|unable to|failed|error|terminating",
      lines,
      ignore.case = TRUE
    ) &
    !grepl("registered to|warning", lines, ignore.case = TRUE)

  if (any(.lic)) {
    return(list(
      cause = "license",
      lines = .nonmemFailureContext(
        lines[which(.lic)[1]:length(lines)],
        ".",
        after = 2L
      )
    ))
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
.nonmemFailureTail <- function(lines, n = 10L) {
  .l <- trimws(lines, which = "right")
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
.nonmemCheckRun <- function(ui, status = NULL, readError = NULL) {
  .exportPath <- ui$nonmemExportPath
  .lst <- ui$nonmemNmlst
  .lstFile <- file.path(.exportPath, .lst)
  .ctlFile <- file.path(.exportPath, ui$nonmemNmctl)
  .dataFile <- file.path(.exportPath, ui$nonmemCsv)
  .cmd <- rxode2::rxGetControl(ui, "runCommand", "")
  .cmdMsg <- if (is.character(.cmd)) {
    paste0(
      "  run command: '",
      paste(.cmd, ui$nonmemNmctl, .lst),
      "' in '",
      .exportPath,
      "'"
    )
  } else {
    paste0(
      "  run command: the function given in nonmemControl(runCommand=) in '",
      .exportPath,
      "'"
    )
  }
  .statusMsg <- if (is.null(status) || identical(as.integer(status), 0L)) {
    NULL
  } else {
    paste0(
      "  exit status: ",
      status,
      if (identical(as.integer(status), 127L)) {
        " (the command was not found)"
      } else {
        ""
      }
    )
  }
  .lines <- .nonmemFailureLines(.exportPath, .lst, ui$nonmemNmctl)
  .fail <- .nonmemClassifyFailure(.lines, ui$nonmemModelName)
  if (!is.null(.fail)) {
    .msg <- switch(
      .fail$cause,
      license = c(
        "NONMEM could not run because of a license problem:",
        paste0("  ", .fail$lines),
        "check NONMEM's license file (typically 'nonmem.lic' in the NONMEM 'license' directory)"
      ),
      data = c(
        "NM-TRAN found an error in the NONMEM data:",
        paste0("  ", .fail$lines),
        paste0("  data: '", .dataFile, "'"),
        paste0("  control stream: '", .ctlFile, "'"),
        "this is likely a problem in how babelmixr2 wrote the data; please report it at https://github.com/nlmixr2/babelmixr2/issues"
      ),
      controlStream = c(
        "NM-TRAN found an error in the NONMEM control stream:",
        paste0("  ", .fail$lines),
        paste0("  control stream: '", .ctlFile, "'"),
        "this is likely a problem in how babelmixr2 translated the model; please report it at https://github.com/nlmixr2/babelmixr2/issues"
      ),
      nmtran = c(
        "NM-TRAN stopped before NONMEM could run:",
        paste0("  ", .fail$lines),
        paste0("  control stream: '", .ctlFile, "'"),
        paste0("  data: '", .dataFile, "'"),
        "this is likely a problem in how babelmixr2 translated the model or data; please report it at https://github.com/nlmixr2/babelmixr2/issues"
      ),
      crash = c(
        "NONMEM stopped during the run:",
        paste0("  ", .fail$lines),
        if (file.exists(file.path(.exportPath, "PRDERR"))) {
          paste0("  solving errors: '", file.path(.exportPath, "PRDERR"), "'")
        },
        paste0("  output: '", .lstFile, "'"),
        "changing the initial estimates or the model may help"
      )
    )
    .nonmemFailureStop(.msg)
  }
  if (!file.exists(.lstFile)) {
    .nonmemFailureStop(
      c(
        paste0("NONMEM did not create its output file '", .lstFile, "'"),
        .cmdMsg,
        .statusMsg,
        .nonmemFailureTail(.lines),
        "likely causes:",
        "  - nonmemControl(runCommand=) is not the right command or path to NONMEM (for example 'nmfe75'); check it runs from a terminal",
        "  - the command cannot run NONMEM on this system (ask your IT support)",
        "  - the NONMEM license is missing or expired (see the NONMEM messages printed above)",
        paste0(
          "  - a runCommand function did not write '",
          .lst,
          "' in the run directory"
        )
      )
    )
  }
  .started <- any(grepl(
    "NONLINEAR MIXED EFFECTS MODEL PROGRAM",
    .lines,
    fixed = TRUE
  ))
  # an evaluation (MAXEVALS=0) or a model with every parameter fixed
  # has no termination message
  .finished <- any(grepl(
    "#TERM:|MINIMIZATION SUCCESSFUL|MINIMIZATION TERMINATED|OPTIMIZATION WAS COMPLETED|OPTIMIZATION WAS NOT COMPLETED|STOCHASTIC PORTION WAS|EXPECTATION ONLY PROCESS|ESTIMATION STEP OMITTED: +YES|NUMBER OF PARAMETERS TO BE ESTIMATED IS 0",
    .nonmemDropModelName(.lines, ui$nonmemModelName)
  ))
  if (!.started) {
    .nonmemFailureStop(
      c(
        paste0("NONMEM did not start estimation (see '", .lstFile, "')"),
        .cmdMsg,
        .statusMsg,
        .nonmemFailureTail(.lines),
        "likely causes: the NONMEM license, the Fortran compiler, or NONMEM's installation; see the NONMEM messages printed above"
      )
    )
  }
  if (!.finished) {
    .nonmemFailureStop(
      c(
        paste0(
          "NONMEM started but did not finish (it may have crashed, run out of memory or been stopped); see '",
          .lstFile,
          "'"
        ),
        .statusMsg,
        "  last output:",
        .nonmemFailureTail(.lines)
      )
    )
  }
  if (!is.null(readError)) {
    if (!is.null(.statusMsg)) {
      # NONMEM finished estimating but then exited abnormally (for
      # example killed while computing the covariance or tables)
      .nonmemFailureStop(
        c(
          paste0(
            "NONMEM exited with an error after estimation and babelmixr2 could not use its output; see '",
            .lstFile,
            "'"
          ),
          .statusMsg,
          paste0("  error: ", conditionMessage(readError)),
          .nonmemFailureTail(.lines),
          "NONMEM may have crashed, run out of memory or been stopped while writing its output"
        )
      )
    }
    .nonmemFailureStop(
      c(
        paste0(
          "babelmixr2 could not read or use NONMEM's output '",
          .lstFile,
          "':"
        ),
        paste0("  ", conditionMessage(readError)),
        .nonmemFailureTail(.lines),
        "if NONMEM finished, please report this at https://github.com/nlmixr2/babelmixr2/issues"
      )
    )
  }
  invisible(NULL)
}

#' Finalize a NONMEM fit, explaining any failure to read its output
#'
#' @param ret The nlmixr2 fit environment
#' @param ui The rxode2 ui being run
#' @param status The exit status of the NONMEM run command, or `NULL`
#' @return The finalized fit; stops explaining why when NONMEM's
#'   output cannot be read
#' @author Matthew L. Fidler
#' @noRd
.nonmemFinalizeOrExplain <- function(ret, ui, status = NULL) {
  tryCatch(.nonmemFinalizeEnv(ret, ui), error = function(e) {
    .nonmemCheckRun(ui, status, readError = e)
  })
}

#' Remove the output of an earlier NONMEM run before running again
#'
#' An earlier run that did not finish leaves its output behind; it is
#' removed so a failure of the new run is never explained with the
#' old run's output.
#'
#' @param ui The rxode2 ui being run
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nonmemRemoveOldOutput <- function(ui) {
  .exportPath <- ui$nonmemExportPath
  unlink(file.path(.exportPath, c(ui$nonmemNmlst, "FMSG", "PRDERR")))
  invisible()
}
