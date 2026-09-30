# babelmixr2 NONMEM/Monolix stress test: finding NONMEM/Monolix and the report
#
# Sourced by stress.R; see README.md in this directory.

#' The command that runs NONMEM
#'
#' @return the babelmixr2.nonmem option, an nmfe7* on the PATH or in
#'   the usual install directories, or "" when none is found
#' @noRd
stressFindNonmem <- function() {
  .o <- getOption("babelmixr2.nonmem", "")
  if (is.character(.o) && nzchar(.o)) {
    return(.o)
  }
  .w <- Sys.which(paste0("nmfe7", 9:0))
  .w <- .w[.w != ""]
  if (length(.w) > 0L) {
    return(unname(.w[1]))
  }
  .globs <- c(
    "/opt/NONMEM/*/run/nmfe7*",
    "/opt/nm*/run/nmfe7*",
    "/usr/local/NONMEM/*/run/nmfe7*",
    "/usr/local/nm*/run/nmfe7*",
    "~/nm*/run/nmfe7*",
    "~/NONMEM/*/run/nmfe7*",
    "C:/nm*/run/nmfe7*.bat",
    "C:/NONMEM/*/run/nmfe7*.bat"
  )
  .f <- Sys.glob(path.expand(.globs))
  .f <- .f[!grepl("\\.(f90|o|obj)$", .f) & file.exists(.f)]
  if (length(.f) == 0L) {
    return("")
  }
  # newest NONMEM first
  sort(.f, decreasing = TRUE)[1]
}

#' How Monolix is run
#'
#' @return description of the Monolix run command or lixoftConnectors,
#'   or "" when Monolix is not found
#' @noRd
stressMonolixStatus <- function() {
  .o <- getOption("babelmixr2.monolix", "")
  if (is.character(.o) && nzchar(.o)) {
    return(paste0("command ", .o))
  }
  if (!requireNamespace("lixoftConnectors", quietly = TRUE)) {
    return("")
  }
  .x <- try(
    suppressMessages(
      lixoftConnectors::initializeLixoftConnectors(
        software = "monolix",
        force = TRUE
      )
    ),
    silent = TRUE
  )
  if (inherits(.x, "try-error") || isFALSE(.x)) {
    return("")
  }
  paste0("lixoftConnectors ", utils::packageVersion("lixoftConnectors"))
}

#' Package versions for the report
#'
#' @return markdown list lines
#' @noRd
stressVersions <- function() {
  .p <- c(
    "babelmixr2",
    "rxode2",
    "nlmixr2est",
    "lotri",
    "nonmem2rx",
    "monolix2rx",
    "nlmixr2lib",
    "lixoftConnectors"
  )
  .v <- vapply(
    .p,
    function(p) {
      if (requireNamespace(p, quietly = TRUE)) {
        as.character(utils::packageVersion(p))
      } else {
        "-"
      }
    },
    character(1)
  )
  .sha <- utils::packageDescription("babelmixr2")$RemoteSha
  c(
    paste0(
      "- ",
      .p,
      " ",
      .v,
      ifelse(
        .p == "babelmixr2" & !is.null(.sha),
        paste0(" (", substr(.sha, 1, 7), ")"),
        ""
      )
    ),
    paste0("- R ", getRversion(), " on ", R.version$platform)
  )
}


#' Markdown summary of a stress test
#'
#' @param res results (from `stressRun()`, with a `failed` column)
#' @param modes modes that were run
#' @param engines engines that were used
#' @param nonmem,monolix how NONMEM/Monolix were run
#' @return markdown lines
#' @noRd
stressSummary <- function(res, modes, engines, nonmem, monolix) {
  .md <- c(
    "# babelmixr2 NONMEM/Monolix stress test",
    "",
    paste0("- date: ", format(Sys.time())),
    paste0("- mode: ", paste(modes, collapse = ", ")),
    stressVersions(),
    if ("run" %in% modes && "nonmem" %in% engines) {
      paste0("- NONMEM: ", nonmem)
    },
    if ("run" %in% modes && "monolix" %in% engines) {
      paste0("- Monolix: ", monolix)
    },
    "",
    "## Summary",
    ""
  )
  .col <- paste(res$mode, res$engine)
  .tab <- table(res$status, .col)
  .md <- c(
    .md,
    paste0("| status | ", paste(colnames(.tab), collapse = " | "), " |"),
    paste0("|---|", paste(rep("---", ncol(.tab)), collapse = "|"), "|"),
    vapply(
      rownames(.tab),
      function(r) {
        paste0("| ", r, " | ", paste(.tab[r, ], collapse = " | "), " |")
      },
      character(1)
    ),
    "",
    "## Failures",
    ""
  )
  .fail <- res[res$failed, ]
  if (nrow(.fail) == 0L) {
    .md <- c(.md, "None.")
  } else {
    .md <- c(
      .md,
      "| case | engine | mode | status | message |",
      "|---|---|---|---|---|",
      sprintf(
        "| %s | %s | %s | %s | %s |",
        .fail$case,
        .fail$engine,
        .fail$mode,
        .fail$status,
        gsub("\\|", "/", substr(paste(.fail$message, .fail$problems), 1, 300))
      )
    )
  }
  if ("run" %in% modes) {
    .ok <- res[res$mode == "run" & res$status %in% c("ok", "problem"), ]
    .num <- function(x, fmt) ifelse(is.na(x), "", sprintf(fmt, x))
    .md <- c(
      .md,
      "",
      "## Fits",
      "",
      paste(
        "| case | engine | status | seconds | objective | IPRED diff % |",
        "PRED diff % | rerun seconds | max rel. diff vs nlmixr2 |"
      ),
      "|---|---|---|---|---|---|---|---|---|",
      sprintf(
        "| %s | %s | %s | %.1f | %s | %s | %s | %s | %s |",
        .ok$case,
        .ok$engine,
        .ok$status,
        .ok$seconds,
        .num(.ok$objf, "%.3f"),
        .num(.ok$ipredRelDiff, "%.3f"),
        .num(.ok$predRelDiff, "%.3f"),
        .num(.ok$rerunSeconds, "%.1f"),
        .num(.ok$maxRelDiffTheta, "%.3f")
      )
    )
  }
  .md
}

#' Zip the output directory to send back
#'
#' @param out output directory
#' @return the zip file (invisibly), or NULL when it cannot be made
#' @noRd
stressBundle <- function(out) {
  .zip <- paste0(out, ".zip")
  .ok <- withr::with_dir(dirname(out), {
    try(utils::zip(.zip, basename(out), flags = "-r9Xq"), silent = TRUE)
  })
  if (inherits(.ok, "try-error") || !file.exists(.zip)) {
    message("could not zip the output; send the directory ", out, " instead")
    return(invisible(NULL))
  }
  message("send this file back: ", .zip)
  invisible(.zip)
}
