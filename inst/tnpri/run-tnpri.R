#!/usr/bin/env Rscript
# babelmixr2 NONMEM $PRIOR TNPRI test kit
#
# Runs babelmixr2's TNPRI support end to end with NONMEM, plus a few
# hand-edited versions of the control stream it writes, and collects
# everything needed to check the results into one archive.  See
# README.md in this directory.
#
# Usage:
#   Rscript run-tnpri.R --nonmem=nmfe75 [options]
#
# Options:
#   --nonmem=COMMAND   command that runs NONMEM, like nmfe75 or
#                      /opt/nm75/run/nmfe75 (required unless --generate)
#   --generate         only write the control streams and data (no NONMEM)
#   --cases=REGEX      only the cases whose name matches REGEX
#   --out=DIR          output directory (default:
#                      babelmixr2-tnpri-YYYYMMDD-HHMMSS)
#   --list             list the cases and exit

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

# The models ---------------------------------------------------------------

one.cmt <- function() {
  ini({
    tka <- 0.45
    tcl <- 1.0
    tv <- 3.45
    add.sd <- 0.7
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl/v * central
    cp <- central / v
    cp ~ add(add.sd)
  })
}

# closed form (ADVAN2 TRANS1, no $MODEL, when rxode2 has linCmtMicro()) with
# an omega block
lin.cmt <- function() {
  ini({
    tka <- 0.45
    tcl <- 1.0
    tv <- 3.45
    prop.sd <- 0.1
    add.sd <- 0.7
    eta.ka ~ 0.6
    eta.cl + eta.v ~ c(0.3,
                       0.01, 0.1)
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    cp <- linCmt()
    cp ~ add(add.sd) + prop(prop.sd)
  })
}

# The prior study is the first half of the subjects, the new study the
# second half
.dat <- nlmixr2data::Oral_1CPT
.datA <- .dat[.dat$ID <= 60, ]
.datB <- .dat[.dat$ID > 60, ]

# The cases -----------------------------------------------------------------
#
# "babelmixr2" cases fit with nlmixr2(est="nonmem") and read the result
# back; "variant" cases take the TNPRI control stream babelmixr2 wrote for
# the `tnpri` case, change it, and run it directly with NONMEM, to learn
# what TNPRI accepts.

.cases <- list(
  list(name="reference", type="babelmixr2",
       description="one.cmt fit to the new data without a prior (the objective to compare to)",
       fit=function(runCommand, generate) {
         nlmixr(one.cmt, .datB, est="nonmem",
                control=nonmemControl(runCommand=runCommand, run=!generate,
                                      modelName="reference"))
       }),
  list(name="tnpri", type="babelmixr2",
       description="one.cmt (ODE): prior from a NONMEM fit of the prior data (dataset prior), then the TNPRI fit",
       fit=function(runCommand, generate) {
         nlmixr(one.cmt, .datB, est="nonmem",
                control=nonmemControl(runCommand=runCommand, run=!generate,
                                      modelName="tnpri",
                                      tnpri=nonmemTnpri(.datA)))
       }),
  list(name="tnpri-lincmt", type="babelmixr2",
       description="lin.cmt (closed form ADVAN2 TRANS1 when rxode2 has linCmtMicro(), otherwise ODEs; omega block; combined error): dataset prior with MODE=1",
       fit=function(runCommand, generate) {
         nlmixr(lin.cmt, .datB, est="nonmem",
                control=nonmemControl(runCommand=runCommand, run=!generate,
                                      modelName="tnprilin",
                                      tnpri=nonmemTnpri(.datA, mode=1)))
       }),
  list(name="tnpri-fit", type="babelmixr2",
       description="one.cmt: prior from an nlmixr2 focei fit (babelmixr2 refits its data with NONMEM)",
       fit=function(runCommand, generate) {
         .f <- nlmixr(one.cmt, .datA, est="focei",
                      control=foceiControl(print=0))
         nlmixr(one.cmt, .datB, est="nonmem",
                control=nonmemControl(runCommand=runCommand, run=!generate,
                                      modelName="tnprifit",
                                      tnpri=nonmemTnpri(.f)))
       }),
  list(name="tnpri-imp", type="refused",
       description="est='imp' with TNPRI is refused before NONMEM runs (NONMEM: do not use TNPRI with the NONMEM 7 methods)",
       fit=function(runCommand, generate) {
         nlmixr(one.cmt, .datB, est="nonmem",
                control=nonmemControl(runCommand=runCommand, run=!generate,
                                      modelName="tnpriimp", est="imp",
                                      tnpri=nonmemTnpri(.datA)))
       }),
  list(name="variant-no-plev", type="variant",
       description="as written by babelmixr2, without PLEV=0 (NONMEM's own default)",
       edit=function(p1, p2) {
         list(sub("(\\$PRIOR TNPRI \\(PROBLEM 2\\))[^\n]*", "\\1", p1), p2)
       }),
  list(name="variant-no-input2", type="variant",
       description="no $INPUT in problem 2 (is it needed?)",
       edit=function(p1, p2) {
         list(p1, .dropRecords(p2, "INPUT"))
       }),
  list(name="variant-code2", type="variant",
       description="the model code ($SUBROUTINES/$MODEL/$PK/$DES/$ERROR) repeated in problem 2 (is it allowed?)",
       edit=function(p1, p2) {
         .code <- .getRecords(p1, c("SUBROUTINES", "MODEL", "PK", "DES", "ERROR"))
         list(p1, .insertAfter(p2, "INPUT", .code))
       }),
  list(name="variant-no-code1", type="variant",
       description="no model code in problem 1, only in problem 2",
       edit=function(p1, p2) {
         .rec <- c("SUBROUTINES", "MODEL", "PK", "DES", "ERROR")
         .code <- .getRecords(p1, .rec)
         list(.dropRecords(p1, .rec), .insertAfter(p2, "INPUT", .code))
       }),
  list(name="variant-extra-theta", type="variant",
       description="problem 2 has one more $THETA than the prior (what does NONMEM do when they do not line up?)",
       edit=function(p1, p2) {
         list(p1, .insertAfter(p2, "THETA", "$THETA 0.5 ; unused extra THETA(5)\n\n"))
       })
)

# Control stream records ----------------------------------------------------

.splitRecords <- function(txt) {
  .l <- strsplit(txt, "\n", fixed=TRUE)[[1]]
  .start <- grepl("^\\$", .l)
  .grp <- cumsum(.start)
  .ret <- split(.l, .grp)
  vapply(.ret, function(x) paste0(paste(x, collapse="\n"), "\n"), character(1),
         USE.NAMES=FALSE)
}
.recName <- function(rec) {
  toupper(sub("^\\$([A-Za-z]+).*", "\\1", strsplit(rec, "\n", fixed=TRUE)[[1]][1]))
}
.getRecords <- function(txt, names) {
  .r <- .splitRecords(txt)
  paste(.r[vapply(.r, .recName, character(1)) %in% names], collapse="")
}
.dropRecords <- function(txt, names) {
  .r <- .splitRecords(txt)
  paste(.r[!(vapply(.r, .recName, character(1)) %in% names)], collapse="")
}
.insertAfter <- function(txt, name, what) {
  .r <- .splitRecords(txt)
  .w <- max(which(vapply(.r, .recName, character(1)) == name))
  paste(c(.r[seq_len(.w)], what, .r[-seq_len(.w)]), collapse="")
}
.splitProblems <- function(ctl) {
  .w <- gregexpr("(^|\n)\\$PROBLEM", ctl)[[1]]
  if (length(.w) != 2L) stop("expected a two problem control stream", call.=FALSE)
  .s <- .w[2] + 1L
  c(substr(ctl, 1L, .s - 1L), substr(ctl, .s, nchar(ctl)))
}

# Collecting the results ----------------------------------------------------

# value (or error) and the warnings of an expression
.try <- function(expr) {
  .w <- character(0)
  .r <- withCallingHandlers(
    tryCatch(list(value=expr), error=function(e) list(error=e)),
    warning=function(w) {
      .w <<- c(.w, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  .r$warnings <- .w
  .r
}

# Everything the NONMEM output of a directory says, the way nonmem2rx
# (and so babelmixr2) reads it
.readers <- function(dir, name) {
  .out <- c(paste0("# nonmem2rx readers for ", file.path(dir, name)), "")
  .one <- function(label, file, fun) {
    .f <- file.path(dir, file)
    if (!file.exists(.f)) return(c(paste0("## ", label, ": no ", file), ""))
    .r <- .try(withr::with_dir(dir, fun(file)))
    c(paste0("## ", label, " (", file, ")"),
      if (!is.null(.r$error)) paste0("ERROR: ", conditionMessage(.r$error))
      else utils::capture.output(utils::str(.r$value, max.level=2, vec.len=8, nchar.max=400)),
      if (length(.r$warnings)) paste0("WARNING: ", .r$warnings),
      "")
  }
  c(.out,
    .one("nminfo", paste0(name, ".lst"), nonmem2rx::nminfo),
    .one("nmext", paste0(name, ".ext"), nonmem2rx::nmext),
    .one("nmtab ext", paste0(name, ".ext"), nonmem2rx::nmtab),
    .one("nmcov", paste0(name, ".cov"), nonmem2rx::nmcov),
    .one("nmtab eta", paste0(name, ".eta"), nonmem2rx::nmtab))
}

# The key lines of a NONMEM listing
.lstSummary <- function(file) {
  if (!file.exists(file)) return(list(status="no output", objf=NA_real_, lines=character(0)))
  .l <- readLines(file, warn=FALSE)
  .key <- grep(paste0("PROBLEM NO|MINIMIZATION|TERMINATED|ERROR|WARNING|#OBJV|#TERM|",
                      "TNPRI|PRIOR|MSF|COVARIANCE STEP|NM-TRAN|AN ERROR|",
                      "NUMBER OF DATA RECORDS|TOT. NO. OF INDIVIDUALS|NONLINEAR MIXED"),
               .l, value=TRUE)
  .objv <- grep("#OBJV:", .l, value=TRUE)
  .objf <- if (length(.objv)) {
    suppressWarnings(as.numeric(gsub("[^0-9.Ee+-]", "", sub(".*#OBJV:", "", .objv[length(.objv)]))))
  } else NA_real_
  .status <- if (any(grepl("MINIMIZATION SUCCESSFUL", .l))) "successful"
             else if (any(grepl("MINIMIZATION TERMINATED", .l))) "terminated"
             else if (any(grepl("AN ERROR WAS FOUND|ERROR", .l))) "error"
             else "unknown"
  list(status=.status, objf=.objf, lines=utils::head(trimws(.key), 60))
}

.runNonmem <- function(cmd, dir, ctl, lst) {
  .full <- paste(cmd, ctl, lst)
  withr::with_dir(dir, system(.full, ignore.stdout=FALSE, ignore.stderr=FALSE))
}

# Running -------------------------------------------------------------------

.out <- .opt("out", paste0("babelmixr2-tnpri-", format(Sys.time(), "%Y%m%d-%H%M%S")))
.generate <- isTRUE(.opt("generate", FALSE))
.cmd <- .opt("nonmem", getOption("babelmixr2.nonmem", ""))
.re <- .opt("cases")
if (!is.null(.re)) {
  .keep <- grepl(.re, vapply(.cases, `[[`, character(1), "name"))
  # the variants need the tnpri case
  .keep <- .keep | (vapply(.cases, `[[`, character(1), "name") == "tnpri" &
                      any(.keep & vapply(.cases, `[[`, character(1), "type") == "variant"))
  .cases <- .cases[.keep]
}
if (isTRUE(.opt("list", FALSE))) {
  for (.c in .cases) cat(sprintf("%-22s %-10s %s\n", .c$name, .c$type, .c$description))
  quit(save="no", status=0)
}
if (!.generate && identical(.cmd, "")) {
  stop("give the NONMEM run command with --nonmem=COMMAND (like --nonmem=nmfe75), or use --generate",
       call.=FALSE)
}

dir.create(.out, showWarnings=FALSE, recursive=TRUE)
.out <- normalizePath(.out)
.runCommand <- if (.generate) NA else .cmd
.versions <- c(
  paste0("date: ", format(Sys.time())),
  paste0("R ", getRversion(), " on ", R.version$platform),
  paste0("babelmixr2 ", utils::packageVersion("babelmixr2"),
         ", nlmixr2est ", utils::packageVersion("nlmixr2est"),
         ", rxode2 ", utils::packageVersion("rxode2"),
         ", nonmem2rx ", utils::packageVersion("nonmem2rx")),
  paste0("NONMEM command: ", if (.generate) "(generate only)" else .cmd))
message(paste(.versions, collapse="\n"))
message("output: ", .out)

.res <- list()
.addRes <- function(case, status, objf=NA_real_, message="", seconds=NA_real_) {
  .res[[length(.res) + 1L]] <<- data.frame(case=case, status=status, objf=objf,
                                           seconds=seconds,
                                           message=substr(gsub("\n", " ", message), 1, 500))
}

.tnpriDir <- NULL
for (.c in .cases) {
  .dir <- file.path(.out, "cases", .c$name)
  dir.create(.dir, showWarnings=FALSE, recursive=TRUE)
  message("\n== ", .c$name, ": ", .c$description)
  .t0 <- proc.time()[["elapsed"]]
  if (.c$type == "refused") {
    # nothing may be run or written
    .try0 <- .try(withr::with_dir(.dir, suppressMessages(.c$fit(.runCommand, .generate))))
    .files <- list.files(.dir, recursive=TRUE)
    .ok <- !is.null(.try0$error) && length(.files) == 0L
    .msg <- if (is.null(.try0$error)) "not refused" else conditionMessage(.try0$error)
    writeLines(c(paste0("case: ", .c$name), .c$description, "",
                 paste0("status: ", if (.ok) "refused" else "problem"),
                 paste0("message: ", .msg),
                 paste0("files written: ", paste(.files, collapse=", "))),
               file.path(.dir, "case.txt"))
    .addRes(.c$name, if (.ok) "refused" else "problem", message=.msg,
            seconds=proc.time()[["elapsed"]] - .t0)
  } else if (.c$type == "babelmixr2") {
    .try0 <- .try(withr::with_dir(.dir, suppressMessages(.c$fit(.runCommand, .generate))))
    .r <- .try0$value
    .sec <- proc.time()[["elapsed"]] - .t0
    .log <- c(paste0("case: ", .c$name), .c$description, "")
    if (!is.null(.try0$error)) {
      .status <- "error"
      .msg <- conditionMessage(.try0$error)
      .objf <- NA_real_
    } else if (inherits(.r, "nlmixr2FitData")) {
      .status <- "ok"
      .msg <- ""
      .objf <- .r$objf
      .log <- c(.log, utils::capture.output(print(.r$objDf)), "",
                utils::capture.output(print(.r$parFixedDf)), "",
                utils::capture.output(print(.r$omega)), "",
                "$message:", .r$env$message)
    } else {
      .status <- if (.generate) "generated" else "not fit"
      .msg <- ""
      .objf <- NA_real_
    }
    if (length(.try0$warnings)) .log <- c(.log, "", paste0("WARNING: ", .try0$warnings))
    .log <- c(.log, "", paste0("status: ", .status), paste0("message: ", .msg))
    writeLines(.log, file.path(.dir, "case.txt"))
    # what nonmem2rx reads from every NONMEM run of the case
    for (.d in list.dirs(.dir, recursive=FALSE)) {
      .ctl <- list.files(.d, pattern="[.]nmctl$")
      for (.x in sub("[.]nmctl$", "", .ctl)) {
        writeLines(.readers(.d, .x), file.path(.d, paste0(.x, "-readers.txt")))
        .s <- .lstSummary(file.path(.d, paste0(.x, ".lst")))
        writeLines(c(paste0("status: ", .s$status), paste0("objf: ", .s$objf), "", .s$lines),
                   file.path(.d, paste0(.x, "-lst-summary.txt")))
      }
    }
    if (.c$name == "tnpri") .tnpriDir <- .dir
    .addRes(.c$name, .status, .objf, .msg, .sec)
  } else {
    # a hand-edited version of the tnpri case's control stream
    .src <- file.path(.tnpriDir, "tnpri-nonmem")
    .msf <- file.path(.tnpriDir, "tnpri_prior-nonmem", "tnpri_prior.msf")
    if (is.null(.tnpriDir) || !file.exists(file.path(.src, "tnpri.nmctl"))) {
      .addRes(.c$name, "skipped", message="the tnpri case did not write its control stream")
      next
    }
    .ctl <- paste(readLines(file.path(.src, "tnpri.nmctl")), collapse="\n")
    .p <- .c$edit(.splitProblems(.ctl)[1], .splitProblems(.ctl)[2])
    writeLines(paste0(.p[[1]], .p[[2]]), file.path(.dir, "tnpri.nmctl"))
    file.copy(file.path(.src, "tnpri.csv"), .dir, overwrite=TRUE)
    if (.generate) {
      .addRes(.c$name, "generated")
      next
    }
    if (!file.exists(.msf)) {
      .addRes(.c$name, "skipped", message="the tnpri case did not create the prior MSF")
      next
    }
    file.copy(.msf, .dir, overwrite=TRUE)
    .r <- .try(.runNonmem(.cmd, .dir, "tnpri.nmctl", "tnpri.lst"))
    .sec <- proc.time()[["elapsed"]] - .t0
    .s <- .lstSummary(file.path(.dir, "tnpri.lst"))
    writeLines(c(paste0("case: ", .c$name), .c$description, "",
                 paste0("status: ", .s$status), paste0("objf: ", .s$objf), "", .s$lines),
               file.path(.dir, "tnpri-lst-summary.txt"))
    writeLines(.readers(.dir, "tnpri"), file.path(.dir, "tnpri-readers.txt"))
    .addRes(.c$name, .s$status, .s$objf,
            if (!is.null(.r$error)) conditionMessage(.r$error) else "", .sec)
  }
}

.res <- do.call(rbind, .res)
utils::write.csv(.res, file.path(.out, "results.csv"), row.names=FALSE)

# summary.md: everything in text, in case only text can leave the machine
.md <- c("# babelmixr2 NONMEM TNPRI test kit", "", paste0("- ", .versions), "",
         "## Cases", "",
         "| case | status | objf | seconds | message |", "|---|---|---|---|---|",
         sprintf("| %s | %s | %s | %s | %s |", .res$case, .res$status,
                 ifelse(is.na(.res$objf), "", sprintf("%.4f", .res$objf)),
                 ifelse(is.na(.res$seconds), "", sprintf("%.0f", .res$seconds)),
                 gsub("\\|", "/", .res$message)))
for (.f in sort(list.files(file.path(.out, "cases"), pattern="(-lst-summary|-readers)[.]txt$|^case[.]txt$",
                           recursive=TRUE, full.names=TRUE))) {
  .md <- c(.md, "", paste0("## ", sub(paste0("^", .out, "/"), "", .f)), "", "```",
           readLines(.f, warn=FALSE), "```")
}
writeLines(.md, file.path(.out, "summary.md"))

# one archive with every control stream, data file and NONMEM output
.tar <- paste0(.out, ".tar.gz")
withr::with_dir(dirname(.out),
                utils::tar(basename(.tar), basename(.out), compression="gzip", tar="internal"))
message("\n", paste(.md[seq_len(min(length(.md), 12 + nrow(.res)))], collapse="\n"))
message("\nsummary: ", file.path(.out, "summary.md"))
message("archive: ", .tar)
