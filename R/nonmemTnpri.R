## NONMEM `$PRIOR TNPRI` from a previous NONMEM run
##
## TNPRI does not take its prior from the control stream.  NONMEM reads
## it from the model specification file (MSF) of an earlier run of the
## same model: the prior means are that run's estimates and the prior
## variances are its covariance matrix of the estimates.  So a TNPRI
## prior is always "the model fit to other data", the way it is used in
## practice, and never a distribution written in `ini({})` (those are
## NWPRI, see `nonmemRxUiGetPrior.R`).
##
## The control stream then has two problems:
##
## - problem 1 reads the MSF (`$MSFI ... ONLYREAD`), says which problem
##   the prior is for (`$PRIOR TNPRI (PROBLEM 2)`) and holds the model's
##   abbreviated code, which NONMEM uses for every problem of the run
##
## - problem 2 is the estimation with the new data (`$DATA ... REWIND`),
##   the initial estimates, `$ESTIMATION`, `$COVARIANCE` and `$TABLE`
##
## Since both problems share the model code, the prior run has to be a
## fit of the same model; TNPRI then lines the parameters up by
## position, which is only right when the two models are the same.

#' NONMEM TNPRI prior from a previous NONMEM run
#'
#' Use the estimates and covariance of an earlier NONMEM fit of the
#' same model as a prior for a new fit (NONMEM's `$PRIOR TNPRI`), by
#' giving this to `nonmemControl(tnpri=)`.
#'
#' NONMEM reads a TNPRI prior from the model specification file (MSF)
#' of the earlier run, so that run needs NONMEM's covariance step.
#' babelmixr2 runs it for you when `prior` is a dataset or an nlmixr2
#' fit: the model is fit to that data with NONMEM (with `$MSFO` and
#' `$COVARIANCE`) in its own directory, `<modelName>_prior-nonmem`, and
#' the new fit reads the resulting MSF.  Like every babelmixr2 NONMEM
#' run it is cached, so it only runs once.
#'
#' The control stream of the new fit then has two problems: the first
#' reads the MSF and holds the model code, and the second estimates the
#' model with the new data and the prior.
#'
#' @param prior Where the prior comes from:
#'
#' - a `data.frame`: the data of the earlier study; babelmixr2 fits the
#'   model to it with NONMEM to create the prior
#'
#' - an nlmixr2 fit of the same model (from any estimation method): its
#'   data are fit with NONMEM to create the prior
#'
#' - a character string: the path to the MSF of a NONMEM run you made
#'   yourself (for example with `nonmemControl(msfo=TRUE)`).  It has to
#'   be a run of the same model, with the same `THETA`s, `OMEGA`s and
#'   `SIGMA`s in the same order, and with the covariance step; babelmixr2
#'   cannot check this.
#'
#' @param plev NONMEM's `PLEV` option of `$PRIOR TNPRI`; `NULL` uses
#'   NONMEM's default
#' @param ivar NONMEM's `IVAR` option of `$PRIOR TNPRI` (for example `1`
#'   to use the variance-covariance matrix of the prior problem); `NULL`
#'   uses NONMEM's default
#' @param ityp NONMEM's `ITYP` option; `NULL` uses NONMEM's default
#' @param ifnd NONMEM's `IFND` option; `NULL` uses NONMEM's default
#' @param mode NONMEM's `MODE` option; `NULL` uses NONMEM's default
#' @param display add NONMEM's `DISPLAY` option
#'
#' @return a `nonmemTnpri` object for `nonmemControl(tnpri=)`
#'
#' @details
#'
#' The prior and the new fit have to use the same NONMEM version.
#'
#' A TNPRI prior cannot be combined with the priors written in the
#' `ini({})` block (those become `$PRIOR NWPRI`); a NONMEM problem has
#' only one `$PRIOR`.
#'
#' NONMEM's objective function includes the prior, so the fit's
#' objective type ends in `tnpri` and is not compared with fits without
#' it.
#'
#' @author Matthew L. Fidler
#' @export
#' @examples
#'
#' nonmemTnpri("prior-run/prior.msf", plev=0.999)
#'
nonmemTnpri <- function(prior, plev=NULL, ivar=NULL, ityp=NULL, ifnd=NULL,
                        mode=NULL, display=FALSE) {
  if (inherits(prior, "nonmemTnpri")) return(prior)
  if (is.character(prior)) {
    checkmate::assertCharacter(prior, len=1, any.missing=FALSE, min.chars=1)
    .type <- "msf"
  } else if (inherits(prior, "nlmixr2FitCore")) {
    .type <- "fit"
  } else if (is.data.frame(prior)) {
    .type <- "data"
  } else {
    stop("'prior' has to be a dataset, an nlmixr2 fit or the path to a NONMEM model specification file (MSF)",
         call.=FALSE)
  }
  if (!is.null(plev)) {
    checkmate::assertNumber(plev, lower=0, upper=1, finite=TRUE)
    if (plev <= 0 || plev >= 1) {
      stop("'plev' has to be between 0 and 1", call.=FALSE)
    }
  }
  for (.n in c("ivar", "ityp", "ifnd", "mode")) {
    .v <- get(.n)
    if (!is.null(.v)) {
      checkmate::assertIntegerish(.v, lower=0, len=1, any.missing=FALSE,
                                  .var.name=.n)
    }
  }
  checkmate::assertLogical(display, len=1, any.missing=FALSE)
  .ret <- list(prior=prior, type=.type, plev=plev,
               ivar=if (is.null(ivar)) NULL else as.integer(ivar),
               ityp=if (is.null(ityp)) NULL else as.integer(ityp),
               ifnd=if (is.null(ifnd)) NULL else as.integer(ifnd),
               mode=if (is.null(mode)) NULL else as.integer(mode),
               display=display)
  class(.ret) <- "nonmemTnpri"
  .ret
}

#' @export
print.nonmemTnpri <- function(x, ...) {
  .what <- switch(x$type,
                  msf=paste0("the NONMEM model specification file '", x$prior, "'"),
                  fit="a NONMEM fit of an nlmixr2 fit's data",
                  data="a NONMEM fit of a dataset")
  cat("NONMEM $PRIOR TNPRI from ", .what, "\n", sep="")
  cat("  ", .nonmemTnpriOptions(x), "\n", sep="")
  invisible(x)
}

#' The options written after `$PRIOR TNPRI`
#'
#' @param tnpri `nonmemTnpri` object
#' @return string
#' @noRd
#' @author Matthew L. Fidler
.nonmemTnpriOptions <- function(tnpri) {
  .ret <- "(PROBLEM 2)"
  if (!is.null(tnpri$plev)) .ret <- paste0(.ret, " PLEV=", format(tnpri$plev, digits=15))
  for (.n in c("ivar", "ityp", "ifnd", "mode")) {
    if (!is.null(tnpri[[.n]])) .ret <- paste0(.ret, " ", toupper(.n), "=", tnpri[[.n]])
  }
  if (isTRUE(tnpri$display)) .ret <- paste0(.ret, " DISPLAY")
  .ret
}

#' The TNPRI prior requested for this model, if any
#'
#' @param ui rxode2 ui
#' @return `nonmemTnpri` object or `NULL`
#' @noRd
#' @author Matthew L. Fidler
.nonmemTnpri <- function(ui) {
  .t <- rxode2::rxGetControl(ui, "tnpri", NULL)
  if (is.null(.t)) return(NULL)
  if (!inherits(.t, "nonmemTnpri")) .t <- nonmemTnpri(.t)
  if (.nonmemHasPriors(ui)) {
    stop("a NONMEM problem has only one $PRIOR: the `ini({})` priors ($PRIOR NWPRI) ",
         "cannot be combined with nonmemControl(tnpri=) ($PRIOR TNPRI)",
         call.=FALSE)
  }
  .t
}

#' The MSF file name `$MSFI` reads
#'
#' This is the name in the directory the control stream is written to;
#' the MSF is copied there before NONMEM is run.
#'
#' @param ui rxode2 ui
#' @return file name
#' @noRd
#' @author Matthew L. Fidler
.nonmemTnpriMsfName <- function(ui) {
  .msf <- rxode2::rxGetControl(ui, ".tnpriMsf", NULL)
  if (!is.null(.msf)) return(basename(.msf))
  .t <- .nonmemTnpri(ui)
  if (.t$type == "msf") return(basename(.t$prior))
  paste0(.nonmemTnpriPriorName(ui), ".msf")
}

#' Model name of the run that creates the prior
#'
#' @param ui rxode2 ui
#' @return model name
#' @noRd
#' @author Matthew L. Fidler
.nonmemTnpriPriorName <- function(ui) {
  paste0(rxUiGet.nonmemModelName(list(ui)), "_prior")
}

#' What has to be the same for the prior run to be a run of this model
#'
#' @param ui rxode2 ui
#' @return list
#' @noRd
#' @author Matthew L. Fidler
.nonmemTnpriModelSig <- function(ui) {
  .iniDf <- ui$iniDf
  list(model=deparse(ui$lstExpr),
       ini=.iniDf[, c("ntheta", "neta1", "neta2", "name", "fix")])
}

#' A copy of the ui that the prior run can change freely
#'
#' Creating the NONMEM files changes the ui's control values (the model
#' name, the data columns, ...), which would otherwise leak into the
#' fit the prior is for.
#'
#' @param ui rxode2 ui
#' @return copy of the ui
#' @noRd
#' @author Matthew L. Fidler
.nonmemCopyUi <- function(ui) {
  .ret <- new.env(parent=emptyenv())
  for (.n in ls(envir=ui, all.names=TRUE)) {
    assign(.n, get(.n, envir=ui), envir=.ret)
  }
  class(.ret) <- class(ui)
  .ret
}

#' Create (or find) the MSF of the prior run
#'
#' For a dataset or an nlmixr2 fit, the model is fit to that data with
#' NONMEM (with `$MSFO` and `$COVARIANCE`) in its own directory.  The
#' result is the path of the MSF, which may not exist yet when NONMEM is
#' not run by babelmixr2 (`runCommand=NA` or `run=FALSE`).
#'
#' @param env nlmixr2 estimation environment
#' @return path to the MSF, invisibly; also recorded in the ui's control
#' @noRd
#' @author Matthew L. Fidler
.nonmemTnpriPrepare <- function(env) {
  .ui <- env$ui
  .t <- .nonmemTnpri(.ui)
  if (is.null(.t)) return(invisible(NULL))
  if (.t$type == "msf") {
    if (!file.exists(.t$prior)) {
      warning("the TNPRI model specification file '", .t$prior, "' does not exist (yet)",
              call.=FALSE)
    }
    .msf <- .t$prior
    rxode2::rxAssignControlValue(.ui, ".tnpriMsf", .msf)
    # a different model specification file is a different fit
    rxode2::rxAssignControlValue(.ui, ".tnpriHash",
                                 if (file.exists(.msf)) digest::digest(.msf, file=TRUE)
                                 else normalizePath(.msf, mustWork=FALSE))
    return(invisible(.msf))
  }
  if (.t$type == "fit") {
    .fit <- .t$prior
    .fitUi <- .fit$ui
    if (!identical(.nonmemTnpriModelSig(.fitUi), .nonmemTnpriModelSig(.ui))) {
      stop("the TNPRI prior has to be a fit of the same model: NONMEM uses one model ",
           "code for both problems and lines the parameters up by position",
           call.=FALSE)
    }
    .data <- .fit$origData
  } else {
    .data <- .t$prior
  }
  .control <- .ui$control
  .priorControl <- .control
  .priorControl$tnpri <- NULL
  .priorControl$msfo <- TRUE
  .priorControl$modelName <- .nonmemTnpriPriorName(.ui)
  # TNPRI reads the covariance of the estimates from the MSF
  if (identical(.priorControl$cov, "")) .priorControl$cov <- "r,s"
  class(.priorControl) <- class(.control)
  .priorUi <- .nonmemCopyUi(.ui)
  assign("control", .priorControl, envir=.priorUi)
  # the prior run's own model number, not the one of this fit
  if (exists(".num", envir=.priorUi, inherits=FALSE)) rm(".num", envir=.priorUi)
  rxode2::rxAssignControlValue(.priorUi, ".modelNumber", 0)
  .priorEnv <- new.env(parent=emptyenv())
  .priorEnv$ui <- .priorUi
  .priorEnv$data <- .data
  .priorEnv$table <- env$table
  .minfo("TNPRI prior: NONMEM fit of the prior data")
  .exp <- .nonmemFamilyExport(.priorEnv)
  .msf <- file.path(.exp$exportPath, paste0(rxUiGet.nonmemModelName(list(.exp$ui)), ".msf"))
  if (.exp$ran) {
    if (!isTRUE(.exp$ui$nonmemSuccessful)) {
      stop("the NONMEM fit that creates the TNPRI prior was not successful; see '",
           .exp$nmctlFile, "'", call.=FALSE)
    }
    .cov <- file.path(.exp$exportPath, rxUiGet.nonmemCovFile(list(.exp$ui)))
    if (!file.exists(.cov)) {
      stop("the covariance step of the NONMEM fit that creates the TNPRI prior failed, ",
           "and TNPRI needs it; see '", .exp$nmctlFile, "'", call.=FALSE)
    }
    if (!file.exists(.msf)) {
      stop("the NONMEM fit that creates the TNPRI prior did not write '", .msf, "'",
           call.=FALSE)
    }
  } else {
    .minfo(paste0("run the TNPRI prior '", .exp$nmctlFile,
                  "' with NONMEM first; it creates '", basename(.msf), "'"))
  }
  rxode2::rxAssignControlValue(.ui, ".tnpriMsf", .msf)
  # the prior run's own control stream and data identify the prior, and
  # are known before NONMEM runs, so a fit exported before the prior was
  # run is still found when NONMEM's results are read back
  rxode2::rxAssignControlValue(.ui, ".tnpriHash", .exp$hash)
  invisible(.msf)
}
