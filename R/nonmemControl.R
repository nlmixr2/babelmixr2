#' NONMEM estimation control
#'
#' @details
#'
#' If \code{runCommand} is given as a string, it will be called with the
#' \code{system()} command like:
#'
#' \code{runCommand controlFile outputFile}.
#'
#' For example, if \code{runCommand="'/path/to/nmfe75'"} then the command line
#' used would look like the following:
#'
#' \code{'/path/to/nmfe75' one.cmt.nmctl one.cmt.lst}
#'
#' If \code{runCommand} is given as a function, it will be called as
#' \code{FUN(ctl, directory, ui)} to run NONMEM.  This allows you to run NONMEM
#' in any way that you may need, as long as you can write it in R.  babelmixr2
#' will wait for the function to return before proceeding.
#'
#' If \code{runCommand} is \code{NA}, \code{nlmixr()} will stop after writing
#' the model files and without starting NONMEM.
#'
#' @param est NONMEM estimation method
#' @param advanOde The ODE solving method for NONMEM
#' @param cov The NONMEM covariance method
#' @param maxeval NONMEM's maxeval (for non posthoc methods)
#' @param tol NONMEM tolerance for ODE solving advan
#' @param atol NONMEM absolute tolerance for ODE solving
#' @param sstol NONMEM tolerance for steady state ODE solving
#' @param ssatol NONMEM absolute tolerance for steady state ODE solving
#' @param sigl NONMEM sigl estimation option
#' @param sigdig the significant digits for NONMEM
#' @param print The print number for NONMEM
#' @param extension NONMEM file extensions
#' @param outputExtension Extension to use for the NONMEM output
#'   listing
#' @param runCommand Command to run NONMEM (typically the path to "nmfe75") or a
#'   function.  See the details for more information.
#' @param iniSigDig How many significant digits are printed in $THETA
#'   and $OMEGA when the estimate is zero.  Also controls the zero
#'   protection numbers
#' @param protectZeros Add methods to protect divide by zero
#' @param muRef Automatically mu-reference the control stream
#' @param rxControl Options to pass to \code{rxode2::rxControl} for
#'   simulations
#' @param mapiter the number of map iterations for IMP method
#' @param niter number of iterations in NONMEM estimation methods
#' @param isample Isample argument for NONMEM  ITS estimation method
#' @param iaccept Iaccept for NONMEM ITS estimation methods
#' @param iscaleMin parameter for IMP NONMEM method (ISCALE_MIN)
#' @param iscaleMax parameter for IMP NONMEM method (ISCALE_MAX)
#' @param df degrees of freedom for IMP method
#' @param seed is the seed for NONMEM methods
#' @param mapinter is the MAPINTER parameter for the IMP method
#' @param  addProp,sumProd,optExpression,calcTables,compress,ci,sigdigTable
#'   Passed to \code{nlmixr2est::foceiControl}
#' @param readRounding Try to read NONMEM output when NONMEM
#'   terminated due to rounding errors
#' @param readBadOpt Try to read NONMEM output when NONMEM terminated
#'   due to an apparent failed optimization
#' @param noabort Add the `NOABORT` option for `$EST`
#' @param modelName Model name used to generate the NONMEM output.  If
#'   `NULL` try to infer from the model name (could be `x` if not
#'   clear).  Otherwise use this character for outputs.
#' @param ... optional \code{genRxControl} argument controlling
#'   automatic \code{rxControl} generation.
#'
#' @return babelmixr2 control option for generating NONMEM control stream and
#'   reading it back into `babelmixr2`/`nlmixr2`
#'
#' @inheritParams nlmixr2est::saemControl
#'
#' @author Matthew L. Fidler
#'
#' @export
#'
#' @examples
#' nonmemControl()
nonmemControl <- function(est=c("focei", "imp", "its", "posthoc"),
                          advanOde=c("advan13", "advan8", "advan6"),
                          cov=c("r,s", "r", "s", ""),
                          maxeval=100000,
                          tol=6,
                          atol=12,
                          sstol=6,
                          ssatol=12,
                          sigl=12,
                          sigdig=3,
                          print=1,
                          extension=getOption("babelmixr2.nmModelExtension", ".nmctl"),
                          outputExtension=getOption("babelmixr2.nmOutputExtension", ".lst"),
                          runCommand=getOption("babelmixr2.nonmem", ""),
                          iniSigDig=5,
                          protectZeros=TRUE,
                          muRef=TRUE,
                          addProp = c("combined2", "combined1"),
                          rxControl=NULL,
                          sumProd = FALSE,
                          optExpression = TRUE,
                          calcTables = TRUE,
                          compress = TRUE,
                          ci = 0.95,
                          sigdigTable=NULL,
                          readRounding=FALSE,
                          readBadOpt=FALSE,
                          niter=100L,
                          isample=1000L,
                          iaccept=0.4,
                          iscaleMin=0.1,
                          iscaleMax=10,
                          df=4,
                          seed=14456,
                          mapiter=1,
                          mapinter=0,
                          noabort=TRUE,
                          modelName=NULL,
                          muRefCovAlg=TRUE,
                          ...) {
  # nonmem manual slides suggest tol=6, sigl=6 sigdig=2
  checkmate::assertIntegerish(maxeval, lower=100, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(sigdig, lower=1, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(print, lower=1, len=1, any.missing=FALSE)
  checkmate::assertLogical(noabort, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(tol, lower=1, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(atol, lower=1, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(sstol, lower=1, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(ssatol, lower=1, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(sigl, lower=1, upper=14, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(iniSigDig, lower=1, len=1, any.missing=FALSE)
  checkmate::assertLogical(protectZeros, len=1, any.missing=FALSE)
  checkmate::assertLogical(muRef, len=1, any.missing=FALSE)
  checkmate::assertLogical(readRounding, len=1, any.missing=FALSE)
  checkmate::assertLogical(readBadOpt, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(niter, lower=1, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(isample, lower=1, len=1, any.missing=FALSE)
  checkmate::assertNumeric(iaccept, lower=0, upper=1, len=1, any.missing=FALSE)
  checkmate::assertNumeric(iscaleMin, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(iscaleMax, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(df, lower=0, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(seed, lower=1, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(mapiter, len=1, any.missing=FALSE)
  checkmate::assertLogical(muRefCovAlg, any.missing=FALSE, len=1)
  if (!is.null(modelName)) {
    checkmate::assertCharacter(modelName, len=1, any.missing=FALSE)
  }
  if (!identical(runCommand, "")) {
    if (!(checkmate::testCharacter(runCommand, len=1) ||
          checkmate::testFunction(runCommand, args=c("ctl", "directory", "ui")))) {
      stop("runCommand must be a character string or a function with arguments 'ctl', 'directory', and 'ui'")
    }
  }
  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% c("genRxControl"))]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste
    (paste0("'", .bad, "'", sep=""), collapse=", "),
    call.=FALSE)
  }

  if (checkmate::testIntegerish(addProp, lower=1, upper=1, len=1)) {
    addProp <- c("combined1", "combined2")[addProp]
  } else {
    addProp <- match.arg(addProp)
  }
  checkmate::assertLogical(compress, any.missing=FALSE, len=1)

  if (!is.null(.xtra$genRxControl)) {
    genRxControl <- .xtra$genRxControl
  } else {
    genRxControl <- FALSE
    if (is.null(rxControl)) {
      rxControl <- rxode2::rxControl(
        rtol=10^(-tol),
        atol=10^(-atol),
        ssRtol=10^(-sstol),
        ssAtol=10^(-ssatol),
        covsInterpolation="nocb",
        method="liblsoda"
      )
      genRxControl <- TRUE
    } else if (is.list(rxControl)) {
      rxControl$rtol <- 10^(-tol)
      rxControl$atol <- 10^(-atol)
      rxControl$ssRtol <- 10^(-sstol)
      rxControl$ssAtol <- 10^(-ssatol)
      rxControl$covsInterpolation <- "nocb"
      rxControl$method <- "liblsoda"
      rxControl <- do.call(rxode2::rxControl, rxControl)
    }
    if (!inherits(rxControl, "rxControl")) {
      stop("rxControl needs to be ode solving options from rxode2::rxControl()",
           call.=FALSE)
    }
  }

  checkmate::assertLogical(sumProd, any.missing=FALSE, len=1)
  checkmate::assertLogical(optExpression, any.missing=FALSE, len=1)
  checkmate::assertNumeric(ci, any.missing=FALSE, len=1, lower=0, upper=1)
  checkmate::assertLogical(calcTables, len=1, any.missing=FALSE)

  .ret <- list(est=match.arg(est),
               cov=match.arg(cov),
               advanOde=match.arg(advanOde),
               maxeval=maxeval,
               print=print,
               noabort=noabort,
               iniSigDig=iniSigDig,
               tol=tol,
               atol=atol,
               sstol=sstol,
               ssatol=ssatol,
               sigl=sigl,
               muRef=muRef,
               sigdig=sigdig,
               protectZeros=protectZeros,
               runCommand=runCommand,
               outputExtension=outputExtension,
               addProp=addProp,
               rxControl=rxControl,
               sumProd = sumProd,
               optExpression=optExpression,
               calcTables = calcTables,
               compress = compress,
               ci = ci,
               sigdigTable=sigdigTable,
               readRounding=readRounding,
               readBadOpt=readBadOpt,
               genRxControl=genRxControl,
               niter=niter,
               isample=isample,
               iaccept=iaccept,
               iscaleMin=iscaleMin,
               iscaleMax=iscaleMax,
               df=df,
               seed=seed,
               mapiter=mapiter,
               modelName=modelName,
               muRefCovAlg=muRefCovAlg
               )
  class(.ret) <- "nonmemControl"
  .ret
}

#' @export
getValidNlmixrCtl.nonmem <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- nonmemControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("nonmemControl", .ctl)
  if (!inherits(.ctl, "nonmemControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- nonmemControl()
  } else {
    .ctl <- do.call(nonmemControl, .ctl)
  }
  .ctl
}

#' @export
nmObjHandleControlObject.nonmemControl <- function(control, env) {
  assign("nonmemControl", control, envir=env)
}

#' @export
nmObjGetControl.nonmem <- function(x, ...) {
  .env <- x[[1]]
  if (exists("nonmemControl", .env)) {
    .control <- get("nonmemControl", .env)
    if (inherits(.control, "nonmemControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "nonmemControl")) return(.control)
  }
  stop("cannot find nonmem related control object", call.=FALSE)
}

#' @export
nmObjHandleControlObject.nonmemControl <- function(control, env) {
  assign("nonmemControl", control, envir=env)
}

#' @export
nmObjGetControl.nonmem <- function(x, ...) {
  .env <- x[[1]]
  if (exists("nonmemControl", .env)) {
    .control <- get("nonmemControl", .env)
    if (inherits(.control, "nonmemControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "nonmemControl")) return(.control)
  }
  stop("cannot find nonmem related control object", call.=FALSE)
}

.nonmemControlToFoceiControl <- function(env, assign=FALSE) {
  .nonmemControl <- env$nonmemControl
  .ui <- env$ui
  .foceiControl <- nlmixr2est::foceiControl(rxControl = env$nonmemControl$rxControl,
                                            maxOuterIterations = 0L, maxInnerIterations = 0L, covMethod = 0L,
                                            etaMat = env$etaMat, sumProd = .nonmemControl$sumProd,
                                            optExpression = .nonmemControl$optExpression, scaleTo = 0,
                                            calcTables = .nonmemControl$calcTables,
                                            addProp = .nonmemControl$addProp,
                                            skipCov = .ui$foceiSkipCov, interaction = 1L,
                                            compress = .nonmemControl$compress,
                                            ci = .nonmemControl$ci,
                                            sigdigTable = .nonmemControl$sigdigTable)
  if (assign)
    env$control <- .foceiControl
  .foceiControl
}
