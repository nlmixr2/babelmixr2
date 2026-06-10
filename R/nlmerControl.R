#' Control for nlmer estimation method in nlmixr2
#'
#' @inheritParams nlmixr2est::foceiControl
#' @inheritParams nlmixr2est::saemControl
#' @inheritParams nlmixr2est::nlmeControl
#'
#' @param optimizer passed to [lme4::nlmerControl()]: optimizer to use
#' @param restart_edge passed to [lme4::nlmerControl()]
#' @param boundary.tol passed to [lme4::nlmerControl()]
#' @param calc.derivs passed to [lme4::nlmerControl()]
#' @param use.last.params passed to [lme4::nlmerControl()]
#' @param sparseX passed to [lme4::nlmerControl()]
#' @param returnNlmer logical; when `TRUE` return the raw lme4 nlmerMod
#'   object instead of the nlmixr2 fit
#' @param muRefCovAlg logical; when `TRUE` apply mu2 covariate
#'   referencing algebraic transformation (same as saem)
#'
#' @return nlmer control structure
#' @export
#' @author Matthew L. Fidler
#' @examples
#' nlmerControl()
#' nlmixr2NlmerControl()
nlmixr2NlmerControl <- function(optimizer = "bobyqa",
                                tolPwrss = 1e-7,
                                optCtrl = list(),
                                returnNlmer = FALSE,
                                muRefCovAlg = TRUE,
                                # nlmixr2 standard options
                                optExpression = TRUE,
                                literalFix = TRUE,
                                sumProd = FALSE,
                                rxControl = NULL,
                                calcTables = TRUE,
                                compress = TRUE,
                                adjObf = TRUE,
                                ci = 0.95,
                                sigdig = 4,
                                sigdigTable = NULL,
                                addProp = c("combined2", "combined1"),
                                ...) {
  checkmate::assertLogical(optExpression, len = 1, any.missing = FALSE)
  checkmate::assertLogical(literalFix, len = 1, any.missing = FALSE)
  checkmate::assertLogical(sumProd, len = 1, any.missing = FALSE)
  checkmate::assertLogical(returnNlmer, len = 1, any.missing = FALSE)
  checkmate::assertLogical(muRefCovAlg, len = 1, any.missing = FALSE)
  checkmate::assertLogical(calcTables, len = 1, any.missing = FALSE)
  checkmate::assertLogical(compress, len = 1, any.missing = TRUE)
  checkmate::assertLogical(adjObf, len = 1, any.missing = TRUE)
  checkmate::assertNumeric(tolPwrss, len = 1, any.missing = FALSE, lower = 0)
  checkmate::assertNumeric(ci, lower = 0, upper = 1, any.missing = FALSE, len = 1)

  addProp <- match.arg(addProp)

  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% c("genRxControl"))]
  if (length(.bad) > 0) {
    stop("unused argument: ",
         paste(paste0("'", .bad, "'"), collapse = ", "),
         call. = FALSE)
  }

  .genRxControl <- FALSE
  if (!is.null(.xtra$genRxControl)) {
    .genRxControl <- .xtra$genRxControl
  }
  if (is.null(rxControl)) {
    if (!is.null(sigdig)) {
      rxControl <- rxode2::rxControl(sigdig = sigdig)
    } else {
      rxControl <- rxode2::rxControl(atol = 1e-4, rtol = 1e-4)
    }
    .genRxControl <- TRUE
  } else if (inherits(rxControl, "rxControl")) {
  } else if (is.list(rxControl)) {
    rxControl <- do.call(rxode2::rxControl, rxControl)
  } else {
    stop("solving options 'rxControl' needs to be generated from 'rxode2::rxControl'",
         call. = FALSE)
  }

  if (!is.null(sigdig)) {
    checkmate::assertNumeric(sigdig, lower = 1, finite = TRUE, any.missing = TRUE, len = 1)
    if (is.null(sigdigTable)) {
      sigdigTable <- round(sigdig)
    }
  }
  if (is.null(sigdigTable)) {
    sigdigTable <- 3
  }
  checkmate::assertIntegerish(sigdigTable, lower = 1, len = 1, any.missing = FALSE)

  .ret <- list(
    optimizer = optimizer,
    tolPwrss = tolPwrss,
    optCtrl = optCtrl,
    returnNlmer = returnNlmer,
    muRefCovAlg = muRefCovAlg,
    optExpression = optExpression,
    literalFix = literalFix,
    sumProd = sumProd,
    rxControl = rxControl,
    calcTables = calcTables,
    compress = compress,
    adjObf = adjObf,
    ci = ci,
    sigdig = sigdig,
    sigdigTable = sigdigTable,
    addProp = addProp,
    genRxControl = .genRxControl
  )
  class(.ret) <- "nlmerControl"
  .ret
}

#' @rdname nlmixr2NlmerControl
#' @export
nlmerControl <- nlmixr2NlmerControl

#' @export
rxUiDeparse.nlmerControl <- function(object, var) {
  .default <- nlmerControl()
  .w <- nlmixr2est::.deparseDifferent(.default, object, "genRxControl")
  nlmixr2est::.deparseFinal(.default, object, .w, var)
}

#' @export
nmObjHandleControlObject.nlmerControl <- function(control, env) {
  assign("nlmerControl", control, envir = env)
}

#' @export
nmObjGetControl.nlmer <- function(x, ...) {
  .env <- x[[1]]
  if (exists("nlmerControl", .env)) {
    .control <- get("nlmerControl", .env)
    if (inherits(.control, "nlmerControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "nlmerControl")) return(.control)
  }
  stop("cannot find nlmer related control object", call. = FALSE)
}

#' @export
getValidNlmixrCtl.nlmer <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- nlmerControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("nlmerControl", .ctl)
  if (!inherits(.ctl, "nlmerControl")) {
    nlmixr2est::.minfo("invalid control for `est=\"nlmer\"`, using default")
    .ctl <- nlmerControl()
  } else {
    .ctl <- do.call(nlmerControl, .ctl)
  }
  .ctl
}

.nlmerFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- nlmerControl()
  }
  if (!inherits(.control, "nlmerControl")) {
    .control <- do.call(nlmerControl, .control)
  }
  assign("control", .control, envir = .ui)
}

.nlmerControlToFoceiControl <- function(env, assign = TRUE) {
  .nlmerControl <- env$nlmerControl
  .ui <- env$ui
  .foceiControl <- nlmixr2est::foceiControl(
    rxControl = env$nlmerControl$rxControl,
    maxOuterIterations = 0L,
    maxInnerIterations = 0L,
    covMethod = 0L,
    etaMat = env$etaMat,
    sumProd = .nlmerControl$sumProd,
    optExpression = .nlmerControl$optExpression,
    literalFix = .nlmerControl$literalFix,
    literalFixRes = FALSE,
    scaleTo = 0,
    calcTables = .nlmerControl$calcTables,
    addProp = .nlmerControl$addProp,
    skipCov = .ui$foceiSkipCov,
    interaction = 1L,
    compress = .nlmerControl$compress,
    ci = .nlmerControl$ci,
    sigdigTable = .nlmerControl$sigdigTable,
    indTolRelax = TRUE
  )
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @export
nmObjGetFoceiControl.nlmer <- function(x, ...) {
  .nlmerControlToFoceiControl(x[[1]])
}
