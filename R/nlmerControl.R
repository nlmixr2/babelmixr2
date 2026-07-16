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
#' @param eventSens method used for the dosing-parameter (alag/F/rate/dur)
#'   sensitivities in the analytic gradient model: `"jump"` routes them
#'   through rxode2's analytic event jumps; `"fd"` falls back to Shi2021
#'   finite differences.  See [nlmixr2est::nlmControl()].
#' @param eventType finite-difference type (`"central"` or `"forward"`)
#'   used for event-related parameters when `eventSens = "fd"`
#' @param shiErr epsilon used when optimizing the ideal finite-difference
#'   step size with the Shi2021 method
#' @param shi21maxFD maximum number of Shi2021 step-size optimization
#'   iterations for the gradient
#' @param stickyRecalcN number of bad ODE solves tolerated before the
#'   per-subject tolerance is stickily relaxed
#' @param maxOdeRecalc maximum number of times to retry a bad ODE solve
#'   with a relaxed tolerance
#' @param odeRecalcFactor factor by which atol/rtol are relaxed on an ODE
#'   solve retry
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
                                # nlm-machinery solving options (per-subject
                                # prediction + analytic gradient are computed in
                                # C via the nlm engine; see .nlmerSolveControl)
                                eventSens = c("jump", "fd"),
                                eventType = c("central", "forward"),
                                shiErr = (.Machine$double.eps)^(1 / 3),
                                shi21maxFD = 20L,
                                stickyRecalcN = 4,
                                maxOdeRecalc = 5,
                                odeRecalcFactor = 10^(0.5),
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
  checkmate::assertNumeric(shiErr, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(shi21maxFD, lower = 1, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(stickyRecalcN, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(maxOdeRecalc, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(odeRecalcFactor, lower = 1, any.missing = FALSE, len = 1)

  eventSens <- match.arg(eventSens)
  eventType <- match.arg(eventType)
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
    eventSens = eventSens,
    eventType = eventType,
    shiErr = shiErr,
    shi21maxFD = as.integer(shi21maxFD),
    stickyRecalcN = as.integer(stickyRecalcN),
    maxOdeRecalc = as.integer(maxOdeRecalc),
    odeRecalcFactor = odeRecalcFactor,
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
    rxode2::.minfo("invalid control for `est=\"nlmer\"`, using default")
    .ctl <- nlmerControl()
  } else {
    .ctl <- do.call(nlmerControl, .ctl)
  }
  .ctl
}

.nlmerFamilyControl <- function(env, ...) {
  nlmixr2est::.nlmFamilyControlGeneric(env, nlmerControl, "nlmerControl")
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
