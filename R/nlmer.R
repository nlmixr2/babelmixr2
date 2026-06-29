# -----------------------------------------------------------------------
# Global surrogate function for lme4::nlmer
# -----------------------------------------------------------------------

#' Surrogate gradient function for lme4::nlmer
#'
#' Computes predictions and the analytical gradient for every observation
#' using the nlm C machinery (a single multithreaded [nlmixr2est::nlmerSolveGrad()]
#' call against the model kept resident by [nlmixr2est::.nlmSetupEnv()]).
#' Parameters arrive in positional order matching
#' `.nlmerGlobal$nlmerEnv$paramNames` (THETA[i] order).
#'
#' lme4 passes one value per observation for each nonlinear parameter; within a
#' subject these are constant (phi = beta + b), so the per-subject parameter
#' vector is read from each subject's first observation row.
#'
#' @param TIME Observation times vector (unused; kept for formula
#'   compatibility -- the solver uses the resident per-subject event data)
#' @param ... Individual-level parameter vectors (one per param in
#'   `saemParamsToEstimate`)
#' @return Numeric vector of predictions with `"gradient"` attribute
#' @export
.nlmerNlmerFun <- function(TIME, ...) {
  .env <- .nlmerGlobal$nlmerEnv
  .paramsList <- list(...)
  .nPar <- .env$nPar
  .nSub <- .env$nSub
  .nObs <- .env$nObs
  .pNames <- .env$paramNames
  .obsIdxBySid <- .env$obsIdxBySid

  # Per-subject theta matrix in solver (first-appearance subject) order
  .thetaMat <- matrix(0.0, .nSub, .nPar)
  for (.i in seq_len(.nSub)) {
    .fi <- .obsIdxBySid[[.i]][1L]
    for (.j in seq_len(.nPar)) {
      .thetaMat[.i, .j] <- .paramsList[[.j]][.fi]
    }
  }

  # Single C solve: nObsTot x (nPar + 1); col 1 = rx_pred_, cols 2.. =
  # d(pred)/d(THETA[j]).  Rows are stacked per subject in solver block order.
  .mat <- nlmixr2est::nlmerSolveGrad(.thetaMat)

  .pred <- numeric(.nObs)
  .grad <- matrix(0.0, .nObs, .nPar)
  colnames(.grad) <- .pNames
  .row <- 0L
  for (.i in seq_len(.nSub)) {
    .idx <- .obsIdxBySid[[.i]]
    .rng <- .row + seq_along(.idx)
    .pred[.idx] <- .mat[.rng, 1L]
    .grad[.idx, ] <- .mat[.rng, -1L, drop = FALSE]
    .row <- .row + length(.idx)
  }
  attr(.pred, "gradient") <- .grad
  .pred
}

# -----------------------------------------------------------------------
# Solve control: derive the nlm-engine solving control from nlmerControl
# -----------------------------------------------------------------------

#' Build the nlm-machinery solving control for nlmer
#'
#' The nlmer engine drives optimization with `lme4::nlmer`, but the inner
#' per-subject prediction + analytic gradient is computed by the nlm C
#' machinery.  `lme4` passes raw per-subject phi (= beta + b), so the nlm
#' scaling must be neutralized (`scaleType = 5L`, "none") to feed those values
#' through unchanged.
#'
#' @param control nlmerControl object
#' @param np number of estimated structural parameters (THETA count)
#' @return A plain list usable as the `control` arg of
#'   [nlmixr2est::.nlmSetupEnv()]
#' @noRd
.nlmerSolveControl <- function(control, np) {
  .sc <- nlmixr2est::nlmControl(
    rxControl       = control$rxControl,
    eventSens       = control$eventSens,
    eventType       = control$eventType,
    shiErr          = control$shiErr,
    shi21maxFD      = control$shi21maxFD,
    stickyRecalcN   = control$stickyRecalcN,
    maxOdeRecalc    = control$maxOdeRecalc,
    odeRecalcFactor = control$odeRecalcFactor,
    optExpression   = control$optExpression,
    sumProd         = control$sumProd,
    solveType       = "grad",
    print           = 0L
  )
  .sc <- unclass(.sc)
  .sc$scaleType <- 5L                # 'none' -> identity (raw phi from lme4)
  .sc$scaleC <- rep(1.0, np)
  .sc
}

# -----------------------------------------------------------------------
# Fit model
# -----------------------------------------------------------------------

.nlmerFitModel <- function(ui, obsData) {
  .ctl <- ui$control
  .lme4Ctl <- lme4::nlmerControl(
    optimizer = .ctl$optimizer,
    tolPwrss  = .ctl$tolPwrss,
    optCtrl   = .ctl$optCtrl
  )
  .obsData <- obsData
  .obsData$ID <- factor(.obsData$ID)
  .start <- ui$nlmerStart
  # lme4::nlmer captures match.call()$formula and re-evaluates it deep in its
  # call stack (where `ui` is not in scope). Storing the formula in the
  # package-level .nlmerGlobal environment and referencing it via :: makes it
  # resolvable from any evaluation context.
  .nlmerGlobal$currentFormula <- ui$nlmerFormula
  on.exit(.nlmerGlobal$currentFormula <- NULL, add = TRUE)
  eval(
    quote(lme4::nlmer(
      formula = babelmixr2:::.nlmerGlobal$currentFormula,
      data    = .obsData,
      start   = .start,
      control = .lme4Ctl
    )),
    envir = environment()
  )
}

# -----------------------------------------------------------------------
# Parameter extraction from lme4::nlmer fit
# -----------------------------------------------------------------------

#' Get the full theta (structural + residual) from a lme4 nlmer fit
#'
#' @param nlmerMod lme4::nlmer fit object
#' @param ui rxode2 ui object
#' @return Named numeric theta vector
#' @noRd
.nlmerGetTheta <- function(nlmerMod, ui) {
  .f <- lme4::fixef(nlmerMod)
  .predDf <- ui$predDf
  .errType <- .predDf$errType
  .sigma <- sigma(nlmerMod)
  if (.errType == "add") {
    .w <- which(ui$iniDf$err == "add")
    if (length(.w) == 1L) return(c(.f, setNames(.sigma, ui$iniDf$name[.w])))
  } else if (.errType == "prop") {
    .w <- which(ui$iniDf$err == "prop")
    if (length(.w) == 1L) return(c(.f, setNames(.sigma, ui$iniDf$name[.w])))
  }
  .f
}

#' Get the random-effects matrix as an eta matrix for nlmixr2
#'
#' @param nlmerMod lme4::nlmer fit object
#' @param ui rxode2 ui object
#' @return Numeric matrix (nSubj x nEta), row-ordered by subject ID
#' @noRd
.nlmerGetEtaMat <- function(nlmerMod, ui) {
  .re <- lme4::ranef(nlmerMod)$ID
  if (is.null(.re)) return(matrix(nrow = 0, ncol = 0))
  .re <- .re[order(as.numeric(row.names(.re))), , drop = FALSE]
  names(.re) <- nlmixr2est:::.nlmeGetNonMuRefNames(names(.re), ui)
  row.names(.re) <- NULL
  as.matrix(.re)
}

#' Get the random-effects covariance as an omega matrix for nlmixr2
#'
#' @param nlmerMod lme4::nlmer fit object
#' @param ui rxode2 ui object
#' @return Named square covariance matrix
#' @noRd
.nlmerGetOmega <- function(nlmerMod, ui) {
  .vc <- lme4::VarCorr(nlmerMod)$ID
  if (is.null(.vc)) {
    .ome <- ui$omega * 0.0
    return(.ome)
  }
  .sdv <- attr(.vc, "stddev")
  .corr <- as.matrix(.vc)
  .ome <- diag(.sdv, nrow = length(.sdv)) %*% .corr %*% diag(.sdv, nrow = length(.sdv))
  # matrix %*% propagates rownames only from left; use corr's names directly
  .pnames <- colnames(.corr)
  if (is.null(.pnames)) .pnames <- names(.sdv)
  .name <- nlmixr2est:::.nlmeGetNonMuRefNames(.pnames, ui)
  dimnames(.ome) <- list(.name, .name)
  # Return only the estimated random effects block (may be a subset of omega)
  .ome
}

#' Covariance matrix of fixed effects from nlmer
#'
#' @param nlmerMod lme4::nlmer fit object
#' @return Square covariance matrix
#' @noRd
.nlmerGetCov <- function(nlmerMod) {
  as.matrix(vcov(nlmerMod))
}

# -----------------------------------------------------------------------
# Family-level fit (assembles the nlmixr2 result environment)
# -----------------------------------------------------------------------

.nlmerFamilyFit <- function(env, ...) {
  .ui <- env$ui
  .control <- .ui$control
  .data <- env$data
  .ret <- new.env(parent = emptyenv())
  .ret$table <- env$table

  nlmixr2est::.foceiPreProcessData(.data, .ret, .ui, .control$rxControl)

  # Compile the nlmer sensitivity model and load it (with its event data) into
  # the nlm C machinery, keeping it resident for every lme4 PWRSS evaluation.
  # nlmerSolveGrad() then services each evaluation with a single multithreaded
  # C solve (per-subject prediction + analytic gradient + jump sensitivities).
  .sm <- .ui$nlmerSensModel
  .modelInfo <- list(
    predOnly   = .sm$predOnly,
    thetaGrad  = .sm$thetaGrad,
    eventTheta = .sm$eventTheta
  )
  .start <- .ui$nlmerStart
  .np <- length(.start)
  .solveCtl <- .nlmerSolveControl(.control, .np)
  nlmixr2est::.nlmSetupEnv(.start, .ui, .ret$dataSav, .modelInfo, .solveCtl)
  on.exit(nlmixr2est::.nlmFreeEnv(), add = TRUE)

  # Observation-row map for scattering the solver output back to lme4's data
  # order.  The solver stacks subjects in first-appearance order (== rxode2's
  # internal id order) and, within a subject, in solve order -- the same
  # within-subject order as the obs frame passed to lme4.
  .dsAll <- .ret$dataSav[.ret$dataSav$EVID != 2, ]
  .obsData <- .dsAll[.dsAll$EVID == 0, ]
  .sidOrder <- unique(.dsAll$ID)
  .obsIdxBySid <- lapply(.sidOrder, function(.s) which(.obsData$ID == .s))

  .nlmerGlobal$nlmerEnv <- new.env(parent = emptyenv())
  .nlmerGlobal$nlmerEnv$paramNames <- .sm$paramNames
  .nlmerGlobal$nlmerEnv$nPar <- .np
  .nlmerGlobal$nlmerEnv$nSub <- length(.sidOrder)
  .nlmerGlobal$nlmerEnv$nObs <- nrow(.obsData)
  .nlmerGlobal$nlmerEnv$obsIdxBySid <- .obsIdxBySid
  on.exit(
    {
      .nlmerGlobal$nlmerEnv <- NULL
    },
    add = TRUE
  )

  .nlmer <- nlmixr2est::.collectWarn(.nlmerFitModel(.ui, .obsData), lst = TRUE)
  .ret$nlmer <- .nlmer[[1]]

  .ret$message <- NULL
  lapply(.nlmer[[2]], function(.w) {
    warning(.w, call. = FALSE)
    if (grepl("fail|did not converge", .w, ignore.case = TRUE)) {
      .ret$message <- c(.ret$message, paste0(.w, " (carefully review results)"))
    }
  })
  if (is.null(.ret$message)) {
    .ret$message <- ""
  } else {
    .ret$message <- paste(.ret$message, collapse = "\n")
  }

  if (rxode2::rxGetControl(.ui, "returnNlmer", FALSE)) {
    return(.ret$nlmer)
  }

  .ret$ui <- .ui
  .ret$adjObf <- rxode2::rxGetControl(.ui, "adjObf", TRUE)
  .ret$fullTheta <- .nlmerGetTheta(.ret$nlmer, .ui)
  .ret$cov <- .nlmerGetCov(.ret$nlmer)
  .ret$covMethod <- "nlmer"
  .ret$etaMat <- .nlmerGetEtaMat(.ret$nlmer, .ui)
  if (nrow(.ret$etaMat) > 0) {
    .ret$etaObf <- data.frame(
      ID   = seq_len(nrow(.ret$etaMat)),
      as.data.frame(.ret$etaMat),
      OBJI = NA_real_
    )
  } else {
    .ret$etaObf <- data.frame(ID = integer(0), OBJI = numeric(0))
  }
  .ret$omega <- .nlmerGetOmega(.ret$nlmer, .ui)
  .ret$control <- .control
  .ret$extra <- paste0(" by ", crayon::bold$yellow("maximum likelihood (lme4::nlmer)"))

  nlmixr2est::.nlmixr2FitUpdateParams(.ret)
  nlmixr2est::nmObjHandleControlObject(.ret$control, .ret)
  if (exists("control", .ui)) rm(list = "control", envir = .ui)

  .ret$est <- "nlmer"
  .ret$objective <- -2.0 * as.numeric(logLik(.ret$nlmer))
  .ret$model <- .ui$ebe
  .ret$ofvType <- "nlmer"
  .nlmerControlToFoceiControl(.ret)
  .ret$theta <- .ret$ui$saemThetaDataFrame

  .ret <- nlmixr2est::nlmixr2CreateOutputFromUi(
    .ret$ui,
    data    = .ret$origData,
    control = .ret$control,
    table   = .ret$table,
    env     = .ret,
    est     = "nlmer"
  )
  .env <- .ret$env
  .env$method <- "nlmer"
  .ret
}

# -----------------------------------------------------------------------
# nlmixr2Est dispatch
# -----------------------------------------------------------------------

#' nlmixr2 estimation using lme4::nlmer
#'
#' @inheritParams nlmixr2est::nlmixr2Est
#' @export
nlmixr2Est.nlmer <- function(env, ...) {
  .ui <- env$ui
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("estimation 'nlmer' requires the 'lme4' package", call. = FALSE)
  }
  rxode2::assertRxUiMixedOnly(
    .ui, " for the estimation routine 'nlmer', try 'focei'",
    .var.name = .ui$modelName
  )
  rxode2::assertRxUiNormal(
    .ui, " for the estimation routine 'nlmer'",
    .var.name = .ui$modelName
  )
  rxode2::assertRxUiSingleEndpoint(
    .ui, " for the estimation routine 'nlmer'",
    .var.name = .ui$modelName
  )
  rxode2::assertRxUiRandomOnIdOnly(
    .ui, " for the estimation routine 'nlmer'",
    .var.name = .ui$modelName
  )
  rxode2::assertRxUiEstimatedResiduals(
    .ui, " for the estimation routine 'nlmer'",
    .var.name = .ui$modelName
  )

  .nlmerFamilyControl(env, ...)
  on.exit(
    {
      if (exists("control", envir = .ui)) rm("control", envir = .ui)
    },
    add = TRUE
  )

  .model <- nlmixr2est::.uiApplyMu2(env)
  nlmixr2est::.uiFinalizeMu2(.nlmerFamilyFit(env, ...), .model)
}

attr(nlmixr2Est.nlmer, "covPresent") <- TRUE
attr(nlmixr2Est.nlmer, "unbounded")  <- TRUE
attr(nlmixr2Est.nlmer, "mu") <- function(control) {
  isTRUE(control$muRefCovAlg)
}
