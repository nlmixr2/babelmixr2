# -----------------------------------------------------------------------
# Global surrogate function for lme4::nlmer
# -----------------------------------------------------------------------

#' Surrogate gradient function for lme4::nlmer
#'
#' Computes predictions and analytical gradient using rxode2 sensitivity
#' equations. Parameters come in positional order matching
#' `.nlmerGlobal$nlmerEnv$paramNames`.
#'
#' @param TIME Observation times vector
#' @param ... Individual-level parameter vectors (one per param in
#'   `saemParamsToEstimate`)
#' @return Numeric vector of predictions with `"gradient"` attribute
#' @export
.nlmerNlmerFun <- function(TIME, ...) {
  .env <- .nlmerGlobal$nlmerEnv
  .paramsList <- list(...)
  .nPar <- length(.paramsList)
  .nObs <- length(TIME)
  .pNames <- .env$paramNames
  names(.paramsList) <- .pNames
  .obsSubjId <- .env$obsSubjId
  .uniqueSubjs <- unique(.obsSubjId)
  .pred <- numeric(.nObs)
  .grad <- matrix(0.0, .nObs, .nPar)
  colnames(.grad) <- .pNames
  .rxCtrl <- .env$rxControl
  .sensModel <- .env$sensModel

  for (.sid in .uniqueSubjs) {
    .idx <- which(.obsSubjId == .sid)
    .firstIdx <- .idx[1]
    .theta <- setNames(
      vapply(seq_len(.nPar), function(.j) .paramsList[[.j]][.firstIdx], double(1)),
      paste0("THETA[", seq_len(.nPar), "]")
    )
    .evData <- .env$subjEvData[[as.character(.sid)]]
    .sol <- do.call(
      rxode2::rxSolve,
      c(list(object = .sensModel$thetaGrad,
             params = .theta,
             events = .evData),
        .rxCtrl)
    )
    .pred[.idx] <- .sol$rx_pred_
    for (.j in seq_len(.nPar)) {
      .sn <- paste0("rx__sens_rx_pred__BY_THETA_", .j, "___")
      .grad[.idx, .j] <- .sol[[.sn]]
    }
  }
  attr(.pred, "gradient") <- .grad
  .pred
}

# -----------------------------------------------------------------------
# Data setup
# -----------------------------------------------------------------------

#' Set up per-subject event data and obs-subject mapping for nlmer
#'
#' @param dataSav Pre-processed data from `.foceiPreProcessData`
#' @return Observation-only data frame (EVID == 0 rows)
#' @noRd
.nlmerFitDataSetup <- function(dataSav) {
  .dsAll <- dataSav[dataSav$EVID != 2, ]
  .ids <- unique(.dsAll$ID)
  .subjEvData <- setNames(
    lapply(as.character(.ids), function(.id) {
      .dsAll[.dsAll$ID == as.integer(.id), ]
    }),
    as.character(.ids)
  )
  .obsData <- .dsAll[.dsAll$EVID == 0, ]
  .nlmerGlobal$nlmerEnv$subjEvData <- .subjEvData
  .nlmerGlobal$nlmerEnv$obsSubjId <- .obsData$ID
  .obsData
}

# -----------------------------------------------------------------------
# Fit model
# -----------------------------------------------------------------------

.nlmerFitModel <- function(ui, dataSav) {
  .obsData <- .nlmerFitDataSetup(dataSav)
  .ctl <- ui$control
  .lme4Ctl <- lme4::nlmerControl(
    optimizer = .ctl$optimizer,
    tolPwrss  = .ctl$tolPwrss,
    optCtrl   = .ctl$optCtrl
  )
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

  # Compile sensitivity model and populate the per-solve global env
  .nlmerGlobal$nlmerEnv <- new.env(parent = emptyenv())
  .sm <- .ui$nlmerSensModel
  .nlmerGlobal$nlmerEnv$sensModel <- .sm
  .nlmerGlobal$nlmerEnv$paramNames <- .sm$paramNames
  .nlmerGlobal$nlmerEnv$rxControl <- as.list(.control$rxControl)
  on.exit(
    {
      .nlmerGlobal$nlmerEnv <- NULL
    },
    add = TRUE
  )

  .nlmer <- nlmixr2est::.collectWarn(.nlmerFitModel(.ui, .ret$dataSav), lst = TRUE)
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
