# -----------------------------------------------------------------------
# Parameter / theta helpers
# -----------------------------------------------------------------------

#' Get the structural params to estimate for nlmer (= saemParamsToEstimate)
#'
#' @param x rxUiGet list(ui)
#' @param ... additional arguments (currently ignored)
#' @return Character vector of parameter names (mu-ref thetas + nonMuEtas)
#' @export
rxUiGet.nlmerThetas <- function(x, ...) {
  x[[1]]$saemParamsToEstimate
}
attr(rxUiGet.nlmerThetas, "rstudio") <- "tka"

#' Starting values for nlmer fixed effects
#'
#' @param x rxUiGet list(ui)
#' @param ... additional arguments (currently ignored)
#' @return Named numeric vector: iniDf estimates for thetas, 0 for nonMuEtas
#' @export
rxUiGet.nlmerStart <- function(x, ...) {
  .ui <- x[[1]]
  .pars <- rxUiGet.nlmerThetas(x, ...)
  .nonMu <- .ui$nonMuEtas
  .iniDf <- .ui$iniDf
  setNames(
    vapply(.pars, function(.p) {
      if (.p %in% .nonMu) return(0.0)
      .w <- which(.iniDf$name == .p)
      if (length(.w) == 1L) return(.iniDf$est[.w])
      0.0
    }, double(1), USE.NAMES = FALSE),
    .pars
  )
}
attr(rxUiGet.nlmerStart, "rstudio") <- c(tka = 0.45)

# -----------------------------------------------------------------------
# Formula builders
# -----------------------------------------------------------------------

#' Build the random-effects part (part 3) of the nlmer 3-part formula
#'
#' @param x rxUiGet list(ui)
#' @param ... additional arguments (currently ignored)
#' @return Formula: param1 + param2 + ... + (rand1 + rand2 | ID)
#' @export
rxUiGet.nlmerRandomFormula <- function(x, ...) {
  .ui <- x[[1]]
  .allPars <- rxUiGet.nlmerThetas(x, ...)
  .muRef <- .ui$muRefDataFrame$theta
  .nonMu <- .ui$nonMuEtas
  .randPars <- c(.muRef, .nonMu)
  .randPars <- .randPars[.randPars %in% .allPars]

  .fixed <- paste(.allPars, collapse = " + ")
  if (length(.randPars) == 0L) {
    as.formula(paste("~", .fixed))
  } else {
    .rand <- paste(.randPars, collapse = " + ")
    as.formula(paste("~", .fixed, "+ (", .rand, "| ID)"))
  }
}
attr(rxUiGet.nlmerRandomFormula, "rstudio") <- tka + tcl + tv ~ 1

#' Build the full 3-part nlmer formula
#'
#' @param x rxUiGet list(ui)
#' @param ... additional arguments (currently ignored)
#' @return 3-part formula for lme4::nlmer
#' @export
rxUiGet.nlmerFormula <- function(x, ...) {
  .pars <- rxUiGet.nlmerThetas(x, ...)
  .part2Args <- paste(c("TIME", .pars), collapse = ", ")
  .part2 <- paste0("babelmixr2::.nlmerNlmerFun(", .part2Args, ")")
  .part3 <- deparse(rxUiGet.nlmerRandomFormula(x, ...)[[2]])
  as.formula(paste("DV ~", .part2, "~", .part3))
}
attr(rxUiGet.nlmerFormula, "rstudio") <-
  DV ~ babelmixr2::.nlmerNlmerFun(TIME, tka, tcl, tv) ~ tka + tcl + tv + (tka | ID)

# -----------------------------------------------------------------------
# Sensitivity model: reuse nlm/nls sensitivity infrastructure
#
# We build an nlmer-specific model in THETA[i] notation so that the
# standard rxExpandFEta_ / nlmHdTheta machinery can compute
# d(rx_pred_)/d(THETA[i]) analytically.
#
# Prefix lines:  tka <- THETA[1]; tcl <- THETA[2]; ...
# Base model:    saem-dropped (mu-ref etas removed from body)
# Sensitivity:   rx__sens_rx_pred__BY_THETA_i___ for each param
# -----------------------------------------------------------------------

#' Build `THETA[i]` prefix substitution lines for nlmer
#'
#' @param ui rxode2 ui object
#' @return List of quoted assignment expressions
#' @noRd
.uiGetNlmerTheta <- function(ui) {
  .params <- ui$saemParamsToEstimate
  .env <- new.env(parent = emptyenv())
  .env$i <- 0L
  lapply(.params, function(.p) {
    .env$i <- .env$i + 1L
    eval(parse(text = paste0("quote(", .p, " <- THETA[", .env$i, "])")))
  })
}

#' nlmer base model (saem-dropped etas, `THETA[i]` prefix, nlme error lines)
#'
#' @param x rxUiGet list(ui)
#' @param ... additional arguments (currently ignored)
#' @export
rxUiGet.nlmerModel0 <- function(x, ...) {
  .f <- x[[1]]
  rxode2::rxCombineErrorLines(
    .f,
    errLines = nlmixr2est::rxGetDistributionNlmeLines(.f),
    prefixLines = .uiGetNlmerTheta(.f),
    paramsLine = NA,
    modelVars = TRUE,
    cmtLines = FALSE,
    dvidLine = FALSE,
    lstExpr = nlmixr2est::.saemDropMuRefFromModel(.f, noCovs = TRUE)
  )
}
attr(rxUiGet.nlmerModel0, "rstudio") <- quote(rxModelVars({}))

#' Prune the nlmer model for symengine loading
#'
#' @param x rxUiGet list(ui)
#' @param ... additional arguments (currently ignored)
#' @noRd
.nlmerPrune <- function(x) {
  .x <- x[[1]]
  .x <- .x$nlmerModel0[[-1]]
  .env <- new.env(parent = emptyenv())
  .env$.if <- NULL
  .env$.def1 <- NULL
  rxode2::.malert("pruning branches ({.code if}/{.code else}) of nlmer model...")
  .ret <- rxode2::.rxPrune(.x, envir = .env,
                           strAssign = rxode2::rxModelVars(x[[1]])$strAssign)
  .mv <- rxode2::rxModelVars(.ret)
  if (rxode2::.rxIsLinCmt() == 1L) {
    .vars <- c(.mv$params, .mv$lhs, .mv$slhs)
    .mv <- rxode2::.rxLinCmtGen(length(.mv$state), .vars)
  }
  rxode2::.msuccess("done")
  rxode2::rxNorm(.mv)
}

#' Load pruned nlmer model into symengine (no sensitivity promotion)
#'
#' @param x rxUiGet list(ui)
#' @param ... additional arguments (currently ignored)
#' @export
rxUiGet.loadPruneNlmer <- function(x, ...) {
  nlmixr2est::.loadSymengine(.nlmerPrune(x), promoteLinSens = FALSE)
}
attr(rxUiGet.loadPruneNlmer, "rstudio") <- emptyenv()

#' Load pruned nlmer model into symengine (with sensitivity promotion)
#'
#' @param x rxUiGet list(ui)
#' @param ... additional arguments (currently ignored)
#' @export
rxUiGet.loadPruneNlmerSens <- function(x, ...) {
  nlmixr2est::.loadSymengine(.nlmerPrune(x), promoteLinSens = TRUE)
}
attr(rxUiGet.loadPruneNlmerSens, "rstudio") <- emptyenv()

#' Compute d(state)/d(`THETA[i]`) sensitivities for nlmer
#'
#' @param x rxUiGet list(ui)
#' @param ... additional arguments (currently ignored)
#' @rdname rxUiGet.nlmerThetaSens
#' @export
rxUiGet.nlmerThetaS <- function(x, ...) {
  .s <- rxUiGet.loadPruneNlmerSens(x, ...)
  nlmixr2est::.sensEtaOrTheta(.s, theta = TRUE)
}
attr(rxUiGet.nlmerThetaS, "rstudio") <- emptyenv()

#' Compute d(f)/d(`THETA[i]`) for all nlmer parameters
#'
#' @param x rxUiGet list(ui)
#' @param ... additional arguments (currently ignored)
#' @export
rxUiGet.nlmerHdTheta <- function(x, ...) {
  .s <- rxUiGet.nlmerThetaS(x)
  .stateVars <- nlmixr2est::rxode2stateOde(.s)
  .predMinusDv <- rxode2::rxGetControl(x[[1]], "predMinusDv", TRUE)
  .grd <- rxode2::rxExpandFEta_(
    .stateVars, .s$..maxTheta,
    ifelse(.predMinusDv, 1L, 2L),
    isTheta = TRUE
  )
  if (rxode2::.useUtf()) {
    rxode2::.malert("calculate \u2202(f)/\u2202(\u03b8)")
  } else {
    rxode2::.malert("calculate d(f)/d(theta)")
  }
  rxode2::rxProgress(dim(.grd)[1])
  on.exit(rxode2::rxProgressAbort())
  .any.zero <- FALSE
  .all.zero <- TRUE
  .ret <- apply(.grd, 1, function(x) {
    .l <- x["calc"]
    .l <- eval(parse(text = .l))
    .ret <- paste0(x["dfe"], "=", rxode2::rxFromSE(.l))
    .zErr <- suppressWarnings(try(as.numeric(get(x["dfe"], .s)), silent = TRUE))
    if (identical(.zErr, 0)) {
      .any.zero <<- TRUE
    } else if (.all.zero) {
      .all.zero <<- FALSE
    }
    rxode2::rxTick()
    .ret
  })
  if (.all.zero) {
    stop("none of the predictions depend on 'THETA'", call. = FALSE)
  }
  if (.any.zero) {
    warning("some of the predictions do not depend on 'THETA'", call. = FALSE)
  }
  .s$..HdTheta <- .ret
  .s$..pred.minus.dv <- .predMinusDv
  rxode2::rxProgressStop()
  .s
}
attr(rxUiGet.nlmerHdTheta, "rstudio") <- emptyenv()

#' params() string for the nlmer `THETA[i]`-indexed model
#'
#' @param x rxUiGet list(ui)
#' @param ... additional arguments (currently ignored)
#' @export
rxUiGet.nlmerParams <- function(x, ...) {
  .ui <- x[[1]]
  .n <- length(.ui$saemParamsToEstimate)
  paste0("params(", paste(c(paste0("THETA[", seq_len(.n), "]"), "DV"),
                          collapse = ", "), ")")
}
attr(rxUiGet.nlmerParams, "rstudio") <- "params(THETA[1], DV)"

#' Finalize nlmer rxode2 model strings from symengine environment
#'
#' @param .s symengine environment (from nlmerHdTheta)
#' @param sum.prod logical
#' @param optExpression logical
#' @noRd
.rxFinalizeNlmer <- function(.s, sum.prod = FALSE, optExpression = TRUE) {
  .prd <- get("rx_pred_", envir = .s)
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.prd))
  .yj <- paste(get("rx_yj_", envir = .s))
  .yj <- paste0("rx_yj_~", rxode2::rxFromSE(.yj))
  .lambda <- paste(get("rx_lambda_", envir = .s))
  .lambda <- paste0("rx_lambda_~", rxode2::rxFromSE(.lambda))
  .hi <- paste(get("rx_hi_", envir = .s))
  .hi <- paste0("rx_hi_~", rxode2::rxFromSE(.hi))
  .low <- paste(get("rx_low_", envir = .s))
  .low <- paste0("rx_low_~", rxode2::rxFromSE(.low))
  .ddt <- .s$..ddt
  if (is.null(.ddt)) .ddt <- character(0)
  .sens <- .s$..sens
  if (is.null(.sens)) .sens <- character(0)

  .s$..nlmerS <- paste(c(
    .s$params,
    .s$..stateInfo["state"],
    .ddt,
    .sens,
    .yj,
    .lambda,
    .hi,
    .low,
    .prd,
    .s$..HdTheta,
    .s$..stateInfo["statef"],
    .s$..stateInfo["dvid"],
    ""
  ), collapse = "\n")

  .lhs0 <- .s$..lhs0
  if (is.null(.lhs0)) .lhs0 <- ""
  .s$..pred.nolhs <- paste(c(
    .s$params,
    .s$..stateInfo["state"],
    .lhs0,
    .ddt,
    .yj,
    .lambda,
    .hi,
    .low,
    .prd,
    .s$..stateInfo["statef"],
    .s$..stateInfo["dvid"],
    ""
  ), collapse = "\n")

  if (sum.prod) {
    rxode2::.malert("stabilizing round off errors in nlmer gradient problem...")
    .s$..nlmerS <- rxode2::rxSumProdModel(.s$..nlmerS)
    rxode2::.msuccess("done")
    rxode2::.malert("stabilizing round off errors in nlmer pred-only problem...")
    .s$..pred.nolhs <- rxode2::rxSumProdModel(.s$..pred.nolhs)
    rxode2::.msuccess("done")
  }
  if (optExpression) {
    .s$..nlmerS <- rxode2::rxOptExpr(.s$..nlmerS, "nlmer gradient")
    .s$..pred.nolhs <- rxode2::rxOptExpr(.s$..pred.nolhs, "nlmer pred-only")
  }
}

#' Full nlmer symengine environment with sensitivities
#'
#' @param x rxUiGet list(ui)
#' @param ... additional arguments (currently ignored)
#' @export
rxUiGet.nlmerEnv <- function(x, ...) {
  .s <- rxUiGet.nlmerHdTheta(x, ...)
  .s$params <- rxUiGet.nlmerParams(x, ...)
  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  .rxFinalizeNlmer(.s, .sumProd, .optExpression)
  .s$..outer <- NULL
  # `eventTheta` flags the structural parameters that enter a dosing
  # expression (alag/F/rate/dur).  Under eventSens="jump" rxode2 injects the
  # analytic event jump into those sensitivity states, so the analytic
  # gradient is correct and the flag is left at 0 (analytic gradient is used).
  # Under "fd" the augmented sensitivity ODE misses the jump, so the parameter
  # is flagged and the nlm engine overrides its gradient column with a Shi2021
  # finite difference.
  if (exists("..maxTheta", .s)) {
    .eventTheta <- rep(0L, .s$..maxTheta)
  } else {
    .eventTheta <- integer(0)
  }
  .eventSens <- rxode2::rxGetControl(x[[1]], "eventSens", "jump")
  if (!identical(.eventSens, "jump")) {
    .reg <- "^THETA\\[([0-9]+)\\]$"
    for (.v in .s$..eventVars) {
      .vars <- as.character(get(.v, envir = .s))
      .vars <- rxode2::rxGetModel(paste0("rx_lhs=", rxode2::rxFromSE(.vars)))$params
      for (.v2 in .vars) {
        if (regexpr(.reg, .v2) != -1) {
          .num <- as.numeric(sub(.reg, "\\1", .v2))
          .eventTheta[.num] <- 1L
        }
      }
    }
  }
  .s$.eventTheta <- .eventTheta
  .s
}
attr(rxUiGet.nlmerEnv, "rstudio") <- emptyenv()

#' Compiled nlmer sensitivity model (prediction + gradients)
#'
#' Returns a list with:
#' - `thetaGrad`: compiled rxode2 model computing rx_pred_ and
#'   rx__sens_rx_pred__BY_THETA_i___ for each estimated parameter.  When
#'   `eventSens = "jump"` this model carries rxode2's analytic event-jump
#'   sensitivities (loaded/activated by [nlmixr2est::.nlmSetupEnv()]).
#' - `predOnly`: compiled rxode2 model for prediction only
#' - `eventTheta`: integer vector flagging `THETA[i]` that need Shi2021 finite
#'   differences for dosing parameters (all 0 under `eventSens = "jump"`)
#' - `paramNames`: character vector of parameter names in `THETA[i]` order
#'
#' The list follows the `modelInfo` contract of [nlmixr2est::.nlmSetupEnv()]
#' so the nlm C machinery can load it once and keep it resident.
#'
#' @param x rxUiGet list(ui)
#' @param ... additional arguments (currently ignored)
#' @export
rxUiGet.nlmerSensModel <- function(x, ...) {
  .s <- rxUiGet.nlmerEnv(x, ...)
  .eventSens <- rxode2::rxGetControl(x[[1]], "eventSens", "jump")
  list(
    thetaGrad  = rxode2::rxode2(.s$..nlmerS, eventSens = .eventSens),
    predOnly   = rxode2::rxode2(.s$..pred.nolhs),
    eventTheta = .s$.eventTheta,
    paramNames = rxUiGet.nlmerThetas(x, ...)
  )
}
attr(rxUiGet.nlmerSensModel, "rstudio") <- list()
