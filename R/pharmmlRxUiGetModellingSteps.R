#' The parameters a PharmML EstimationStep should list
#'
#' These are exactly the population parameters the ParameterModel and
#' ObservationModel declare: thetas, covariate coefficients, residual-error
#' parameters, and the omega variances/covariances.
#'
#' Residual-error parameters are declared in the ObservationModel of their
#' endpoint (`om1`, `om2`, ...), everything else in the ParameterModel, so each
#' row also carries the block its `SymbRef` has to name.
#'
#' @param ui rxode2 UI
#' @return data frame with `name`, `blk`, `est`, `lower`, `upper` and `fix`
#' @noRd
.pharmmlEstimatedParameters <- function(ui) {
  .iniDf <- ui$iniDf
  .theta <- .iniDf[!is.na(.iniDf$ntheta), , drop = FALSE]
  .blk <- rep(.pmlBlk[["parameter"]], nrow(.theta))
  .err <- !is.na(.theta$err)
  .blk[.err] <- paste0(
    "om",
    match(.theta$condition[.err], paste(ui$predDf$cond))
  )
  .ret <- data.frame(
    name = .theta$name,
    blk = .blk,
    est = .theta$est,
    lower = .theta$lower,
    upper = .theta$upper,
    fix = .theta$fix,
    stringsAsFactors = FALSE
  )

  .eta <- .iniDf[!is.na(.iniDf$neta1), , drop = FALSE]
  if (nrow(.eta) > 0L) {
    .diag <- .eta[.eta$neta1 == .eta$neta2, , drop = FALSE]
    .ret <- rbind(
      .ret,
      data.frame(
        name = .pharmmlOmegaName(.diag$name),
        blk = .pmlBlk[["parameter"]],
        est = .diag$est,
        lower = 0,
        upper = Inf,
        fix = .diag$fix,
        stringsAsFactors = FALSE
      )
    )
    .off <- .eta[.eta$neta1 != .eta$neta2, , drop = FALSE]
    if (nrow(.off) > 0L) {
      .ret <- rbind(
        .ret,
        data.frame(
          name = .pharmmlCovarianceNames(ui),
          blk = .pmlBlk[["parameter"]],
          est = .off$est,
          lower = -Inf,
          upper = Inf,
          fix = .off$fix,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  .ret
}

#' One ParameterEstimation entry
#'
#' @param row a row of `.pharmmlEstimatedParameters()`
#' @return character(1)
#' @noRd
.pharmmlParameterEstimation <- function(row) {
  .children <- .pmlNode(
    "ct:SymbRef",
    attrs = c(blkIdRef = row$blk, symbIdRef = row$name)
  )
  .init <- .pmlNode(
    "mstep:InitialEstimate",
    attrs = if (isTRUE(row$fix)) c(fixed = "true") else NULL,
    children = .pmlText("ct:Real", row$est)
  )
  .bounds <- character(0)
  if (is.finite(row$lower)) {
    .bounds <- c(
      .bounds,
      .pmlNode("mstep:LowerBound", children = .pmlText("ct:Real", row$lower))
    )
  }
  if (is.finite(row$upper)) {
    .bounds <- c(
      .bounds,
      .pmlNode("mstep:UpperBound", children = .pmlText("ct:Real", row$upper))
    )
  }
  # The schema's sequence is SymbRef, InitialEstimate, LowerBound, UpperBound.
  .pmlNode("mstep:ParameterEstimation", children = c(.children, .init, .bounds))
}

#' PharmML ModellingSteps
#'
#' A single `EstimationStep` describing the estimation the model is set up for.
#' No algorithm is named: `babelmixr2` writes the model, and which estimator a
#' receiving tool uses is that tool's business.
#'
#' @param ui rxode2 UI
#'
#' @param indent Indent depth
#'
#' @return character(1) holding the `ModellingSteps` element
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.pharmmlModellingSteps <- function(ui, indent = 0L) {
  .par <- .pharmmlEstimatedParameters(ui)
  .est <- vapply(
    seq_len(nrow(.par)),
    function(.i) {
      .pharmmlParameterEstimation(.par[.i, ])
    },
    character(1),
    USE.NAMES = FALSE
  )

  .step <- .pmlNode(
    "mstep:EstimationStep",
    attrs = c(oid = "estStep"),
    children = c(
      .pmlNode(
        "mstep:ExternalDataSetReference",
        children = .pmlNode("ct:OidRef", attrs = c(oidRef = "nmOid"))
      ),
      .pmlNode("mstep:ParametersToEstimate", children = .est),
      .pmlNode("mstep:Operation", attrs = c(order = "1", opType = "estPop"))
    )
  )

  .pmlNode(
    "mstep:ModellingSteps",
    children = c(
      .step,
      .pmlNode(
        "mstep:StepDependencies",
        children = .pmlNode(
          "mstep:Step",
          children = .pmlNode("ct:OidRef", attrs = c(oidRef = "estStep"))
        )
      )
    ),
    indent = indent
  )
}

#' @export
rxUiGet.pharmmlModellingSteps <- function(x, ...) {
  .pharmmlModellingSteps(x[[1]])
}
attr(rxUiGet.pharmmlModellingSteps, "rstudio") <- "character"
