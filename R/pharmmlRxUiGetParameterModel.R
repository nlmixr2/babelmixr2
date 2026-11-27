#' Map an nlmixr2 mu-referenced transformation to a PharmML transformation
#'
#' This is the PharmML counterpart of Monolix's
#' `.mlxTranCurEvalToDistribution()`: the same `muRefCurEval$curEval` values,
#' rendered as PharmML `Transformation/@type` rather than as a Monolix
#' distribution name.  `add` (an untransformed parameter) has no
#' `Transformation` element at all.
#'
#' @param curEval Current evaluation from `ui$muRefCurEval`
#'
#' @return PharmML transformation type, or `NA_character_` when the parameter
#'   is untransformed
#'
#' @noRd
.pharmmlCurEvalToTransformation <- function(curEval) {
  .ret <- switch(ifelse(curEval %in% c("", "*", "**", "/", "^", "+", "-"),
                        "add", curEval),
                 exp = "log",
                 expit = "logit",
                 probitInv = "probit",
                 add = NA_character_,
                 "unsupported")
  if (identical(.ret, "unsupported")) {
    stop("PharmML translation of the transformation '", curEval,
         "' is not supported", call. = FALSE)
  }
  .ret
}

#' The transformation applied to a mu-referenced theta
#'
#' @param est theta name
#' @param ui rxode2 UI
#' @return PharmML transformation type or NA_character_
#' @noRd
.pharmmlThetaTransformation <- function(est, ui) {
  .cur <- ui$muRefCurEval
  .w <- which(.cur$parameter == est)
  if (length(.w) == 0L) return(NA_character_)
  if (length(.w) > 1L) {
    stop("duplicate parameter '", est, "' in `muRefCurEval`", call. = FALSE)
  }
  .pharmmlCurEvalToTransformation(.cur$curEval[.w])
}

#' The eta paired with a mu-referenced theta, if any
#'
#' @param est theta name
#' @param ui rxode2 UI
#' @return eta name, or NA_character_ when the theta has no random effect
#' @noRd
.pharmmlEtaFor <- function(est, ui) {
  .tab <- ui$muRefTable
  if (is.null(.tab)) return(NA_character_)
  .w <- which(.tab$theta == est)
  if (length(.w) != 1L) return(NA_character_)
  .tab$eta[.w]
}

#' Name of the population parameter holding an eta's variance
#'
#' @param eta eta name
#' @return population parameter name
#' @noRd
.pharmmlOmegaName <- function(eta) {
  paste0("omega_", eta)
}

#' Name of the population parameter holding an eta pair's covariance
#'
#' @param eta1 first eta name
#' @param eta2 second eta name
#' @return population parameter name
#' @noRd
.pharmmlCovName <- function(eta1, eta2) {
  paste0("cov_", eta1, "_", eta2)
}

#' A ct:Assign wrapping a single already-emitted child
#'
#' @param child emitted child node
#' @param indent indent depth
#' @return character(1)
#' @noRd
.pharmmlAssign <- function(child, indent = 0L) {
  .pmlNode("ct:Assign", children = child, indent = indent)
}

#' A ct:Assign wrapping a bare symbol reference
#'
#' @param sym symbol name
#' @param indent indent depth
#' @return character(1)
#' @noRd
.pharmmlAssignSym <- function(sym, indent = 0L) {
  .pharmmlAssign(.pmlNode("ct:SymbRef", attrs = c(symbIdRef = sym)), indent)
}

#' The RandomVariable block for one eta
#'
#' nlmixr2 stores eta variances, so the ProbOnto `Normal2` (mean/var)
#' parameterisation is used rather than `Normal1` (mean/stdev); that keeps the
#' written value identical to the `iniDf` estimate with no square root.
#'
#' @param eta eta name
#' @param indent indent depth
#' @return character(1)
#' @noRd
.pharmmlRandomVariable <- function(eta, indent = 0L) {
  .pmlNode(
    "mdef:RandomVariable",
    attrs = c(symbId = eta),
    children = c(
      .pmlNode("ct:VariabilityReference",
               children = .pmlNode("ct:SymbRef",
                                   attrs = c(blkIdRef = .pmlBlk[["variabilityParameter"]],
                                             symbIdRef = "id"))),
      .pmlNode("mdef:Distribution",
               children = .pmlNode(
                 "po:ProbOnto", attrs = c(name = "Normal2"),
                 children = c(
                   .pmlNode("po:Parameter", attrs = c(name = "mean"),
                            children = .pharmmlAssign(.pmlText("ct:Real", 0))),
                   .pmlNode("po:Parameter", attrs = c(name = "var"),
                            children = .pharmmlAssignSym(.pharmmlOmegaName(eta))))))),
    indent = indent)
}

#' The IndividualParameter block for one mu-referenced parameter
#'
#' @param var model variable name (the left-hand side in the nlmixr2 model)
#' @param est theta name
#' @param ui rxode2 UI
#' @param indent indent depth
#' @return character(1)
#' @noRd
.pharmmlIndividualParameter <- function(var, est, ui, indent = 0L) {
  .trans <- .pharmmlThetaTransformation(est, ui)
  .eta <- .pharmmlEtaFor(est, ui)
  .cov <- .pharmmlMuRefCovariates(est, ui)

  .children <- character(0)
  if (!is.na(.trans)) {
    .children <- c(.children,
                   .pmlNode("mdef:Transformation", attrs = c(type = .trans)))
  }
  if (length(.cov) == 0L) {
    .children <- c(.children,
                   .pmlNode("mdef:PopulationValue",
                            children = .pharmmlAssignSym(est)))
  } else {
    .children <- c(.children,
                   .pmlNode("mdef:LinearCovariate",
                            children = c(
                              .pmlNode("mdef:PopulationValue",
                                       children = .pharmmlAssignSym(est)),
                              .cov)))
  }
  if (!is.na(.eta)) {
    .children <- c(.children,
                   .pmlNode("mdef:RandomEffects",
                            children = .pharmmlAssignSym(.eta)))
  }

  .pmlNode("mdef:IndividualParameter", attrs = c(symbId = var),
           children = .pmlNode("mdef:StructuredModel", children = .children),
           indent = indent)
}

#' Covariate terms for one mu-referenced theta
#'
#' @param est theta name
#' @param ui rxode2 UI
#' @return character vector of emitted `mdef:Covariate` nodes, possibly empty
#' @noRd
.pharmmlMuRefCovariates <- function(est, ui) {
  .df <- ui$saemMuRefCovariateDataFrame
  if (is.null(.df) || nrow(.df) == 0L) return(character(0))
  .w <- which(.df$theta == est)
  if (length(.w) == 0L) return(character(0))
  vapply(.w, function(.i) {
    .pmlNode("mdef:Covariate",
             children = c(
               .pmlNode("ct:SymbRef",
                        attrs = c(blkIdRef = .pmlBlk[["covariate"]],
                                  symbIdRef = .df$covariate[.i])),
               .pmlNode("mdef:FixedEffect",
                        children = .pharmmlAssignSym(.df$covariateParameter[.i]))))
  }, character(1), USE.NAMES = FALSE)
}

#' The Correlation block for off-diagonal omega entries
#'
#' This is the PharmML counterpart of Monolix's `.mlxtranIndividualCor()`.
#'
#' @param ui rxode2 UI
#' @return character(1), or `character(0)` when the omega matrix is diagonal
#' @noRd
.pharmmlCorrelation <- function(ui) {
  .iniDf <- ui$iniDf
  .eta <- .iniDf[!is.na(.iniDf$neta1), , drop = FALSE]
  if (nrow(.eta) == 0L) return(character(0))
  .off <- .eta[.eta$neta1 != .eta$neta2, , drop = FALSE]
  if (nrow(.off) == 0L) return(character(0))

  .nameOf <- function(n) {
    .w <- which(.eta$neta1 == n & .eta$neta2 == n)
    if (length(.w) != 1L) {
      stop("cannot find the eta with index ", n, call. = FALSE)
    }
    .eta$name[.w]
  }

  .pairs <- lapply(seq_len(nrow(.off)), function(.i) {
    # Normalise the ordering: iniDf may list the pair either way round, and the
    # covariance parameter name has to be stable regardless.
    .lo <- min(.off$neta1[.i], .off$neta2[.i])
    .hi <- max(.off$neta1[.i], .off$neta2[.i])
    .e1 <- .nameOf(.lo)
    .e2 <- .nameOf(.hi)
    .pmlNode("mdef:Pairwise",
             children = c(
               .pmlNode("mdef:RandomVariable1",
                        children = .pharmmlAssignSym(.e1)),
               .pmlNode("mdef:RandomVariable2",
                        children = .pharmmlAssignSym(.e2)),
               .pmlNode("mdef:Covariance",
                        children = .pharmmlAssignSym(.pharmmlCovName(.e1, .e2)))))
  })

  .pmlNode("mdef:Correlation",
           children = c(
             .pmlNode("ct:VariabilityReference",
                      children = .pmlNode("ct:SymbRef",
                                          attrs = c(blkIdRef = .pmlBlk[["variabilityParameter"]],
                                                    symbIdRef = "id"))),
             unlist(.pairs, use.names = FALSE)))
}

#' PharmML ParameterModel block
#'
#' Built from exactly the information the Monolix `[INDIVIDUAL]` writer uses:
#' `$getSplitMuModel`, `$muRefTable`, `$muRefCurEval` and
#' `$saemMuRefCovariateDataFrame`.
#'
#' @param ui rxode2 UI
#'
#' @param indent Indent depth
#'
#' @return character(1) holding the `ParameterModel` block
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.pharmmlParameterModel <- function(ui, indent = 0L) {
  .split <- ui$getSplitMuModel
  .muRef <- c(.split$pureMuRef, .split$taintMuRef)

  .popParams <- character(0)
  .randomVars <- character(0)
  .indiv <- character(0)

  for (.i in seq_along(.muRef)) {
    .est <- names(.muRef)[.i]
    .var <- setNames(.muRef[.i], NULL)
    .popParams <- c(.popParams,
                    .pmlNode("mdef:PopulationParameter", attrs = c(symbId = .est)))
    .eta <- .pharmmlEtaFor(.est, ui)
    if (!is.na(.eta)) {
      .popParams <- c(.popParams,
                      .pmlNode("mdef:PopulationParameter",
                               attrs = c(symbId = .pharmmlOmegaName(.eta))))
      .randomVars <- c(.randomVars, .pharmmlRandomVariable(.eta))
    }
    for (.cp in .pharmmlMuRefCovariateParams(.est, ui)) {
      .popParams <- c(.popParams,
                      .pmlNode("mdef:PopulationParameter", attrs = c(symbId = .cp)))
    }
    .indiv <- c(.indiv, .pharmmlIndividualParameter(.var, .est, ui))
  }

  .corr <- .pharmmlCorrelation(ui)
  for (.cn in .pharmmlCovarianceNames(ui)) {
    .popParams <- c(.popParams,
                    .pmlNode("mdef:PopulationParameter", attrs = c(symbId = .cn)))
  }

  .pmlNode("mdef:ParameterModel",
           attrs = c(blkId = .pmlBlk[["parameter"]]),
           children = c(.popParams, .randomVars, .indiv, .corr),
           indent = indent)
}

#' Covariate coefficient parameter names for one theta
#'
#' @param est theta name
#' @param ui rxode2 UI
#' @return character vector, possibly empty
#' @noRd
.pharmmlMuRefCovariateParams <- function(est, ui) {
  .df <- ui$saemMuRefCovariateDataFrame
  if (is.null(.df) || nrow(.df) == 0L) return(character(0))
  .df$covariateParameter[.df$theta == est]
}

#' Covariance population-parameter names for the off-diagonal omegas
#'
#' @param ui rxode2 UI
#' @return character vector, possibly empty
#' @noRd
.pharmmlCovarianceNames <- function(ui) {
  .iniDf <- ui$iniDf
  .eta <- .iniDf[!is.na(.iniDf$neta1), , drop = FALSE]
  if (nrow(.eta) == 0L) return(character(0))
  .off <- .eta[.eta$neta1 != .eta$neta2, , drop = FALSE]
  if (nrow(.off) == 0L) return(character(0))
  .nameOf <- function(n) .eta$name[which(.eta$neta1 == n & .eta$neta2 == n)]
  vapply(seq_len(nrow(.off)), function(.i) {
    .pharmmlCovName(.nameOf(min(.off$neta1[.i], .off$neta2[.i])),
                    .nameOf(max(.off$neta1[.i], .off$neta2[.i])))
  }, character(1), USE.NAMES = FALSE)
}

#' @export
rxUiGet.pharmmlParameterModel <- function(x, ...) {
  .pharmmlParameterModel(x[[1]])
}
attr(rxUiGet.pharmmlParameterModel, "rstudio") <- "character"
