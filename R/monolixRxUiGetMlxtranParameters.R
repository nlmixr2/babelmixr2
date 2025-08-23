#' Get the omega correlation matroix
#'
#'
#' @param ui rxode2 ui object
#' @return correlation matrix for the omega parameters
#' @noRd
#' @author Matthew L. Fidler
.getOmegaR <- function(ui) {
  .cov <- ui$omega
  .sd2 <- sqrt(diag(.cov))
  .cor <- stats::cov2cor(.cov)
  dimnames(.cor) <- dimnames(.cov)
  diag(.cor) <- .sd2
  .cor
}
#' This gets the monolix parameter ini based on monolix's naming structure
#'
#' This really is not where the saem starts, but it is the "standard" parameter
#'
#' @param est estimate value from nlmixr2
#' @param name name of the parameter from nlmixr2
#' @param curEval current evaluation data-frame to understand the transformation.
#' @param ui rxode2 ui object to get if the parmeter is an variance error term or not.
#' @return estimate as put into the monolix .mlxtran file
#' @noRd
#' @author Matthew L. Fidler
.getNonMonolixParameterIni <- function(est, name, curEval, ui) {
  .ui <- ui
  .w <- which(curEval$parameter == name)
  if (length(.w) == 1) {
    .ce <- paste(curEval$curEval[.w])
    .low <- curEval$low[.w]
    if (is.na(.low)) .low <- 0
    .high <- curEval$hi[.w]
    if (is.na(.high)) .high <- 0
    return(switch(.ce,
                  exp=exp(est),
                  expit=rxode2::expit(est, .low, .high),
                  probitInv=rxode2::probitInv(est, .low, .high),
                  est))
  }
  .predDf <- .ui$predDf
  if (any(names(.predDf) == "variance")) {
    # now lets check if the term is from a variance error term
    .iniDf <- .ui$iniDf
    .w <- which(.iniDf$name == name)
    if (length(.w) == 1) {
      .cond <- .iniDf$condition[.w]
      .w2 <- which(.predDf$cond == .cond)
      if (length(.w2) == 1) {
        if (.predDf$variance[.w2] &&
              .iniDf$err[.w] %in% c("add", "prop", "propT",
                                   "pow", "powT", "logn",
                                   "dlogn", "lnorm",
                                   "dlnorm", "logitNorm",
                                   "probitNorm")) {
          # monolix estimates sd instead of variance (like nlmixr2 saem)
          return(sqrt(est))
        }
      }
    }
  }
  return(est)
}

#' @export
rxUiGet.mlxtranParameter <- function(x, ...) {
  .ui <- x[[1]]
  .r <- .getOmegaR(.ui)
  .split <- .ui$getSplitMuModel
  .muRef <- c(.split$pureMuRef, .split$taintMuRef)
  .iniDf <- .ui$iniDf
  .covDataFrame <- .ui$saemMuRefCovariateDataFrame
  .curEval <- .ui$muRefCurEval
  paste0("<PARAMETER>\n",paste(vapply(seq_along(.iniDf$name), function(i) {
    .cur <- .iniDf[i, ]
    if (is.na(.cur$neta1)) {
      if (is.na(.cur$err)) {
        .w <- which(.covDataFrame$covariateParameter == .cur$name)
        if (length(.w) == 1) {
          .par <- paste0("beta_", .muRef[.covDataFrame$theta[.w]], "_",
                         .covDataFrame$covariate[.w])
        } else {
          .monolixName <- .muRef[.cur$name]
          if (is.na(.monolixName)) stop("babelmixr2 can't figure out '", .cur$name, "' for monolix",
                                        call.=FALSE)
          .par <- paste0(.monolixName, "_pop")
        }
        .est <- .getNonMonolixParameterIni(.cur$est, .cur$name, .curEval, .ui)
        return(paste0(.par, "={value=", .est, ", method=",
                      ifelse(.cur$fix, "FIXED", "MLE"), "}"))
      } else {
        .par <- eval(str2lang(paste0("rxToMonolix(", .cur$name, ", ui=.ui)")))
        .est <- .getNonMonolixParameterIni(.cur$est, .cur$name, .curEval, .ui)
        return(paste0(.par, "={value=", .est, ", method=",
                      ifelse(.cur$fix, "FIXED", "MLE"), "}"))
      }
    } else if (.cur$neta1 == .cur$neta2) {
      .par <- paste0("omega_", .mlxtranGetIndividualMuRefEtaMonolixName(.ui, .cur$neta1, .muRef))
      .est <- .r[.cur$name, .cur$name]
      return(paste0(.par, "={value=", .est, ", method=",
                    ifelse(.cur$fix, "FIXED", "MLE"), "}"))
    } else {
      .par <- paste0("corr_", .mlxtranGetIndividualMuRefEtaMonolixName(.ui, .cur$neta1, .muRef),"_",
                     .mlxtranGetIndividualMuRefEtaMonolixName(.ui, .cur$neta2, .muRef))
      .n1 <- .iniDf[which(.iniDf$neta1 == .cur$neta1 & .iniDf$neta2 == .cur$neta1), "name"]
      .n2 <- .iniDf[which(.iniDf$neta1 == .cur$neta2 & .iniDf$neta2 == .cur$neta2), "name"]
      .est <- .r[.n1, .n2]
      return(paste0(.par, "={value=", .est, ", method=",
                    ifelse(.cur$fix, "FIXED", "MLE"), "}"))
    }
  }, character(1), USE.NAMES=FALSE), collapse="\n"))
}

attr(rxUiGet.mlxtranParameter, "rstudio") <- "mlxtranParameter"
