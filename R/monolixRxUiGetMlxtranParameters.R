.getOmegaR <- function(ui) {
  .cov <- ui$omega
  .sd2 <- sqrt(diag(.cov))
  .cor <- stats::cov2cor(.cov)
  dimnames(.cor) <- dimnames(.cov)
  diag(.cor) <- .sd2
  .cor
}

#' @export
rxUiGet.mlxtranParameter <- function(x, ...) {
  .ui <- x[[1]]
  .r <- .getOmegaR(.ui)
  .split <- rxUiGet.getSplitMuModel(x, ...)
  .muRef <- c(.split$pureMuRef, .split$taintMuRef)
  .iniDf <- .ui$iniDf
  .covDataFrame <- .ui$saemMuRefCovariateDataFrame

  paste0("\n\n<PARAMETER>\n",paste(vapply(seq_along(.iniDf$name), function(i) {
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
        .est <- .cur$est
        return(paste0(.par, "={value=", .est, ", method=",
                      ifelse(.cur$fix, "FIXED", "MLE"), "}"))
      } else {
        .par <- eval(str2lang(paste0("rxToMonolix(", .cur$name, ")")))
        .est <- .cur$est
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

