#' @export
rxUiGet.monolixOutputVersion <- function(x, ...) {
  .exportPath <- rxUiGet.monolixExportPath(x, ...)
  .summary <- file.path(.exportPath, "summary.txt")
  if (file.exists(.summary)) {
    .lines <- readLines(.summary, n=5)
    .w <- which(regexpr(".*[vV]ersion *: *[^ ]*.*", .lines) != -1)
    if (length(.w) == 1) {
      .line <- .lines[.w]
      # 2019 is 5.1.1
      return(sub(".*[vV]ersion *: *([^ ]*).*", "\\1", .line))
    }
  }
  NULL
}


#' @export
rxUiGet.monolixPopulationParameters <- function(x, ...) {
  .exportPath <- rxUiGet.monolixExportPath(x, ...)
  .popParameters <- file.path(.exportPath, "populationParameters.txt")
  if (file.exists(.popParameters)) {
    return(read.csv(.popParameters))
  }
  NULL
}

#' @export
rxUiGet.monolixOmega <- function(x, ...) {
  .ui <- x[[1]]
  .pop <- rxUiGet.monolixPopulationParameters(x, ...)
  if (is.null(.pop)) return(NULL)
  .iniDf <- .ui$iniDf
  .split <- rxUiGet.getSplitMuModel(x, ...)
  .muRef <- c(.split$pureMuRef, .split$taintMuRef)
  .muRefCurEval <- .ui$muRefCurEval
  .eta <- .iniDf[is.na(.iniDf$ntheta), ]
  .r <- .getOmegaR(.ui)
  for (.i in seq_along(.eta$name)) {
    .n1 <- .eta$neta1[.i]
    .n2 <- .eta$neta2[.i]
    if (.n1 == .n2) {
      .par <- paste0("omega_", .mlxtranGetIndividualMuRefEtaMonolixName(.ui, .n1, .muRef))
      .w <- which(.pop$parameter == .par)
      if (length(.w) == 1){
        .par <- .pop$value[.w]
      } else {
        .par <- NA_real_
      }
      .r[.n1, .n2] <- .r[.n2, .n1] <- .par
    } else {
      .par <- paste0("corr_", .mlxtranGetIndividualMuRefEtaMonolixName(.ui, .n1, .muRef),"_",
                     .mlxtranGetIndividualMuRefEtaMonolixName(.ui, .n2, .muRef))
      .w <- which(.pop$parameter == .par)
      if (length(.w) == 1){
        .par <- .pop$value[.w]
      } else {
        .par <- NA_real_
      }
      .r[.n1, .n2] <- .r[.n2, .n1] <- .par
    }
  }
  .sd <- diag(.r)
  diag(.r) <- 1
  .d <- diag(.sd)
  .omega <- .d %*% .r %*% .d
  .omega
}

.monolixGetPopParValue <- function(name, muRefCurEval, muRef, covDataFrame, pop) {
  .w <- which(covDataFrame$covariateParameter == name)
  if (length(.w) == 1) {
    .par <- paste0("beta_", .muRef[.covDataFrame$theta[.w]], "_",
                   .covDataFrame$covariate[.w])
    .w <- which(pop$parameter == .par)
    if (length(.w) != 1) return(NA_real_)
    return(pop$value[.w])
  }
  .w <- which(muRefCurEval$parameter == name)
  if (length(.w) == 1) {
    .curEval <- as.character(muRefCurEval$curEval[.w])
    .low <- muRefCurEval$low[.w]
    .hi <- muRefCurEval$hi[.w]
    if (is.na(.low)) .low <- 0
    if (is.na(.hi)) .hi <- 1
  } else {
    .curEval <- ""
    .low <- 0
    .hi <- 1
  }
  .par <- paste0(muRef[name], "_pop")
  .w <- which(pop$parameter == .par)
  if (length(.w) != 1) return(NA_real_)
  .par <- pop$value[.w]
  return(switch(.curEval,
                exp=log(.par),
                expit=rxode2::logit(.par, .low, .hi),
                probitInv=rxode2::probit(.par, .low, .hi),
                .par))
}

#' @export
rxUiGet.monolixNewIniDf <- function(x, ...) {
  .ui <- x[[1]]
  .pop <- rxUiGet.monolixPopulationParameters(x, ...)
  if (is.null(.pop)) return(NULL)
  .omega <- rxUiGet.monolixOmega(x, ...)
  .iniDf <- .ui$iniDf
  .muRefCurEval <- .ui$muRefCurEval
  .split <- rxUiGet.getSplitMuModel(x, ...)
  .muRef <- c(.split$pureMuRef, .split$taintMuRef)
  .covDataFrame <- .ui$saemMuRefCovariateDataFrame
  .iniDf$est <- vapply(seq_along(.iniDf$name),
                  function(i) {
                    .n <- .iniDf$name[i]
                    if (is.na(.iniDf$ntheta[i])) {
                      return(.omega[.iniDf$neta1[i], .iniDf$neta2[i]])
                    }
                    .isPop <- is.na(.iniDf$err[i])
                    if (.isPop) {
                      return(.monolixGetPopParValue(.n, .muRefCurEval, .muRef, .covDataFrame, .pop))
                    }
                    .isErr <- !is.na(.iniDf$ntheta[i]) & !is.na(.iniDf$err[i])
                    if (.isErr) {
                      .par <- eval(str2lang(paste0("rxToMonolix(", .n, ", ui=.ui)")))
                      .w <- which(.pop$parameter == .par)
                      if (length(.w) != 1) return(NA_real_)
                      return(.pop$value[.w])
                    }
                    NA_real_
                  }, double(1), USE.NAMES=FALSE)
  .iniDf
}

#' @export
rxUiGet.monolixIndividualParameters <- function(x, ...) {
  .exportPath <- rxUiGet.monolixExportPath(x, ...)
  .individualParameters <- file.path(.exportPath, "IndividualParameters", "estimatedRandomEffects.txt")
  if (file.exists(.individualParameters)) {
    return(read.csv(.individualParameters))
  }
  NULL
}

 #' @export
xUiGet.monolixCovarianceEstimatesLin <- function(x, ...) {
  .exportPath <- rxUiGet.monolixExportPath(x, ...)
  .covLin <- file.path(.exportPath, "FisherInformation", "covarianceEstimatesLin.txt")
  if (file.exists(.covLin)) {
    .c <- read.csv(.covLin, header=FALSE)
    .n <- .c[, 1]
    .c <- as.matrix(.c[, -1])
    dimnames(.c) <- list(.n, .n)
    return(.c)
  }
  NULL
}


 #' @export
rxUiGet.monolixCovarianceEstimatesSA <- function(x, ...) {
  .exportPath <- rxUiGet.monolixExportPath(x, ...)
  .covSA <- file.path(.exportPath, "FisherInformation", "covarianceEstimatesSA.txt")
  if (file.exists(.covSA)) {
    .c <- read.csv(.covSA, header=FALSE)
    .n <- .c[, 1]
    .c <- as.matrix(.c[, -1])
    dimnames(.c) <- list(.n, .n)
    return(.c)
  }
  NULL
}

.monolixCovarianceNeedsConversion <- function(x, sa) {
  .ver <- rxUiGet.monolixOutputVersion(x)
  if (!is.character(.ver)) return(TRUE)
  .reg <- ".*([0-9][0-9][0-9][0-9]).*"
  if (regexpr(.reg, .ver) != -1) {
    .num <- as.numeric(sub(.reg, "\\1", .ver))
    if (sa) {
      return(.num < 2020)
    } else {
      return(.num < 2021)
    }
  }
  TRUE
}

#' @export
rxUiGet.monolixCovariance <- function(x, ...) {
  .cov <- rxUiGet.monolixCovarianceEstimatesSA(x, ...)
  .split <- rxUiGet.getSplitMuModel(x, ...)
  .muRef <- c(.split$pureMuRef, .split$taintMuRef)
  .sa <- TRUE
  if (is.null(.cov)) {
    .cov <- rxUiGet.monolixCovarianceEstimatesLin(x, ...)
    .sa <- FALSE
  }
  .j <- rxUiGet.monolixJacobian(x, ...)
  .n <- vapply(dimnames(.j)[[1]], function(n){
    paste0(.muRef[n], "_pop")
  }, character(1), USE.NAMES=FALSE)
  .cov <- .cov[.n, .n]
  if (.monolixCovarianceNeedsConversion(x, .sa)) {
    .jInv <- diag(1/diag(.j))
    # 2020+ SA returns correct matrix
    # 2021+ Lin returns correct matrix
    # Otherwise use inverse matrix to get the correct covariance
    .cov <- .jInv %*% .cov %*% .jInv
    dimnames(.cov) <- dimnames(.j)
  }
  .cov
}
