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
  .ret <- rxode2::rxGetControl(.ui, ".monolixOmega", NULL)
  if (!is.null(.ret)) return(.ret)
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
  dimnames(.omega) <- list(.eta$name, .eta$name)
  rxode2::rxAssignControlValue(.ui, ".monolixOmega", .omega)
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
rxUiGet.monolixFullTheta <- function(x, ...) {
  .ui <- x[[1]]
  .full <- rxode2::rxGetControl(.ui, ".monolixFullTheta", NULL)
  if (!is.null(.full)) return(.full)
  .pop <- rxUiGet.monolixPopulationParameters(x, ...)
  if (is.null(.pop)) return(NULL)
  .iniDf <- .ui$iniDf
  .theta <- .iniDf[!is.na(.iniDf$ntheta), ]
  .muRefCurEval <- .ui$muRefCurEval
  .split <- rxUiGet.getSplitMuModel(x, ...)
  .muRef <- c(.split$pureMuRef, .split$taintMuRef)
  .covDataFrame <- .ui$saemMuRefCovariateDataFrame
  .fullTheta <- setNames(vapply(seq_along(.theta$name),
                                function(i) {
                                  .n <- .theta$name[i]
                                  .isPop <- is.na(.theta$err[i])
                                  if (.isPop) {
                                    return(.monolixGetPopParValue(.n, .muRefCurEval, .muRef, .covDataFrame, .pop))
                                  }
                                  .isErr <- !is.na(.theta$ntheta[i]) & !is.na(.theta$err[i])
                                  if (.isErr) {
                                    .par <- eval(str2lang(paste0("rxToMonolix(", .n, ", ui=.ui)")))
                                    .w <- which(.pop$parameter == .par)
                                    if (length(.w) != 1) return(NA_real_)
                                    return(.pop$value[.w])
                                  }
                                  NA_real_
                                }, double(1), USE.NAMES=FALSE),
                         .theta$name)
  rxode2::rxAssignControlValue(.ui, ".monolixFullTheta", .fullTheta)
  .fullTheta
}

#' @export
rxUiGet.monolixIniDf <- function(x, ...) {
  .omega <- rxUiGet.monolixOmega(x, ...)
  .theta <- rxUiGet.monolixFullTheta(x, ...)
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .etas <- .ui$iniDf[is.na(.ui$iniDf$ntheta), ]
  .est <- c(.theta,
            vapply(seq_along(.etas$neta1), function(i) {
              .n1 <- .etas$neta1[i]
              .n2 <- .etas$neta2[i]
              .omega[.n1, .n2]
            }, double(1), USE.NAMES=FALSE))
  .iniDf$est <- .est
  .iniDf
}

#' @export
rxUiGet.monolixTheta <- function(x, ...) {
  .fullTheta <- rxUiGet.monolixFullTheta(x, ...)
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .n <- names(.fullTheta)
  .fullTheta <- setNames(.fullTheta, NULL)
  .theta <- .iniDf[!is.na(.iniDf$ntheta), ]
  data.frame(lower=.theta$lower, theta=.fullTheta,
             fixed=.theta$fix, upper=.theta$upper,
             row.names=.n)
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
rxUiGet.monolixIndividualLL <- function(x, ...) {
  .exportPath <- rxUiGet.monolixExportPath(x, ...)
  .individualParameters <- file.path(.exportPath, "LogLikelihood", "individualLL.txt")
  if (file.exists(.individualParameters)) {
    return(read.csv(.individualParameters))
  }
  NULL
}

#' @export
rxUiGet.monolixEtaObf <- function(x, ...) {
  .ui <- x[[1]]
  .etas <- .ui$iniDf[!is.na(.ui$iniDf$neta1), ]
  .etas <- .etas[.etas$neta1 == .etas$neta2, ]
  .split <- rxUiGet.getSplitMuModel(x, ...)
  .muRef <- c(.split$pureMuRef, .split$taintMuRef)
  .etaMonolix <- rxUiGet.monolixIndividualParameters(x, ...)
  if (is.null(.etaMonolix)) return(NULL)
  .n <- c("id", vapply(.etas$neta1, function(i) {
    paste0("eta_",   .mlxtranGetIndividualMuRefEtaMonolixName(.ui, i, .muRef), "_SAEM")
  }, character(1), USE.NAMES=FALSE))
  .etaObf <- .etaMonolix[, .n]
  names(.etaObf) <- c("ID", .etas$name)
  .indLL <- rxUiGet.monolixIndividualLL(x, ...)
  if (is.null(.indLL)) {
    .etaObf$OBJI <- NA_real_
  } else {
    names(.indLL) <- c("ID", "OBJI")
    .etaObf <- merge(.etaObf, .indLL, by="ID")
    .etaObf <- .etaObf[, c("ID", .etas$name, "OBJI")]
  }
  .etaObf
}

#' @export
rxUiGet.monolixLL <- function(x, ...) {
  .exportPath <- rxUiGet.monolixExportPath(x, ...)
  .individualParameters <- file.path(.exportPath, "LogLikelihood", "logLikelihood.txt")
  if (file.exists(.individualParameters)) {
    return(read.csv(.individualParameters))
  }
  NULL
}


 #' @export
rxUiGet.monolixCovarianceEstimatesLin <- function(x, ...) {
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
  .ui <- x[[1]]
  rxode2::rxAssignControlValue(.ui, ".covMethod", ifelse(.sa, "MonolixSA", "MonolixLin"))
  if (.monolixCovarianceNeedsConversion(x, .sa)) {
    .jInv <- diag(1/diag(.j))
    # 2020+ SA returns correct matrix
    # 2021+ Lin returns correct matrix
    # Otherwise use inverse matrix to get the correct covariance
    .cov <- .jInv %*% .cov %*% .jInv
    dimnames(.cov) <- dimnames(.j)
    rxode2::rxAssignControlValue(.ui, ".covMethod", ifelse(.sa, "MonolixSA*", "MonolixLin*"))
  }

  .cov
}
