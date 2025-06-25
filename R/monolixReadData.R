.monolixWaitForFile <- function(file, maxWait=10) {
  i <- 1
  while (i <= maxWait) {
    if (file.exists(file)) {
      return(TRUE)
    }
    Sys.sleep(1)
    i <- i + 1
  }
  stop("the file '", file, "' does not exist even though monolix export path is present",
       call.=FALSE)
}

#' @export
rxUiGet.monolixOutputVersion <- function(x, ...) {
  .exportPath <- rxUiGet.monolixExportPath(x, ...)
  if (!file.exists(.exportPath)) return(NULL)
  .summary <- file.path(.exportPath, "summary.txt")
  .monolixWaitForFile(.summary)
  .lines <- readLines(.summary, n=5)
  .w <- which(regexpr(".*[vV]ersion *: *[^ ]*.*", .lines) != -1)
  if (length(.w) == 1) {
    .line <- .lines[.w]
    # 2019 is 5.1.1
    return(sub(".*[vV]ersion *: *([^ ]*).*", "\\1", .line))
  }
  NULL
}

#' @export
rxUiGet.monolixHasChartData <- function(x, ...) {
  .chart <- rxUiGet.monolixCvParam(x, ...)
  if (file.exists(.chart)) return(TRUE)
  .mlxtran <- rxUiGet.monolixMlxtranFile(x, ...)
  if (!file.exists(.mlxtran)) return(FALSE)
  if (!.hasLixoftConnectors()) {
    return(FALSE)
  }
  .minfo("trying to create the chart data with lixoftConnectors::computeChartsData()")
  .x <- try(lixoftConnectors::loadProject(.mlxtran), silent=TRUE)
  if (inherits(.x, "try-error")) return(FALSE)
  .x <- try(lixoftConnectors::computeChartsData(), silent=TRUE)
  .minfo("done")
  file.exists(.chart)
}

#' @export
rxUiGet.monolixParHistoryRaw <- function(x, ...) {
  if (rxUiGet.monolixHasChartData(x, ...)) {
    return(read.csv(rxUiGet.monolixCvParam(x, ...)))
  }
  NULL
}

#' @export
rxUiGet.monolixParHistory <- function(x, ...) {
  .raw <- rxUiGet.monolixParHistoryRaw(x, ...)
  if (!is.null(.raw)) {
    .ui <- x[[1]]
    .iniDf <- .ui$iniDf
    .split <- .ui$getSplitMuModel
    .muRef <- c(.split$pureMuRef, .split$taintMuRef)
    .muRefCurEval <- .ui$muRefCurEval
    .eta <- .iniDf[is.na(.iniDf$ntheta), ]
    .theta <- .iniDf[!is.na(.iniDf$ntheta), ]
    .r <- .getOmegaR(.ui)
    .n <- names(.raw)
    .en <- dimnames(.r)[[1]]
    .covDataFrame <- .ui$saemMuRefCovariateDataFrame
    for (.i in seq_along(.eta$name)) {
      .n1 <- .eta$neta1[.i]
      .n2 <- .eta$neta2[.i]
      if (.n1 == .n2) {
        .par <- paste0("omega_", .mlxtranGetIndividualMuRefEtaMonolixName(.ui, .n1, .muRef))
        .n <- sub(.par, paste0("O(", .en[.n1], ")"), .n, fixed=TRUE)
      } else {
        .par <- paste0("corr_", .mlxtranGetIndividualMuRefEtaMonolixName(.ui, .n1, .muRef),"_",
                       .mlxtranGetIndividualMuRefEtaMonolixName(.ui, .n2, .muRef))
        .n <- sub(.par, paste0("C(", .en[.n1], ",", .en[.n2], ")"), .n, fixed=TRUE)
      }
    }
    for (.name in .theta$name) {
      .w <- which(.covDataFrame$covariateParameter == .name)
      if (length(.w) == 1) {
        .par <- paste0("beta_", .muRef[.covDataFrame$theta[.w]], "_",
                       .covDataFrame$covariate[.w])
        .n <- sub(.par, .name, .n, fixed=TRUE)
      }  else {
        .par <- paste0(.muRef[.name], "_pop")
        .w <- which(.n == .par)
        if (length(.w) == 1) {
          .n[.w] <- .par
        } else {
          .i <- which(.name == .theta$name)
          .isErr <- !is.na(.theta$ntheta[.i]) & !is.na(.theta$err[.i])
          if (.isErr) {
            .par <- eval(str2lang(paste0("rxToMonolix(", .name, ", ui=.ui)")))
            .n <- sub(.par, .name, .n, fixed=TRUE)
          }
        }
      }
    }
    .n <- sub("iteration", "iter", .n)
    .ret <- setNames(.raw, .n)
    .w <- which(.ret$phase == 2)[1]
    .niter <- .ret$iter[.w]
    .ret <- .ret[, names(.ret) != "phase"]
    .cls <- class(.ret)
    attr(.cls, "niter") <- .niter
    class(.ret) <- .cls
    return(.ret)
  }
  NULL
}

#' @export
rxUiGet.monolixPopulationParameters <- function(x, ...) {
  .exportPath <- rxUiGet.monolixExportPath(x, ...)
  if (!file.exists(.exportPath)) return(NULL)
  .popParameters <- file.path(.exportPath, "populationParameters.txt")
  .monolixWaitForFile(.popParameters)
  read.csv(.popParameters)
}

#' @export
rxUiGet.monolixOmega <- function(x, ...) {
  .ui <- x[[1]]
  .ret <- rxode2::rxGetControl(.ui, ".monolixOmega", NULL)
  if (!is.null(.ret)) return(.ret)
  .pop <- rxUiGet.monolixPopulationParameters(x, ...)
  if (is.null(.pop)) return(NULL)
  .iniDf <- .ui$iniDf
  .split <- .ui$getSplitMuModel
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
  if (!inherits(.sd, "matrix") && length(.sd) == 1) {
    .omega <- matrix(.sd * .sd, 1, 1, dimnames=list(names(.sd), names(.sd)))
  } else {
    diag(.r) <- 1
    .d <- diag(.sd)
    .omega <- .d %*% .r %*% .d
    .eta <- .eta[.eta$neta1 == .eta$neta2, ]
    dimnames(.omega) <- list(.eta$name, .eta$name)
  }
  rxode2::rxAssignControlValue(.ui, ".monolixOmega", .omega)
  .omega
}

.monolixGetPopParValue <- function(name, muRefCurEval, muRef, covDataFrame, pop) {
  .w <- which(covDataFrame$covariateParameter == name)
  if (length(.w) == 1) {
    .par <- paste0("beta_", muRef[covDataFrame$theta[.w]], "_",
                   covDataFrame$covariate[.w])
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

# This method allows rxUi objects to get the full theta values with
# `$monolixFullTheta`
#
# @param x argument list (as in all rxUiGet.* functions) which has
#   x[[1]] be the ui object
# @param ... other arguments (currently ignored)
# @return A named vector of the full theta values as needed for the UI object
# @author Matthew L. Fidler
# @keywords internal
#' @export
rxUiGet.monolixFullTheta <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  .hasVar <- any(names(.predDf) == "variance")
  .full <- rxode2::rxGetControl(.ui, ".monolixFullTheta", NULL)
  if (!is.null(.full)) return(.full)
  .pop <- rxUiGet.monolixPopulationParameters(x, ...)
  if (is.null(.pop)) return(NULL)
  .iniDf <- .ui$iniDf
  .theta <- .iniDf[!is.na(.iniDf$ntheta), ]
  .muRefCurEval <- .ui$muRefCurEval
  .split <- .ui$getSplitMuModel
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
                                    if (.hasVar) {
                                      ## Here it is the standard deviation,
                                      ##
                                      ## If it needs to be the variance
                                      ## it should be changed
                                      # First figure out the endpoint
                                      .cond <- .theta$condition[i]
                                      .w2 <- which(.predDf$cond == .cond)
                                      if (length(.w2) == 1) {
                                        # Now if the endpoint is a variance and
                                        # one of the parameters that is of variance type
                                        # return value^2
                                        if (.predDf$variance[.w2] &&
                                              .theta$err[i] %in% c("add", "prop", "propT",
                                                                   "pow", "powT", "logn",
                                                                   "dlogn", "lnorm",
                                                                   "dlnorm", "logitNorm",
                                                                   "probitNorm")) {
                                          return(.pop$value[.w]^2)
                                        }
                                      }
                                    }
                                    return(.pop$value[.w])
                                  }
                                  NA_real_
                                }, double(1), USE.NAMES=FALSE),
                         .theta$name)
  rxode2::rxAssignControlValue(.ui, ".monolixFullTheta", .fullTheta)
  .fullTheta
}

.bblIniDf <- function(theta, omega, ui) {
  .iniDf <- ui$iniDf
  .etas <- ui$iniDf[is.na(ui$iniDf$ntheta), ]
  .est <- c(theta,
            vapply(seq_along(.etas$neta1), function(i) {
              .n1 <- .etas$neta1[i]
              .n2 <- .etas$neta2[i]
              omega[.n1, .n2]
            }, double(1), USE.NAMES=FALSE))
  .iniDf$est <- .est
  .iniDf
}

#' @export
rxUiGet.monolixIniDf <- function(x, ...) {
  .omega <- rxUiGet.monolixOmega(x, ...)
  .theta <- rxUiGet.monolixFullTheta(x, ...)
  .ui <- x[[1]]
  .bblIniDf(.theta, .omega, .ui)
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
  if (!file.exists(.exportPath)) return(NULL)
  .individualParameters <- file.path(.exportPath, "IndividualParameters", "estimatedRandomEffects.txt")
  .monolixWaitForFile(.individualParameters)
  read.csv(.individualParameters)
}

#' @export
rxUiGet.monolixIndividualLL <- function(x, ...) {
  .ui <- x[[1]]
  .ret <- rxode2::rxGetControl(.ui, ".monolixIndividualLL", NULL)
  if (!is.null(.ret)) return(.ret)
  .exportPath <- rxUiGet.monolixExportPath(x, ...)
  if (!file.exists(.exportPath)) return(NULL)
  .individualParameters <- file.path(.exportPath, "LogLikelihood", "individualLL.txt")
  .monolixWaitForFile(.individualParameters)
  .ret <- read.csv(.individualParameters)
  rxode2::rxAssignControlValue(.ui, ".monolixIndividualLL", .ret)
  .ret
}

#' @export
rxUiGet.monolixEtaObf <- function(x, ...) {
  .ui <- x[[1]]
  .etas <- .ui$iniDf[!is.na(.ui$iniDf$neta1), ]
  .etas <- .etas[.etas$neta1 == .etas$neta2, ]
  .split <- .ui$getSplitMuModel
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
  .ui <- x[[1]]
  .ret <- rxode2::rxGetControl(.ui, ".monolixLL", NULL)
  if (!is.null(.ret)) return(.ret)
  .exportPath <- rxUiGet.monolixExportPath(x, ...)
  if (!file.exists(.exportPath)) return(NULL)
  .individualParameters <- file.path(.exportPath, "LogLikelihood", "logLikelihood.txt")
  .monolixWaitForFile(.individualParameters)
  .ret <- read.csv(.individualParameters)
  rxode2::rxAssignControlValue(.ui, ".monolixLL", .ret)
  .ret
}

#' @export
rxUiGet.monolixObjf <- function(x, ...) {
  .ll <- rxUiGet.monolixLL(x, ...)
  names(.ll)[2] <- "val"
  .w <- which(.ll$criteria == "OFV")
  if (length(.w) == 1L) {
    return(.ll[.w, "val"])
  } else {
    .w <- which(.ll$criteria == "-2LL")
    if (length(.w) == 1L) {
      return(.ll[.w, "val"])
    } else {
      stop("cannot figure out the objective function value from monolix",
           call.=FALSE)
    }
  }
}

#' @export
rxUiGet.monolixObjfType <- function(x, ...) {
  .ll <- rxUiGet.monolixObjf(x, ...)
  paste("monolix", names(.ll)[2])
}

 #' @export
rxUiGet.monolixCovarianceEstimatesLin <- function(x, ...) {
  .exportPath <- rxUiGet.monolixExportPath(x, ...)
  if (!file.exists(.exportPath)) return(NULL)
  .covLin <- file.path(.exportPath, "FisherInformation", "covarianceEstimatesLin.txt")
  if (!file.exists(.covLin)) return(NULL)
  .c <- read.csv(.covLin, header=FALSE)
  .n <- .c[, 1]
  .c <- as.matrix(.c[, -1])
  dimnames(.c) <- list(.n, .n)
  .c
}


#' @export
rxUiGet.monolixCovarianceEstimatesSA <- function(x, ...) {
  .exportPath <- rxUiGet.monolixExportPath(x, ...)
  if (!file.exists(.exportPath)) return(NULL)
  .covSA <- file.path(.exportPath, "FisherInformation", "covarianceEstimatesSA.txt")
  if (!file.exists(.covSA)) return(NULL)
  .c <- read.csv(.covSA, header=FALSE)
  .n <- .c[, 1]
  .c <- as.matrix(.c[, -1])
  dimnames(.c) <- list(.n, .n)
  .c
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
  .ui <- x[[1]]
  .split <- .ui$getSplitMuModel

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
  if (any(is.na(.cov)) || any(is.nan(.cov))) {
    warning("there are NaNs in the covariance in monolix, not included in nlmixr2 fit; check model results",
            call.=FALSE)
    return(NULL)
  }
  .cov
}
#' @export
rxUiGet.monolixPreds <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  .exportPath <- rxUiGet.monolixExportPath(x, ...)
  if (!file.exists(.exportPath)) return(NULL)
  .mlxtran <- monolix2rx::.monolixGetMlxtran(.ui)
  if (inherits(.mlxtran, "monolix2rxMlxtran")) {
    if (length(.predDf$var) > 1) {
      do.call("rbind", lapply(seq_along(.predDf$cond),
                              function(i) {
                                .var <- .predDf$cond[i]
                                .file <- file.path(.exportPath,
                                                   paste0("predictions_", .var, ".txt"))
                                .monolixWaitForFile(.file)
                                .ret <- read.csv(.file)
                                .ret$CMT <- .predDf$cond[i]
                                names(.ret) <- sub("id", "ID",
                                                   sub("time", "TIME",
                                                       sub(.var, "DV", names(.ret))))
                                .ret
                              }))
    } else {
      .var <- .predDf$cond
      .file <- file.path(.exportPath,"predictions.txt")
      .monolixWaitForFile(.file)
      .ret <- read.csv(.file)
      names(.ret) <- sub("id", "ID",
                         sub("time", "TIME",
                             sub(.var, "DV", names(.ret))))
      .ret

    }
  } else {
    if (length(.predDf$var) > 1) {
      do.call("rbind", lapply(seq_along(.predDf$var),
                              function(i) {
                                .var <- .predDf$var[i]
                                .file <- file.path(.exportPath,
                                                   paste0("predictions_rx_prd_", .var, ".txt"))
                                .monolixWaitForFile(.file)
                                .ret <- read.csv(.file)
                                .ret$CMT <- .predDf$cond[i]
                                names(.ret) <- sub("id", "ID",
                                                   sub("time", "TIME",
                                                       sub(paste0("rx_prd_", .var), "DV", names(.ret))))
                                .ret
                              }))
    } else {
      .var <- .predDf$var
      .file <- file.path(.exportPath,"predictions.txt")
      .monolixWaitForFile(.file)
      .ret <- read.csv(.file)
      names(.ret) <- sub("id", "ID",
                         sub("time", "TIME",
                             sub(paste0("rx_prd_", .var), "DV", names(.ret))))
      .ret
    }
  }
}

.monolixMergePredsAndCalcRelativeErr <- function(fit) {
  .singleEndpoint <- length(fit$ui$predDf$cond) == 1
  .tmp <- as.data.frame(fit)
  .tmp$ID <-as.integer(.tmp$ID)
  .by <- c("ID", "TIME", "CMT", "DV")
  if (.singleEndpoint) {
    .by <- c("ID", "TIME", "DV")
  } else {
    .tmp$CMT <- paste(.tmp$CMT)
  }
  .ret <- merge(fit$ui$monolixPreds, .tmp, by=.by)
  .ci <- (1 - fit$monolixControl$ci) / 2
  .q <- c(0, .ci, 0.5, 1 - .ci, 1)
  .qi <- stats::quantile(with(.ret, 100*abs((IPRED-indivPred_SAEM)/indivPred_SAEM)), .q, na.rm=TRUE)
  .qp <- stats::quantile(with(.ret, 100*abs((PRED-popPred)/popPred)), .q, na.rm=TRUE)
  .qai <- stats::quantile(with(.ret, abs(IPRED-indivPred_SAEM)), .q, na.rm=TRUE)
  .qap <- stats::quantile(with(.ret, abs((PRED-popPred)/popPred)), .q, na.rm=TRUE)
  .sigdig <- 3
  .msg <- c(paste0("IPRED relative difference compared to Monolix IPRED: ", round(.qi[3], 2),
                 "%; ", fit$monolixControl$ci * 100,"% percentile: (",
                 round(.qi[2], 2), "%,", round(.qi[4], 2), "%); rtol=", signif(.qi[3] / 100, digits=.sigdig)),
            paste0("PRED relative difference compared to Monolix PRED: ", round(.qp[3], 2),
                   "%; ", fit$monolixControl$ci * 100,"% percentile: (",
                   round(.qp[2], 2), "%,", round(.qp[4], 2), "%); rtol=", signif(.qp[3] / 100, digits=.sigdig)),
            paste0("IPRED absolute difference compared to Monolix IPRED: atol=",
                   signif(.qai[3], digits=.sigdig),
                 "; ", fit$monolixControl$ci * 100,"% percentile: (",
                 signif(.qai[2], digits=.sigdig), ", ", signif(.qai[4], digits=.sigdig), ")"),
            paste0("PRED absolute difference compared to Monolix PRED: atol=",
                   signif(.qap[3], digits=.sigdig),
                   "; ", fit$monolixControl$ci * 100,"% percentile: (",
                   signif(.qap[2], digits=.sigdig), ", ", signif(.qp[4], digits=.sigdig), ")"))
  list(individualRel=.qi , popRel=.qp,
       individualAbs=.qai, popAbs=.qap,
       message=.msg)
}
