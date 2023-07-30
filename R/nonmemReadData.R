#' @export
rxUiGet.nonmemOutputLst <- function(x, ...) {
  .ui <- x[[1]]
  .info <- rxode2::rxGetControl(.ui, ".lstInfo", NULL)
  if (!is.null(.info)) return(.info)
  .f <-  .ui$file
  .exportPath <- rxUiGet.nonmemExportPath(x, ...)
  if (!is.null(.f)) {
    .info <- nonmem2rx::nminfo(.f)
  } else {
    .lst <- rxUiGet.nonmemLst(x, ...)
    if (!file.exists(file.path(.exportPath, .lst))) return(NULL)
    .info <- withr::with_dir(.exportPath, {
      nonmem2rx::nminfo(.lst)
    })
  }
  rxode2::rxAssignControlValue(.ui, ".lstInfo", .info)
  .info
}

#' @export
rxUiGet.nonmemOutputVersion <- function(x, ...) {
  .info <- rxUiGet.nonmemOutputLst(x, ...)
  if (is.null(.info)) return(NULL)
  .info$nonmem
}

#' @export
rxUiGet.nonmemOutputExt <- function(x, ...) {
  .ui <- x[[1]]
  .info <- rxode2::rxGetControl(.ui, ".ext", NULL)
  if (!is.null(.info)) return(.info)
  .exportPath <- rxUiGet.nonmemExportPath(x, ...)
  .ext <- rxUiGet.nonmemExt(x, ...)
  if (!file.exists(file.path(.exportPath, .ext)))  return(NULL)
  .ext <- withr::with_dir(.exportPath,
                          nonmem2rx::nmext(.ext))
  rxode2::rxAssignControlValue(.ui, ".ext", .ext)
  .ext
}

.getThetaNames <- function(ui) {
  if (exists("file", envir=ui)) {
    # here we are unsure of the theta name order since it could have
    # been rearranged
    .ext <- ui$nonmemOutputLst
    .theta <- .ext$theta
    .iniTheta <- ui$iniDf
    .iniTheta <- .iniTheta[!is.na(.iniTheta$ntheta), ]
    vapply(.theta, function(x) {
      .w <- which(abs(.iniTheta$est-x) < 1e-6)
      if (length(.w) == 1L) return(.iniTheta$name[.w])
      stop("there is a mismatch between estimates for nonmem2rx conversion",
           call.=FALSE)
    }, character(1), USE.NAMES=FALSE)
  } else {
    # here we control the theta name order
    .iniDf <- ui$iniDf
    .iniDf$name[!is.na(.iniDf$ntheta)]
  }
}

#' @export
rxUiGet.nonmemFullTheta <- function(x, ...) {
  .ui <- x[[1]]
  .ext <- rxUiGet.nonmemOutputLst(x, ...)
  if (is.null(.ext$theta)) {
    stop('cannot locate NONMEM output', call.=FALSE)
  }
  .ret <- setNames(.ext$theta, .getThetaNames(.ui))
  if (exists("file", envir=.ui)) {
    .iniDf <- .ui$iniDf
    # use same ordering of ui
    .ret <- .ret[.iniDf$name[!is.na(.iniDf$ntheta)]]
  }
  .ret
}

#' @export
rxUiGet.nonmemThetaDf <- function(x, ...) {
  .fullTheta <- rxUiGet.nonmemFullTheta(x, ...)
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .n <- names(.fullTheta)
  .fullTheta <- setNames(.fullTheta, NULL)
  .theta <- .iniDf[!is.na(.iniDf$ntheta), ]
  data.frame(lower=.theta$lower, theta=.fullTheta,
             fixed=.theta$fix, upper=.theta$upper,
             row.names=.n)
}

.getEtaNames <- function(ui) {
  if (exists("file", ui)) {
    # here we are unsure of the eta name order since it could have
    # been rearranged
    .ext <- ui$nonmemOutputLst
    .omega <- .ext$omega
    .iniEta <- ui$iniDf
    .iniEta <- .iniEta[which(.iniEta$neta1 == .iniEta$neta2), ]
    vapply(seq_along(diag(.omega)), function(x) {
      .w <- which(abs(.iniEta$est- .omega[x, x]) < 1e-6)
      if (length(.w) == 1L) return(.iniEta$name[.w])
      stop("there is a mismatch between eta estimates for nonmem2rx conversion",
           call.=FALSE)
    }, character(1), USE.NAMES=FALSE)
  } else {
    .iniDf <- ui$iniDf
    .iniDf <- .iniDf[is.na(.iniDf$ntheta), ]
    .iniDf <- .iniDf[.iniDf$neta1 == .iniDf$neta2, ]
    .iniDf$name
  }
}

#' @export
rxUiGet.nonmemOutputOmega <- function(x, ...) {
  .ui <- x[[1]]
  .n <- .getEtaNames(.ui)
  .ext <- rxUiGet.nonmemOutputLst(x, ...)
  .omega <- .ext$omega
  dimnames(.omega) <- list(.n, .n)
  if (exists("file", envir=.ui)) {
    # here we are unsure of the eta name order since it could have
    # been rearranged;
    # match ui
    .iniDf <- .ui$iniDf
    .iniDf <- .iniDf[is.na(.iniDf$ntheta), ]
    .iniDf <- .iniDf[.iniDf$neta1 == .iniDf$neta2, ]
    .omega <- .omega[.iniDf$name,.iniDf$name]
  }
  .omega
}

#' @export
rxUiGet.nonmemIniDf <- function(x, ...) {
  .omega <- rxUiGet.nonmemOutputOmega(x, ...)
  .theta <- rxUiGet.nonmemFullTheta(x, ...)
  .ui <- x[[1]]
  .bblIniDf(.theta, .omega, .ui)
}

#' @export
rxUiGet.nonmemEtaObf <- function(x, ...) {
  .ui <- x[[1]]
  if (exists("etaData", envir=.ui)) {
    .dat <- .ui$etaData
    .iniDf <- .ui$iniDf
    .iniDf <- .iniDf[is.na(.iniDf$ntheta), ]
    .iniDf <- .iniDf[.iniDf$neta1 == .iniDf$neta2, ]

    .dat <- .dat[,c("ID",.iniDf$name)]
    .dat$OBJI <- NA # not captured, could look...
    .dat$ID <- as.integer(.dat$ID)
    return(.dat)
  } else {
    .exportPath <- rxUiGet.nonmemExportPath(x, ...)
    .etaTable <- rxUiGet.nonmemEtaTableName(x, ...)
    if (!file.exists(file.path(.exportPath, .etaTable))) return(NULL)
    .ret <- withr::with_dir(.exportPath,
                            nonmem2rx::nmtab(.etaTable))
    .ret <- .ret[.ret$NMREP ==1, names(.ret) != "NMREP"]
    .n <- c("ID", .getEtaNames(.ui), "OBJI")
    names(.ret) <- .n
    .ret$ID <- as.integer(.ret$ID)
    .ret
  }
}

.getNonmemOrderNames <- function(ui) {
  .t <- .getThetaNames(ui)
  .s <- "_sigma"
  .e0 <- .getEtaNames(ui)
  .ef <- NULL
  .iniDf <- ui$iniDf
  for (.i in seq_along(.e0)) {
    for (.j in seq(1, .i)) {
      if (.i == .j) {
        .ef <- c(.ef, .e0[.i])
      } else {
        .v <- paste0("(", .e0[.j], ",", .e0[.i], ")")
        if (!(.v %in% .iniDf$name)) .v <- paste0("(", .e0[.i], ",", .e0[.j], ")")
        .ef <- c(.ef, .v)
      }
    }
  }
  c(.t, .s, .ef)
}

#' @export
rxUiGet.nonmemCovariance <- function(x, ...) {
  .ui <- x[[1]]
  if (exists("thetaMat", envir=.ui)) {
    .iniDf <- .ui$iniDf
    .n <- .iniDf$name[!is.na(.iniDf$ntheta)]
    .ui$thetaMat[.n, .n]
  } else {
    .exportPath <- rxUiGet.nonmemExportPath(x, ...)
    .covFile <- rxUiGet.nonmemCovFile(x, ...)
    if (!file.exists(file.path(.exportPath, .covFile))) return(NULL)
    .ret <- as.matrix(withr::with_dir(.exportPath,
                                      nonmem2rx::nmcov(.covFile)))
    .d <- .getNonmemOrderNames(.ui)
    dimnames(.ret) <- list(.d, .d)
    .t <- .getThetaNames(.ui)
    .ret[.t, .t]
  }
}

#' @export
rxUiGet.nonmemObjf <- function(x, ...) {
  .ret <- rxUiGet.nonmemOutputLst(x, ...)
  .ret$objf
}

#' @export
rxUiGet.nonmemParHistory <- function(x, ...) {
  .ui <- x[[1]]
  .exportPath <- rxUiGet.nonmemExportPath(x, ...)
  .ext <- rxUiGet.nonmemExt(x, ...)
  .ret <- withr::with_dir(.exportPath,
                          nonmem2rx::nmtab(.ext))
  .ret <- .ret[.ret$NMREP ==1, names(.ret) != "NMREP"]
  .d <- c("iter", .getNonmemOrderNames(.ui), "objf")
  names(.ret) <- .d
  .ret <- .ret[.ret$iter > 0, names(.ret) != "_sigma"]
  .n <- c("iter", .ui$iniDf$name, "objf")
  .ret <- .ret[, .n]
  .ret$type <- "Unscaled"
  .ret
}

#' @export
rxUiGet.nonmemObjfType <- function(x, ...) {
  .ui <- x[[1]]
  if (exists("nonmemData", .ui)) {
    return("nonmem2rx")
  }
  .est <- rxode2::rxGetControl(.ui, "est", "focei")
  if (.est %in% c("focei", "posthoc")) {
    return("nonmem focei")
  } else if (.est %in% "imp") {
    return("nonmem imp")
  } else if (.est %in% "its") {
    return("nonmem its")
  } else {
    stop("unknown objective type", call.=FALSE)
  }
}

#' @export
rxUiGet.nonmemRunTime <- function(x, ...) {
  .ui <- x[[1]]
  .info <- rxUiGet.nonmemOutputLst(x, ...)
  if (is.null(.info)) return(NULL)
  .info$time
}

#' @export
rxUiGet.nonmemPreds <- function(x, ...) {
  .ui <- x[[1]]
  if (exists("nonmemData", envir=.ui)) {
    .ipredData <- .ui$ipredData
    .w <- which(tolower(names(.ipredData)) == "id")
    if (length(.w) != 1L) return(NULL)
    names(.ipredData)[.w] <- "ID"
    .w <- which(tolower(names(.ipredData)) == "time")
    if (length(.w) != 1L) return(NULL)
    names(.ipredData)[.w] <- "TIME"
    .w <- which(tolower(names(.ipredData)) == "ipred")
    if (length(.w) != 1L) return(NULL)
    names(.ipredData)[.w] <- "IPRED"
    #[,c("ID", "TIME", "IPRED")]
    .ipredData <- .ipredData[,c("ID", "TIME", "IPRED")]
    .predData <- .ui$predData[,"PRED", drop=FALSE]
    .w <- which(tolower(names(.predData)) == "pred")
    if (length(.w) != 1L) return(NULL)
    names(.predData)[.w] <- "PRED"
    .predData <- .predData[,"PRED", drop=FALSE]
    if (length(.predData$PRED) == length(.ipredData$ID)) {
      .ret <- cbind(.ipredData, .predData)
      .ret$RXROW <- seq_along(.ipredData$ID)
      names(.ret) <- c("ID", "TIME", "nonmemIPRED", "nonmemPRED", "RXROW")
      return(.ret)
    }
    return(NULL)
  } else {
    .exportPath <- rxUiGet.nonmemExportPath(x, ...)
    .sdTable <- rxUiGet.nonmemSdTableName(x, ...)
    if (!file.exists(file.path(.exportPath, .sdTable))) return(NULL)
    .ret <- withr::with_dir(.exportPath,
                            nonmem2rx::nmtab(.sdTable))
    .ret <- .ret[.ret$NMREP ==1, names(.ret) != "NMREP"]
    setNames(.ret,
             c("ID", "TIME", "nonmemIPRED", "nonmemPRED", "RXROW"))
    
  }
}

#' @export
rxUiGet.nonmemTransMessage <- function(x, ...) {
  .lst <- rxUiGet.nonmemOutputLst(x, ...)
  if (is.null(.lst)) return(NULL)
  .lst$nmtran
}

#' @export
rxUiGet.nonmemTermMessage <- function(x, ...) {
  .lst <- rxUiGet.nonmemOutputLst(x, ...)
  if (is.null(.lst)) return(NULL)
  .lst$termInfo
}

#' @export
rxUiGet.nonmemSuccessful <- function(x, ...) {
  .term <- rxUiGet.nonmemTermMessage(x, ...)
  (regexpr("0MINIMIZATION SUCCESSFUL", .term) != -1)
}

#' @export
rxUiGet.nonmemRoundingErrors <- function(x, ...) {
  .term <- rxUiGet.nonmemTermMessage(x, ...)
  (regexpr("DUE TO ROUNDING ERRORS", .term) != -1)
}

.nonmemMergePredsAndCalcRelativeErr <- function(fit) {
  .np <- fit$ui$nonmemPreds
  if (is.null(.np)) {
    warning("without predictions output, absolute/relative difference between nlmixr2 and NONMEM predictions cannot be calculated",
            call.=FALSE)
    return(NULL)
  }
  .tmp <- as.data.frame(fit)
  .tmp$ID <-as.integer(.tmp$ID)
  .tmp$RXROW <- fit$env$.rownum
  .by <- c("ID", "TIME", "RXROW")
  .ret <- merge(.np, .tmp, by=.by)
-  if (!is.numeric(fit$nonmemControl$ci)) {
    .ci0 <- 0.95
  } else {
    .ci0 <- fit$nonmemControl$ci
  }
  .ci <- (1 - .ci0) / 2
  .q <- c(0, .ci, 0.5, 1 - .ci, 1)
  .qi <- stats::quantile(with(.ret, 100*abs((IPRED-nonmemIPRED)/nonmemIPRED)), .q, na.rm=TRUE)
  .qp <- stats::quantile(with(.ret, 100*abs((PRED-nonmemPRED)/nonmemPRED)), .q, na.rm=TRUE)
  .qai <- stats::quantile(with(.ret, abs(IPRED-nonmemIPRED)), .q, na.rm=TRUE)
  .qap <- stats::quantile(with(.ret, abs((PRED-nonmemPRED)/nonmemPRED)), .q, na.rm=TRUE)
  .sigdig <- 3
  .msg <- c(paste0("IPRED relative difference compared to Nonmem IPRED: ", round(.qi[3], 2),
                 "%; ", .ci0 * 100,"% percentile: (",
                 round(.qi[2], 2), "%,", round(.qi[4], 2), "%); rtol=", signif(.qi[3] / 100,
                                                                               digits=.sigdig)),
            paste0("PRED relative difference compared to Nonmem PRED: ", round(.qp[3], 2),
                   "%; ", .ci0 * 100,"% percentile: (",
                   round(.qp[2], 2), "%,", round(.qp[4], 2), "%); rtol=", signif(.qp[3] / 100,
                                                                                 digits=.sigdig)),
            paste0("IPRED absolute difference compared to Nonmem IPRED: ", .ci0 * 100,"% percentile: (",
                 signif(.qai[2], .sigdig), ", ", signif(.qai[4], .sigdig), "); atol=",
                   signif(.qai[3], .sigdig)),
            paste0("PRED absolute difference compared to Nonmem PRED: ", .ci0 * 100,"% percentile: (",
                   signif(.qap[2], .sigdig), ",", signif(.qp[4], .sigdig), "); atol=",
                   signif(.qap[3], .sigdig)))
  list(individualRel=.qi , popRel=.qp,
       individualAbs=.qai, popAbs=.qap,
       message=.msg)
}
