#' @export
rxUiGet.nonmemOutputXml <- function(x, ...) {
  .ui <- x[[1]]
  .info <- rxode2::rxGetControl(.ui, ".xml", NULL)
  if (!is.null(.info)) return(.info)
  .exportPath <- rxUiGet.nonmemExportPath(x, ...)
  .xml <- rxUiGet.nonmemXml(x, ...)
  if (!file.exists(file.path(.exportPath, .xml))) return(NULL)
  .info <- withr::with_dir(.exportPath, {
    pmxTools::read_nm(.xml, quiet=TRUE)
  })
  rxode2::rxAssignControlValue(.ui, ".xml", .info)
  .info
}

#' @export
rxUiGet.nonmemOutputVersion <- function(x, ...) {
  .info <- rxUiGet.nonmemOutputXml(x, ...)
  if (is.null(.info)) return(NULL)
  .ver <- .info$nonmem$program_information[[1]]
  gsub(".*VERSION +([^ ]*).*", "\\1", gsub("\n", " ", .ver))
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
                          pmxTools::read_nmext(.ext, quiet=TRUE))
  rxode2::rxAssignControlValue(.ui, ".ext", .ext)
  .ext
}

.getThetaNames <- function(ui) {
  .iniDf <- ui$iniDf
  .iniDf$name[!is.na(!.iniDf$ntheta)]
}

#' @export
rxUiGet.nonmemFullTheta <- function(x, ...) {
  .ui <- x[[1]]
  .ext <- rxUiGet.nonmemOutputExt(x, ...)
  setNames(.ext$Thetas, .getThetaNames(.ui))
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
  .iniDf <- ui$iniDf
  .iniDf <- .iniDf[is.na(.iniDf$ntheta), ]
  .iniDf <- .iniDf[.iniDf$neta1 == .iniDf$neta2, ]
  .iniDf$name
}

#' @export
rxUiGet.nonmemOutputOmega <- function(x, ...) {
  .ui <- x[[1]]
  .n <- .getEtaNames(.ui)
  .ext <- rxUiGet.nonmemOutputExt(x, ...)
  .omegaLst <- .ext$Omega
  .len <- length(.omegaLst)
  .v <- NULL
  for (.i in seq_along(.omegaLst)) {
    .v <- c(.v, .omegaLst[[.i]])
  }
  eval(str2lang(paste("lotri::lotri(",
                      paste(paste(.n, collapse="+"),
                            "~", deparse1(.v)),
                      ")")))
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
  .exportPath <- rxUiGet.nonmemExportPath(x, ...)
  .etaTable <- rxUiGet.nonmemEtaTableName(x, ...)
  if (!file.exists(file.path(.exportPath, .etaTable))) return(NULL)
  .ret <- withr::with_dir(.exportPath,
                          pmxTools::read_nm_multi_table(.etaTable))
  .n <- c("ID", .getEtaNames(.ui), "OBJI")
  names(.ret) <- .n
  .ret$ID <- as.integer(.ret$ID)
  .ret
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
  .exportPath <- rxUiGet.nonmemExportPath(x, ...)
  .covFile <- rxUiGet.nonmemCovFile(x, ...)
  if (!file.exists(file.path(.exportPath, .covFile))) return(NULL)
  .ret <- as.matrix(withr::with_dir(.exportPath,
                                    pmxTools::read_nmcov(.ui$modelName, quiet=TRUE)))
  .d <- .getNonmemOrderNames(.ui)
  dimnames(.ret) <- list(.d, .d)
  .t <- .getThetaNames(.ui)
  .ret[.t, .t]
}

#' @export
rxUiGet.nonmemObjf <- function(x, ...) {
  .ret <- rxUiGet.nonmemOutputExt(x, ...)
  .ret$OFV
}

#' @export
rxUiGet.nonmemParHistory <- function(x, ...) {
  .ui <- x[[1]]
  .exportPath <- rxUiGet.nonmemExportPath(x, ...)
  .ext <- rxUiGet.nonmemExt(x, ...)
  .ret <- withr::with_dir(.exportPath,
                          pmxTools::read_nm_multi_table(.ext))
  .d <- c("iter", .getNonmemOrderNames(.ui), "objf")
  names(.ret) <- .d
  .ret <- .ret[.ret$iter > 0, names(.ret) != "_sigma"]
  .n <- c("iter", .ui$iniDf$name, "objf")
  .ret[, .n]
}

#' @export
rxUiGet.nonmemObjfType <- function(x, ...) {
  .ui <- x[[1]]
  .est <- rxode2::rxGetControl(.ui, "est", "focei")
  if (.est %in% c("focei", "posthoc")) {
    return("nonmem focei")
  } else {
    stop("unknown objective type", call.=FALSE)
  }
}

#' @export
rxUiGet.nonmemRunTime <- function(x, ...) {
  .ui <- x[[1]]
  .xml <- rxUiGet.nonmemOutputXml(x, ...)
  if (is.null(.xml)) return(NULL)
  .start <- as.POSIXct(.xml$start_datetime[[1]],format="%Y-%m-%dT%H:%M:%S",tz=Sys.timezone())
  .stop <- as.POSIXct(.xml$stop_datetime[[1]],format="%Y-%m-%dT%H:%M:%S",tz=Sys.timezone())
  as.numeric(difftime(.stop, .start, tz=Sys.timezone(), units = "secs"))
}

#' @export
rxUiGet.nonmemPreds <- function(x, ...) {
  .exportPath <- rxUiGet.nonmemExportPath(x, ...)
  .sdTable <- rxUiGet.nonmemSdTableName(x, ...)
  if (!file.exists(file.path(.exportPath, .sdTable))) return(NULL)
  setNames(withr::with_dir(.exportPath,
                           pmxTools::read_nm_multi_table(.sdTable)),
           c("ID", "TIME", "nonmemIPRED", "nonmemPRED", "RXROW"))
}

#' @export
rxUiGet.nonmemTransMessage <- function(x, ...) {
  .xml <- rxUiGet.nonmemOutputXml(x, ...)
  if (is.null(.xml)) return(NULL)
  .xml$nmtran[[1]]
}

#' @export
rxUiGet.nonmemTermMessage <- function(x, ...) {
  .xml <- rxUiGet.nonmemOutputXml(x, ...)
  if (is.null(.xml)) return(NULL)
  .xml$nonmem$problem$estimation$termination_information[[1]]
}

#' @export
rxUiGet.nonmemSuccessful <- function(x, ...) {
  .xml <- rxUiGet.nonmemOutputXml(x, ...)
  if (is.null(.xml)) return(NULL)
  .term <- as.integer(.xml$nonmem$problem$termination_status[[1]])
  (.term == 0)
}

.nonmemMergePredsAndCalcRelativeErr <- function(fit) {
  .tmp <- as.data.frame(fit)
  .tmp$ID <-as.integer(.tmp$ID)
  .tmp$RXROW <- fit$env$.rownum
  .by <- c("ID", "TIME", "RXROW")
  .ret <- merge(fit$ui$nonmemPreds, .tmp, by=.by)
  .ci <- (1 - fit$nonmemControl$ci) / 2
  .q <- c(0, .ci, 0.5, 1 - .ci, 1)
  .qi <- quantile(with(.ret, 100*abs((IPRED-nonmemIPRED)/nonmemIPRED)), .q, na.rm=TRUE)
  .qp <- quantile(with(.ret, 100*abs((PRED-nonmemPRED)/nonmemPRED)), .q, na.rm=TRUE)
  .qai <- quantile(with(.ret, abs(IPRED-nonmemIPRED)), .q, na.rm=TRUE)
  .qap <- quantile(with(.ret, abs((PRED-nonmemPRED)/nonmemPRED)), .q, na.rm=TRUE)
  .sigdig <- 3
  .msg <- c(paste0("IPRED relative difference compared to Nonmem IPRED: ", round(.qi[3], 2),
                 "%; ", fit$nonmemControl$ci * 100,"% percentile: (",
                 round(.qi[2], 2), "%,", round(.qi[4], 2), "%); rtol=", signif(.qi[3] / 100, .sigdig)),
            paste0("PRED relative difference compared to Nonmem PRED: ", round(.qp[3], 2),
                   "%; ", fit$nonmemControl$ci * 100,"% percentile: (",
                   round(.qp[2], 2), "%,", round(.qp[4], 2), "%); rtol=", signif(.qp[3] / 100, , .sigdig)),
            paste0("IPRED absolute difference compared to Nonmem IPRED: atol=",
                   signif(.qai[3], .sigdig),
                 "; ", fit$nonmemControl$ci * 100,"% percentile: (",
                 signif(.qai[2], .sigdig), ", ", signif(.qai[4], .sigdig), ")"),
            paste0("PRED absolute difference compared to Nonmem PRED: atol=",
                   signif(.qap[3], .sigdig),
                   "; ", fit$nonmemControl$ci * 100,"% percentile: (",
                   signif(.qap[2], .sigdig), ",
", signif(.qp[4], .sigdig), ")"))
  list(individualRel=.qi , popRel=.qp,
       individualAbs=.qai, popAbs=.qap,
       message=.msg)
}
