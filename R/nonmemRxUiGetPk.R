.nonmemGetMuNum <- function(theta, ui) {
  .muRefDf <- ui$muRefDataFrame
  .iniDf <- ui$iniDf
  .w <- which(.muRefDf$theta == theta)
  if (length(.w) != 1) return(NA_character_)
  .eta <- .muRefDf$eta[.w]
  .w <- which(.iniDf$name == .eta)
  if (length(.w) != 1) return(NA_character_)
  paste0("MU_", .iniDf$neta1[.w])
}

.nonmemGetThetaNum <- function(theta, ui) {
  .iniDf <- ui$iniDf
  .w <- which(.iniDf$name == theta)
  if (length(.w) != 1) return(NA_character_)
  if (is.na(.iniDf$ntheta[.w])) return(NA_character_)
  paste0("THETA(", .iniDf$ntheta[.w], ")")
}

.nonmemGetEtaNum <- function(theta, ui) {
  .iniDf <- ui$iniDf
  .w <- which(.iniDf$name == theta)
  if (length(.w) != 1) return(NA_character_)
  if (is.na(.iniDf$neta1[.w])) return(NA_character_)
  paste0("ETA(", .iniDf$neta1[.w], ")")
}

.nonmemGetThetaMuCov <- function(theta, ui, covRefDf) {
  .w <- which(covRefDf$theta == theta)
  if (length(.w) == 0) return(NA_character_)
  paste(paste0(covRefDf$covariate[.w], "*",
               vapply(covRefDf$covariateParameter[.w], .nonmemGetThetaNum, character(1),
                      ui=ui)),
        collapse="+")
}

#'@export
rxUiGet.nonmemThetaRep <- function(x, ...) {
  .split <- rxUiGet.getSplitMuModel(x, ...)
  .muRef <- c(.split$pureMuRef, .split$taintMuRef)
  .thetas <- names(.muRef)
  .ui <- x[[1]]
  .covRefDf <- .ui$saemMuRefCovariateDataFrame
  .ret <- data.frame(theta=.thetas,
                     nmTheta=vapply(.thetas, .nonmemGetThetaNum, character(1), ui=.ui,
                                    USE.NAMES=FALSE),
             mu=vapply(.thetas, .nonmemGetMuNum, character(1), ui=.ui,
                       USE.NAMES=FALSE),
             cov=vapply(.thetas, .nonmemGetThetaMuCov, character(1),
                        ui=.ui, covRefDf=.covRefDf, USE.NAMES=FALSE))
  .ret$nmEta <- paste0("ETA(",substr(.ret$mu,4, 10),")")
  .ret
}

#'@export
rxUiGet.nonmemPkDes <- function(x, ...) {
  .ui <- x[[1]]
  .split <- rxUiGet.getSplitMuModel(x, ...)
  .mu <- rxUiGet.nonmemThetaRep(x, ...)
  .ret <- vapply(seq_along(.mu$mu), function(i) {
    paste0("  ", .mu$mu[i], "=", .mu$nmTheta[i],
           ifelse(is.na(.mu$cov[i]), "",
                  paste0("+", .mu$cov[i])))
  }, character(1), USE.NAMES=FALSE)
  .ret <- paste(.ret, collapse="\n")
  .mu2 <- setNames(paste0(.mu$mu, "+", .mu$nmEta), .mu$theta)
  assign(".thetaMu", .mu2, envir=.ui)
  on.exit({
    if (exists(".thetaMu", envir=.ui)) {
      rm(".thetaMu", envir=.ui)
    }
  })
  .pk <- paste0("$PK\n",
                 .ret,"\n",
                 paste(vapply(seq_along(.split$muRefDef),
                              function(i) {
                                .rxToNonmem(.split$muRefDef[[i]], ui=.ui)
                              }, character(1), USE.NAMES=FALSE),
                       collapse="\n"))

  .mainModel <- rxode2::rxCombineErrorLines(.ui,
                                            errLines=nmGetDistributionNonmemLines(.ui),
                                            paramsLine=NA,
                                            modelVars=TRUE,
                                            cmtLines=FALSE,
                                            dvidLine=FALSE,
                                            lstExpr=.split$modelWithDrop,
                                            useIf=TRUE)
  .norm <- rxode2::rxNorm(eval(.mainModel))
  .des <- rxToNonmem(.norm, ui=.ui)
  .prop <- .nonmemGetCmtProperties(.ui)
  .pk2 <- vapply(seq_along(.prop$cmt),
                 function(i) {
                   .cmt <- .prop$cmt[i]
                   .ret <- NULL
                   if (!is.na(.prop$f[i])) {
                     .ret <- c(.ret,
                               paste0("  F", .cmt, "=", .prop$f[i]))
                   }
                   if (!is.na(.prop$dur[i])) {
                     .ret <- c(.ret,
                               paste0("  DUR", .cmt, "=", .prop$dur[i]))
                   }
                   if (!is.na(.prop$lag[i])) {
                     .ret <- c(.ret,
                               paste0("  ALAG", .cmt, "=", .prop$lag[i]))
                   }
                   if (!is.na(.prop$init[i])) {
                     .ret <- c(.ret,
                               paste0("  A_0(", .cmt, ")=", .prop$init[i]))
                   }
                   if (is.null(.ret)) return(NA_character_)
                   paste(.ret, collapse="\n")
                 }, character(1), USE.NAMES=FALSE)
  .pk2 <- .pk2[!is.na(.pk2)]
  .pk2 <- ifelse(length(.pk2) > 0, paste0("\n", .pk2), "")
  rm(".thetaMu", envir=.ui)
  paste0(.pk, .pk2,
         "\n\n$DES\n",
         .des)
}
