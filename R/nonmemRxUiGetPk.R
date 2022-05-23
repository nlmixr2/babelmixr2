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
                       collapse="\n")
                 )
  rm(".thetaMu", envir=.ui)
  .pk
}


