.popedInfo <- new.env(parent=emptyenv())

#' @export
rxUiGet.popedFfFun <- function(x, ...) {
  .ui <- x[[1]]
}

.popedGetBpopNum <- function(theta, ui) {
  .iniDf <- ui$iniDf
  .w <- which(.iniDf$name == theta)
  if (length(.w) != 1) return(NA_character_)
  if (is.na(.iniDf$ntheta[.w])) return(NA_character_)
  paste0("bpop[", .iniDf$ntheta[.w], "]")
}


#'@export
rxUiGet.popedBpopRep <- function(x, ...) {
  .ui <- x[[1]]
  .split <- .ui$getSplitMuModel
  .muRef <- c(.split$pureMuRef, .split$taintMuRef)
  .thetas <- names(.muRef)
  .covRefDf <- .ui$saemMuRefCovariateDataFrame
  .ret <- data.frame(theta=.thetas,
                     bpop=vapply(.thetas, .popedGetBpopNum,
                                    character(1), ui=.ui, USE.NAMES=FALSE),
                     mu=vapply(.thetas, .nonmemGetMuNum, character(1), ui=.ui,
                               USE.NAMES=FALSE),
                     cov=vapply(.thetas, .nonmemGetThetaMuCov, character(1),
                                ui=.ui, covRefDf=.covRefDf, USE.NAMES=FALSE))
  .ret$b <- ifelse(is.na(.ret$mu), NA_character_,
                       paste0("b[",substr(.ret$mu,4, 10),"]"))
  .ret
}


rxUiGet.popedFgFun <- function(x, ...) {
  # function(x, a, bpop, b, bocc)
  # x=?
  # a=covariates (could be dose, tau etc)
  # bpop = population variables
  # b = eta variables
  # bocc = occasion variables
  .ui <- x[[1]]
  .split <- .ui$getSplitMuModel
  .mu <- rxUiGet.popedBpopRep(x, ...)
  .ret <- vapply(seq_along(.mu$mu), function(i) {
    if (is.na(.mu$mu[i])) return(NA_character_)
    paste0(.mu$mu[i], " <- ", .mu$bpop[i],
           ifelse(is.na(.mu$cov[i]), "",
                  paste0("+", .mu$cov[i])))
  }, character(1), USE.NAMES=FALSE)
  .ret <- .ret[!is.na(.ret)]
  .mu2 <- ifelse(is.na(.mu$mu),
                          paste0(.mu$theta, "<-", .mu$bpop),
                 paste0(.mu$theta, "<-", .mu$mu, "+", .mu$b))
  .covDef <- .ui$allCovs
  .covDefLst <- lapply(seq_along(.covDef),
                       function(i) {
                         str2lang(paste0(.covDef[i], "<- a[", i, "]"))
                       })
  .v <- c(.split$pureMuRef, .split$taintMuRef, .covDef)
  .body1 <- c(list(quote(`{`)),
              lapply(c(.ret, .mu2),
                     function(x) {
                       str2lang(x)
                     }),
              .covDefLst,
              .split$muRefDef,
              list(str2lang(paste("c(", paste(paste0(.v, "=", .v), collapse=","), ")"))))
  #
  .body1 <- as.call(.body1)
  .f <- function(x, a, bpop, b, bocc) {}
  body(.f) <- .body1
  .f
}
attr(rxUiGet.popedFfFun, "desc") <- "PopED parameter model (fg_fun)"

rxUiGet.popedBpop <- function(x, ...) {

}

rxUiGet.popedNotfixedBpop <- function(x, ...) {

}
