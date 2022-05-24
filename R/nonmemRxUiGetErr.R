.getErr <- function(err) {
  .err <- readLines(file.path(system.file(package="babelmixr2"), err))
  paste0("$ERROR\n",
         paste(gsub("^      ", "  ", .err), collapse="\n"))
}

#'@export
rxUiGet.nonmemErr <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  if (length(.predDf$cond) == 1L) {
    .iniDf <- .ui$iniDf
    if (.predDf$transform == "untransformed") return(.getErr("err.txt"))
    .transform <- paste(.predDf$transform)
    .low <- paste(.predDf$trLow)
    if (regexpr("^[0-9]+$", .low) != -1) .low <- paste0(.low, ".0")
    .hi <- paste(.predDf$trHi)
    if (regexpr("^[0-9]+$", .hi) != -1) .hi <- paste0(.hi, ".0")
    if (.transform == "boxCox") {
      .w <- which(.iniDf$err == "boxCox")
      .t <- (.nonmemGetThetaNum(.iniDf$name[.w],.ui))
      return(gsub("LAMBDA", .t, .getErr("err-boxCox.txt")))
    } else if (.transform == "yeoJohnson") {
      # The rest are not
      .w <- which(.iniDf$err == "yeoJohnson")
      .t <- (.nonmemGetThetaNum(.iniDf$name[.w],.ui))
      return(gsub("LAMBDA", .t, .getErr("err-yeoJohnson.txt")))
    } else if (.transform == "lnorm") {
      return(.getErr("err-lnorm.txt"))
    } else if (.transform == "logit") {
      return(gsub("LOW", .low,
                  gsub("HIGH", .hi,
                       .getErr("err-logit.txt"))))
    } else if (.transform == "logit + yeoJohnson") {
      .w <- which(.iniDf$err == "yeoJohnson")
      .t <- (.nonmemGetThetaNum(.iniDf$name[.w],.ui))
      return(gsub("LAMBDA", .t,
                  gsub("LOW", .low,
                       gsub("HIGH", .hi,
                            .getErr("err-logitYeoJohnson.txt")))))
    }
  }
  return(.getErr("err.txt"))
}
