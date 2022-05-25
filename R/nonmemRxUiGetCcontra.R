#' @export
rxUiGet.nonmemContr <- function(x, ...) {
  #contr.txt
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  if (length(.predDf$cond) != 1) return(NULL)
  if (.predDf$transform == "untransformed") return(NULL)
  return(paste(readLines(file.path(system.file(package="babelmixr2"), "contr.txt")),
               collapse="\n"))
}

#'@export
rxUiGet.nonmemCcontra <- function(x, ...) {
  # ccontra.txt
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  .iniDf <- .ui$iniDf
  if (length(.predDf$cond) != 1) return(NULL)
  if (.predDf$transform == "untransformed") return(NULL)
  .transform <- paste(.predDf$transform)
  .low <- paste(.predDf$trLow)
  if (regexpr("^[0-9]+$", .low) != -1) .low <- paste0(.low, ".0")
  .hi <- paste(.predDf$trHi)
  if (regexpr("^[0-9]+$", .hi) != -1) .hi <- paste0(.hi, ".0")
  if (.transform == "boxCox") {
    # this one is included in PsN
    .w <- which(.iniDf$err == "boxCox")
    .t <- tolower(.nonmemGetThetaNum(.iniDf$name[.w],.ui))
    .ret <- paste(readLines(file.path(system.file(package="babelmixr2"), "ccontra-boxCox.txt")),
                  collapse="\n")
    return(gsub("LAMBDA", .t, .ret))
  } else if (.transform == "yeoJohnson") {
    # The rest are not
    .w <- which(.iniDf$err == "yeoJohnson")
    .t <- tolower(.nonmemGetThetaNum(.iniDf$name[.w],.ui))
    .ret <- paste(readLines(file.path(system.file(package="babelmixr2"), "ccontra-yeoJohnson.txt")),
                  collapse="\n")
    return(gsub("LAMBDA", .t, .ret))
  } else if (.transform == "lnorm") {
    return(paste(readLines(file.path(system.file(package="babelmixr2"), "ccontra-lnorm.txt")),
                  collapse="\n"))
  } else if (.transform == "logit") {
    return(gsub("LOW", .low,
                gsub("HIGH", .hi,
                     paste(readLines(file.path(system.file(package="babelmixr2"), "ccontra-logit.txt")),
                           collapse="\n"))))
  } else if (.transform == "logit + yeoJohnson") {
    .w <- which(.iniDf$err == "yeoJohnson")
    .t <- tolower(.nonmemGetThetaNum(.iniDf$name[.w],.ui))
    return(gsub("LAMBDA", .t,
                gsub("LOW", .low,
                     gsub("HIGH", .hi,
                          paste(readLines(file.path(system.file(package="babelmixr2"),
                                                    "ccontra-logitYeoJohnson.txt")),
                                collapse="\n")))))
  } else if (.transform == "logit + boxCox") {

  }
  NULL
}

