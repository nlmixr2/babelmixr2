#' @export
rxUiGet.nonmemOutputXml <- function(x, ...) {
  .ui <- x[[1]]
  .info <- rxode2::rxGetControl(.ui, ".xml", NULL)
  if (!is.null(.info)) return(.info)
  .exportPath <- rxUiGet.nonmemExportPath(x, ...)
  .xml <- rxUiGet.nonmemXml(x, ...)
  if (file.exists(file.path(.exportPath, .xml))) {
    .info <- withr::with_dir(.exportPath, {
      pmxTools::read_nm(.xml, quiet=TRUE)
    })
    rxode2::rxAssignControlValue(.ui, ".xml", .info)
    return(.info)
  }
  NULL
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

#' @export
rxUiGet.nonmemFullTheta <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .ext <- rxUiGet.nonmemOutputExt(x, ...)
  setNames(.ext$Thetas, .iniDf$name[!is.na(!.iniDf$ntheta)])
}

#' @export
rxUiGet.nonmemOutputOmega <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .iniDf <- .iniDf[is.na(.iniDf$ntheta), ]
  .iniDf <- .iniDf[.iniDf$neta1 == .iniDf$neta2, ]
  .n <- .iniDf$name
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

