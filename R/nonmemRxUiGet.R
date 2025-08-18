#' @export
rxUiGet.nonmemModelName <- function(x, ...) {
  .ui <- x[[1]]
  if (exists("file", envir=.ui)) {
    return(sub("[.][^.]*$", "",basename(.ui$file)))
  }
  .modelName <- rxode2::rxGetControl(.ui, "modelName", NULL)
  if (is.null(.modelName)) {
    .modelName <- .ui$modelName
  } else {
    assign("modelName", .modelName, .ui)
  }
  if (isTRUE(checkmate::checkCharacter(.modelName, len=1, any.missing=FALSE))) {
    return(.modelName)
  }
  "x"
}
attr(rxUiGet.nonmemModelName, "rstudio") <- ""

#' @export
rxUiGet.nonmemExportPath <- function(x, ...) {
  .ui <- x[[1]]
  .f <- .ui$file
  if (!is.null(.f)) {
    return(dirname(.f))
  }
  .extra <- ""
  if (exists(".num", .ui)) {
    .num <- get(".num", .num, .ui)
  } else {
    .num <- rxode2::rxGetControl(.ui, ".modelNumber", 0)
  }
  if (.num > 0) {
    .extra <- sprintf("-%03d", .num)
    assign(".num", .num, .ui)
  }
  paste0(rxUiGet.nonmemModelName(x, ...), .extra, "-nonmem")
}
attr(rxUiGet.nonmemExportPath, "rstudio") <- ""

#' @export
rxUiGet.nonmemEtaTableName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".eta")
}
attr(rxUiGet.nonmemEtaTableName, "rstudio") <- ""

#' @export
rxUiGet.nonmemSdTableName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".pred")
}
attr(rxUiGet.nonmemSdTableName, "rstudio") <- ""

#' @export
rxUiGet.nonmemContraName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".contra")
}
attr(rxUiGet.nonmemContraName, "rstudio") <- ""

#' @export
rxUiGet.nonmemCcontraName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".ccontra")
}
attr(rxUiGet.nonmemCcontraName, "") <- ""

#' @export
rxUiGet.nonmemCsv <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".csv")
}
attr(rxUiGet.nonmemCsv, "rstudio") <- ""

#' @export
rxUiGet.nonmemNmctl <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), rxode2::rxGetControl(.ui, "extension", ".nmctl"))
}
attr(rxUiGet.nonmemNmctl, "rstudio") <- ""

#' @export
rxUiGet.nonmemNmlst <- function(x, ...) {
  .ui <- x[[1]]
  if (exists("outputExtension", envir=.ui)) {
    .lst <- .ui$outputExtension
  } else {
    .lst <- ".lst"
  }
  paste0(rxUiGet.nonmemModelName(x, ...), rxode2::rxGetControl(.ui, "outputExtension", .lst))
}
attr(rxUiGet.nonmemNmlst, "rstudio") <- ""

#' @export
rxUiGet.nonmemHashFile <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".md5")
}
attr(rxUiGet.nonmemHashFile, "rstudio") <- ""

#' @export
rxUiGet.nonmemQs <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...),
         ifelse(rxode2::rxGetControl(.ui, "readRounding", FALSE), "-rounding",
                ifelse(rxode2::rxGetControl(.ui, "readBadOpt", FALSE), "-bad-opt", "")),
         ".qs")
}
attr(rxUiGet.nonmemQs, "rstudio") <- ""

#' @export
rxUiGet.nonmemXml <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".xml")
}
attr(rxUiGet.nonmemXml, "rstudio") <- ""


#' @export
rxUiGet.nonmemLst <- function(x, ...) {
  .ui <- x[[1]]
  if (exists("outputExtension", envir=.ui)) {
    .lst <- .ui$outputExtension
  } else {
    .lst <- ".lst"
  }
  .lst <- rxode2::rxGetControl(.ui, "outputExtension", .lst)
  paste0(rxUiGet.nonmemModelName(x, ...), .lst)
}
attr(rxUiGet.nonmemLst, "rstudio") <- ""

#' @export
rxUiGet.nonmemExt <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".ext")
}
attr(rxUiGet.nonmemExt, "rstudio") <- ""


#' @export
rxUiGet.nonmemCovFile <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".cov")
}
attr(rxUiGet.nonmemCovFile, "rstudio") <- ""
