#' @export
rxUiGet.nonmemModelName <- function(x, ...) {
  .ui <- x[[1]]
  .modelName <- .ui$modelName
  if (isTRUE(checkmate::checkCharacter(.modelName, len=1, any.missing=FALSE))) {
    return(.modelName)
  }
  "x"
}

#' @export
rxUiGet.nonmemExportPath <- function(x, ...) {
  .ui <- x[[1]]
  .extra <- ""
  .num <- rxode2::rxGetControl(.ui, ".modelNumber", 0)
  if (.num > 0) {
    .extra <- sprintf("-%03d", .num)
  }
  paste0(rxUiGet.nonmemModelName(x, ...), .extra, "-nonmem")
}

#' @export
rxUiGet.nonmemEtaTableName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".eta")
}

#' @export
rxUiGet.nonmemSdTableName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".pred")
}

#' @export
rxUiGet.nonmemContraName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".contra")
}

#' @export
rxUiGet.nonmemCcontraName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".ccontra")
}

#' @export
rxUiGet.nonmemCsv <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".csv")
}

#' @export
rxUiGet.nonmemNmctl <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), rxode2::rxGetControl(.ui, "extension", ".nmctl"))
}

#' @export
rxUiGet.nonmemNmlst <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), rxode2::rxGetControl(.ui, "outputExtension", ".lst"))
}

#' @export
rxUiGet.nonmemHashFile <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".md5")
}

#' @export
rxUiGet.nonmemQs <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...),
         ifelse(rxode2::rxGetControl(.ui, "readRounding", FALSE), "-rounding",
                ifelse(rxode2::rxGetControl(.ui, "readBadOpt", FALSE), "-bad-opt", "")),
         ".qs")
}


#' @export
rxUiGet.nonmemXml <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".xml")
}

#' @export
rxUiGet.nonmemExt <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".ext")
}


#' @export
rxUiGet.nonmemCovFile <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemModelName(x, ...), ".cov")
}
