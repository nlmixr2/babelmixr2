#' @export
rxUiGet.nonmemExportPath <- function(x, ...) {
  .ui <- x[[1]]
  .extra <- ""
  .num <- rxode2::rxGetControl(.ui, ".modelNumber", 0)
  if (.num > 0) {
    .extra <- sprintf("-%03d", .num)
  }
  paste0(.ui$modelName, .extra, "-nonmem")
}

#' @export
rxUiGet.nonmemEtaTableName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemExportPath(x, ...), ".eta")
}

#' @export
rxUiGet.nonmemSdTableName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemExportPath(x, ...), ".pred")
}

#' @export
rxUiGet.nonmemContraName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemExportPath(x, ...), ".contra")
}

#' @export
rxUiGet.nonmemCcontraName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemExportPath(x, ...), ".ccontra")
}

#' @export
rxUiGet.nonmemCsv <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemExportPath(x, ...), ".csv")
}

#' @export
rxUiGet.nonmemNmctl <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemExportPath(x, ...), rxode2::rxGetControl(.ui, "extension", ".nmctl"))
}

#' @export
rxUiGet.nonmemQs <- function(x, ...) {
  .ui <- x[[1]]
  paste0(rxUiGet.nonmemExportPath(x, ...), ".qs")
}
