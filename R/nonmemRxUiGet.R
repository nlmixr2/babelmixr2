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
  paste0(.ui$modelName, ".eta")
}

#' @export
rxUiGet.nonmemSdTableName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(.ui$modelName, ".pred")
}

#' @export
rxUiGet.nonmemContraName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(.ui$modelName, ".contra")
}

#' @export
rxUiGet.nonmemCcontraName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(.ui$modelName, ".ccontra")
}

#' @export
rxUiGet.nonmemCsv <- function(x, ...) {
  .ui <- x[[1]]
  paste0(.ui$modelName, ".csv")
}

#' @export
rxUiGet.nonmemNmctl <- function(x, ...) {
  .ui <- x[[1]]
  paste0(.ui$modelName, rxode2::rxGetControl(.ui, "extension", ".nmctl"))
}

#' @export
rxUiGet.nonmemNmlst <- function(x, ...) {
  .ui <- x[[1]]
  paste0(.ui$modelName, rxode2::rxGetControl(.ui, "outputExtension", ".lst"))
}

#' @export
rxUiGet.nonmemHashFile <- function(x, ...) {
  .ui <- x[[1]]
  paste0(.ui$modelName, ".md5")
}

#' @export
rxUiGet.nonmemQs <- function(x, ...) {
  .ui <- x[[1]]
  paste0(.ui$modelName, ".qs")
}
