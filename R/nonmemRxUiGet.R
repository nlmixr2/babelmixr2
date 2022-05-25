#' @export
rxUiGet.nonmemEtaTableName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(.ui$modelName, ".eta")
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
