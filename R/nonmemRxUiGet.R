#' @export
rxUiGet.nonmemEtaTableName <- function(x, ...) {
  .ui <- x[[1]]
  paste0("eta_", .ui$modelName)
}
