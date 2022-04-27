
#' @export
rxUiGet.getMonolixModel <- function(x, ...) {
  .ui <- x[[1]]
  .split <- rxUiGet.getSplitMuModel(x, ...)
  .predDf <- .ui$predDf
  # first drop the error lines
  .mainModel <- .split$modelWithDrop[-.predDf$line]
  .drop <- which(vapply(seq_along(.mainModel), function(x) {
    identical(.mainModel[[x]], quote(`_drop`))
  }, logical(1)))
  .mainModel <- .mainModel[-.drop]
}
