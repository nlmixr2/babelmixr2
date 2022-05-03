
#' @export
rxUiGet.monolixModelFileName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(.ui$modelName, ".txt")
}


#' @export
rxUiGet.mlxtranModel <- function(x, ...) {
  .ui <- x[[1]]
  paste0("<MODEL>\n\n",
         rxUiGet.mlxtranModelIndividual(x, ...),"\n\n",
         rxUiGet.mlxtranModelLongitudinal(x, ...))
}



