#' @export
rxUiGet.mlxtranModelIndividual <- function(x, ...) {
  .ui <- x[[1]]
  paste0("[INDIVIDUAL]\n",
         "input={}\n\n"
         "DEFINITION:\n",
         )
}

#' @export
rxUiGet.monolixModelFileName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(.ui$modelName, ".txt")
}



