
#' @export
rxUiGet.monolixModelFileName <- function(x, ...) {
  .ui <- x[[1]]
  file.path(getwd(), paste0(.ui$modelName, ".txt"))
}

#' @export
rxUiGet.monolixExportPath <- function(x, ...) {
  .ui <- x[[1]]
  file.path(getwd(), .ui$modelName)
}


#' @export
rxUiGet.mlxtranModel <- function(x, ...) {
  .ui <- x[[1]]
  # note there is some categorical covariates that are not taken care of here...
  paste0("<MODEL>\n\n",
         rxUiGet.mlxtranModelCovariate(x, ...),
         rxUiGet.mlxtranModelIndividual(x, ...),"\n\n",
         rxUiGet.mlxtranModelLongitudinal(x, ...))
}



