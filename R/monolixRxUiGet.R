
#' @export
rxUiGet.monolixExportPath <- function(x, ...) {
  .ui <- x[[1]]
  file.path(getwd(), .ui$modelName)
}

#' @export
rxUiGet.monolixModelFileName <- function(x, ...) {
  paste0(rxUiGet.monolixExportPath(x, ...), ".txt")
}
#' @export
rxUiGet.monolixDataFile <- function(x, ...) {
  paste0(rxUiGet.monolixExportPath(x, ...), ".csv")
}

#' @export
rxUiGet.monolixQs <- function(x, ...) {
  file.path(rxUiGet.monolixExportPath(x, ...), "nlmixr.qs")
}

#' @export
rxUiGet.monolixRunLock <- function(x, ...) {
  paste0(rxUiGet.monolixExportPath(x, ...), ".run-lock")
}

#' @export
rxUiGet.monolixMlxtranFile <- function(x, ...) {
  paste0(rxUiGet.monolixExportPath(x, ...), ".mlxtran")
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

#' @export
rxUiGet.mlxtran <- function(x, ...) {
  paste(c(rxUiGet.mlxtranDatafile(x, ...),"",
          rxUiGet.mlxtranModel(x, ...),"",
          rxUiGet.mlxtranFit(x, ...),"",
          rxUiGet.mlxtranParameter(x, ...),"",
          rxUiGet.mlxtranMonolix(x, ...)),
        collapse="\n")
}
