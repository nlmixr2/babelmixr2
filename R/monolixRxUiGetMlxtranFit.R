#' @export
rxUiGet.mlxtranFit <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  if (length(.predDf$cond) == 1L) {
    .data <- paste0("data={", paste(paste0("rx_prd_", .predDf$var), collapse=", "), "}")
  } else {
    .data <- paste0("data={", paste(paste0("y", .predDf$dvid), collapse=", "), "}")
  }
  .model <- paste0("model={", paste(paste0("rx_prd_", .predDf$var), collapse=", "), "}")
  paste(c("<FIT>", .data, .model), collapse="\n")
}
attr(rxUiGet.mlxtranFit, "rstudio") <- "character"
