#' @export
rxUiGet.mlxtranFit <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  .prd <- paste(paste0("rx_prd_", .monolixVar(.predDf$var)), collapse = ", ")
  if (length(.predDf$cond) == 1L) {
    .data <- paste0("data={", .prd, "}")
  } else {
    .data <- paste0("data={", paste(paste0("y", .predDf$dvid), collapse=", "), "}")
  }
  .model <- paste0("model={", .prd, "}")
  paste(c("<FIT>", .data, .model), collapse="\n")
}
attr(rxUiGet.mlxtranFit, "rstudio") <- "character"
