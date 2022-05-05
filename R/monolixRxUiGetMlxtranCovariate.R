#' @export
rxUiGet.mlxtranModelCovariate <- function(x, ...) {
  .ui <- x[[1]]
  .both <- .ui$saemInParsAndMuRefCovariates
  if (length(.both$covars) > 0) {
    return(paste0("[COVARIATE]\ninput = {",
                  paste(.both$covars, collapse=", "), "}\n\n"))
  }
  ""
}
