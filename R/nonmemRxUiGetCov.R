#'@export
rxUiGet.nonmemCov <- function(x, ...) {
  .ui <- x[[1]]
  .est <- rxode2::rxGetControl(.ui, "cov", "r,s")
  if (.est == "") {
    return("")
  } else if (.est == "r,s") {
    return("$COVARIANCE")
  } else if (.est == "r") {
    return("$COVARIANCE MATRIX=R")
  } else {
    return("$COVARIANCE MATRIX=S")
  }
}
