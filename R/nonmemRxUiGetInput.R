#' @export
rxUiGet.nonmemInput <- function(x, ...) {
  .ui <- x[[1]]
  paste(paste(c("$INPUT ID TIME EVID AMT II DV CMT SS",
         vapply(.ui$allCovs, .nmGetVar, character(1), ui=.ui,
                USE.NAMES=FALSE), "ROW_"),
        collapse=" "), "\n")
}
