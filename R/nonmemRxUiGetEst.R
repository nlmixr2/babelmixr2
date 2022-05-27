#' @export
rxUiGet.nonmemEst <- function(x, ...) {
  .ui <- x[[1]]
  .est <- rxode2::rxGetControl(.ui, "est", "focei")
  if (.est == "focei") {
    paste0("$ESTIMATION METHOD=1 INTER MAXEVALS=",
           sprintf("%d", rxode2::rxGetControl(.ui, "maxeval", 10000)),
           " PRINT=", sprintf("%d", rxode2::rxGetControl(.ui, "print", 1)),
           ifelse(rxode2::rxGetControl(.ui, "noabort", TRUE),
                  " NOABORT", ""), "\n")
  } else if (.est == "posthoc") {
    paste0("$ESTIMATION METHOD=0 MAXEVALS=0 POSTHOC PRINT=", rxode2::rxGetControl(.ui, "print", 1),
           ifelse(rxode2::rxGetControl(.ui, "noabort", TRUE),
                  " NOABORT", ""), "\n")
  } else {
    stop("unsupported NONMEM est method",
         call.=FALSE)
  }
}
