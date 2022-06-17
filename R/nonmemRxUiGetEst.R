#' @export
rxUiGet.nonmemEst <- function(x, ...) {
  .ui <- x[[1]]
  .est <- rxode2::rxGetControl(.ui, "est", "focei")
  if (.est == "focei") {
    paste0("$ESTIMATION METHOD=1 INTER MAXEVALS=",
           sprintf("%d", rxode2::rxGetControl(.ui, "maxeval", 10000)),
           " SIGDIG=",sprintf("%d", rxode2::rxGetControl(.ui, "sigdig", 3)),
           " SIGL=",sprintf("%d", rxode2::rxGetControl(.ui, "sigl", 12)),
           " PRINT=", sprintf("%d", rxode2::rxGetControl(.ui, "print", 1)),
           ifelse(rxode2::rxGetControl(.ui, "noabort", TRUE),
                  " NOABORT", ""), "\n")
  } else if (.est == "posthoc") {
    paste0("$ESTIMATION METHOD=0 MAXEVALS=0 POSTHOC PRINT=", rxode2::rxGetControl(.ui, "print", 1),
           " SIGDIG=",sprintf("%d", rxode2::rxGetControl(.ui, "sigdig", 10000)),
           " SIGL=",sprintf("%d", rxode2::rxGetControl(.ui, "sigl", 10000)),
           ifelse(rxode2::rxGetControl(.ui, "noabort", TRUE),
                  " NOABORT", ""), "\n")
  } else if (.est == "imp") {
    paste0("$ESTIMATION METHOD=IMP INTERACTION PRINT=", rxode2::rxGetControl(.ui, "print", 1),
           " SEED=",sprintf("%d", rxode2::rxGetControl(.ui, "seed", 14456)),
           " NITER=",sprintf("%d", rxode2::rxGetControl(.ui, "niter", 100)),
           " ISAMPLE=",sprintf("%d", rxode2::rxGetControl(.ui, "isample", 1000)),
           "\n  ",
           " IACCEPT=",sprintf("%f", rxode2::rxGetControl(.ui, "iaccept", 0.4)),
           " ISCALE_MIN=", sprintf("%f", rxode2::rxGetControl(.ui, "iscaleMin", 0.1)),
           " ISCALE_MAX=", sprintf("%f", rxode2::rxGetControl(.ui, "iscaleMax", 10.0)),
           "\n  ",
           " DF=", sprintf("%f", rxode2::rxGetControl(.ui, "df", 4)),
           " MAPITER=", sprintf("%f", rxode2::rxGetControl(.ui, "mapiter", 1)),
           ifelse(rxode2::rxGetControl(.ui, "noabort", TRUE),
                  " NOABORT", ""), "\n")
  } else {
    stop("unsupported NONMEM est method",
         call.=FALSE)
  }
}
