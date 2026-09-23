#' @export
rxUiGet.nonmemTable <- function(x, ...) {
  .etas <- "LAST"
  if (!is.null(rxUiGet.nonmemPriorSpec(x, ...))) {
    # NWPRI's own prior variances are OMEGAs after the model's, so
    # LAST would output etas that are not in the model
    .iniDf <- x[[1]]$iniDf
    .etas <- paste(sum(!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2))
  }
  paste(c("$TABLE ID ETAS(1:", .etas, ") OBJI FIRSTONLY ONEHEADER NOPRINT\n",
          "     FORMAT=s1PE17.9 NOAPPEND FILE=",
          rxUiGet.nonmemEtaTableName(x, ...),"\n\n",
          "$TABLE ID TIME IPRED PRED RXROW ONEHEADER NOPRINT\n",
          "    FORMAT=s1PE17.9 NOAPPEND FILE=",
          rxUiGet.nonmemSdTableName(x, ...),"\n"
          ), collapse="")
}
attr(rxUiGet.nonmemTable, "rstudio") <- "nonmemTable"
