#' @export
rxUiGet.nonmemTable <- function(x, ...) {
  paste(c("$TABLE ID ETAS(1:LAST) OBJI FIRSTONLY ONEHEADER NOPRINT\n",
          "     FORMAT=s1PE17.9 NOAPPEND FILE=",
          rxUiGet.nonmemEtaTableName(x, ...),"\n\n",
          "$TABLE ID TIME IPRED PRED RXROW ONEHEADER NOPRINT\n",
          "    FORMAT=s1PE17.9 NOAPPEND FILE=",
          rxUiGet.nonmemSdTableName(x, ...),"\n"
          ), collapse="")
}
