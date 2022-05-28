#' @export
rxUiGet.nonmemTable <- function(x, ...) {
  paste(c("$TABLE ID ETAS(1:LAST) FIRSTONLY ONEHEADER NOPRINT\n",
          "    NOAPPEND FILE=",
          rxUiGet.nonmemEtaTableName(x, ...),"\n\n",
          "$TABLE ID TIME IPRED ONEHEADER NOPRINT\n",
          "    NOAPPEND FILE=",
          rxUiGet.nonmemSdTableName(x, ...),"\n"
          ), collapse="")
}
