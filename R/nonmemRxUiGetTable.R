#' @export
rxUiGet.nonmemTable <- function(x, ...) {
  paste(c("$TABLE ID ETAS(1:LAST) FIRSTONLY ONEHEADER NOPRINT\n",
          "    NOAPPEND FILE=",
          rxUiGet.nonmemEtaTableName(x, ...)), collapse="")
}
