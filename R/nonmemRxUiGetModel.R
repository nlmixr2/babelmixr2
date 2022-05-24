#' @export
rxUiGet.nonmemModel <- function(x, ...) {
  paste0(rxUiGet.nonmemPkDes(x, ...),"\n",
         rxUiGet.nonmemErr(x, ...),"\n",
         rxUiGet.nonmemTheta(x, ...),"\n\n",
         rxUiGet.nonmemOmega(x, ...),"\n",
         "$SIGMA 1 FIX\n\n",
         rxUiGet.nonmemEst(x, ...),"\n",
         rxUiGet.nonmemCov(x, ...), "\n\n",
         rxUiGet.nonmemTable(x, ...))
}
