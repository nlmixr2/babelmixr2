#'@export
rxUiGet.nonmemErr <- function(x, ...) {
  paste(c("$ERROR\n",
           "  IPRED = RX_PRED\n",
           "  IRES = DV = IPRED\n",
           "  W     = DSQRT(W)\n",
           "  IF (W .EQ. 0)  W = 1\n",
           "  IWRES = IRES/W\n",
           "  Y = IPRED + EPS(1)*W\n"), collapse="")
}
