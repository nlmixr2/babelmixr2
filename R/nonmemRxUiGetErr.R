#'@export
rxUiGet.nonmemErr <- function(x, ...) {
  paste0(c("$ERROR\n",
           "  IPRED = RX_PRED\n",
           "  IRES = DV = IPRED\n",
           "  W     = DSQRT(W)\n",
           "  IF (W .EQ. 0)  W = 1\n",
           "  IWRES = IRES/W\n",
           "  Y = IPRED + EPS(1)*W\n\n",
           "$SIGMA 1 FIX\n"))
}
