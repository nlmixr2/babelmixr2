#' @export
rxUiGet.nonmemSub <- function(x, ...) {
  .ui <- x[[1]]
  .advan <- toupper(rxode2::rxGetControl(.ui, "advanOde", "advan13"))
  .ret <- paste0("$SUBROUTINES ", .advan, " TOL=", rxode2::rxGetControl(.ui, "tol", 6),
                 " SIGL=", rxode2::rxGetControl(.ui, "sigl", 6))
  .c <- rxUiGet.nonmemContr(x, ...)
  if (!is.null(.c)) {
    .ret <- paste0(.ret, "\n  ",
                   "CONTR=",
                   rxUiGet.nonmemContraName(x, ...),
                   " CCONTR=",
                   rxUiGet.nonmemCcontraName(x, ...))
  }
  .ret
}