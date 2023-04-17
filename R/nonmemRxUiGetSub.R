#' @export
rxUiGet.nonmemSub <- function(x, ...) {
  .ui <- x[[1]]
  .state <- rxode2::rxModelVars(.ui)$state
  if (length(.state) == 0) return("")
  .advan <- toupper(rxode2::rxGetControl(.ui, "advanOde", "advan13"))
  .ret <- paste0("$SUBROUTINES ", .advan, " TOL=", rxode2::rxGetControl(.ui, "tol", 6),
                 " ATOL=", rxode2::rxGetControl(.ui, "atol", 12),
                 " SSTOL=", rxode2::rxGetControl(.ui, "sstol", 6),
                 " SSATOL=", rxode2::rxGetControl(.ui, "ssatol", 12))
  .c <- rxUiGet.nonmemContra(x, ...)
  if (!is.null(.c)) {
    .ret <- paste0(.ret, "\n  ",
                   "CONTR=",
                   rxUiGet.nonmemContraName(x, ...),
                   " CCONTR=",
                   rxUiGet.nonmemCcontraName(x, ...))
  }
  .ret
}
