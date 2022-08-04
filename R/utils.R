#' Get compartment number from name
#'
#' @param nm Name of the compartment
#' @param ui rxode2 ui to query against
#' @param error Should this error if not found (if `FALSE` return `NA_integer`)
#' @return compartment number (or possibily NA)
#' @author Bill Denney and Matthew L. Fidler
#' @noRd
.rxGetCmtNumber <- function(nm, ui, error=TRUE) {
  if (is.numeric(nm)) {
    nm
  } else {
    if (is.name(nm)) {
      nm <- as.character(nm)
    } else if (error) {
      stop("invalid attempt to find a compartment number: ", as.character(nm))
    } else {
      return(NA_integer_)
    }
    .cmt <- which(rxode2::rxModelVars(ui)$state %in% nm)
    if (length(.cmt) == 0) {
      if (error) {
        stop("cannot find compartment '", nm, "' in the model",
             call.=FALSE)
      } else {
        return(NA_integer_)
      }
    }
    .cmt
  }
}
