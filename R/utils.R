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
    } else if (is.character(nm)) {
      checkmate::assertCharacter(nm, len=1, any.missing=FALSE)
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

#' Determine if the expression is d/dt() 
#'
#' @param expr expression
#' @return boolean determining if expression is a d/dt() expression
#' @author Matthew L. Fidler
#' @noRd
.rxIsDdt <- function(expr) {
  if (length(expr) != 3) return(FALSE)
  if (!identical(expr[[1]], quote(`/`))) return(FALSE)
  if (!identical(expr[[2]], quote(`d`))) return(FALSE)
  if (length(expr[[3]]) != 2) return(FALSE)
  if (!identical(expr[[3]][[1]], quote(`dt`))) return(FALSE)
  TRUE
}
#' Is a variable known to be non-zero
#'
#' @param variable Varible name to check
#' @param ui rxode2 ui
#' @return logical to say if the variable is known to be non-zero
#' @author Matthew L. Fidler
#' @noRd
.rxIsKnownNonZeroVariable <- function(variable, ui) {
  if (!is.null(names(variable))) {
    if (.rxIsKnownNonZeroVariable(names(variable), ui)) return(TRUE)
  }
  variable <- as.character(variable)
  if (variable %in% ui$allCovs) return(TRUE)
  .split <- ui$getSplitMuModel
  .mus <- c(.split$pureMuRef, .split$taintMuRef)
  .w <- which(variable == .mus)
  if (length(.w) != 1) return(FALSE)
  .tv <- names(.w)
  .w <- which(ui$muRefCurEval$parameter == .tv)
  if (length(.w) != 1) return(FALSE)
  .curEval <- ui$muRefCurEval$parameter$curEval[.w]
  if (.curEval == "exp") return(TRUE)
  if (any(.curEval == c("expit", "probitInv"))) {
    .low <- ui$muRefCurEval$parameter$low[.w]
    if (.low >= 0) return(TRUE)
    .hi <- ui$muRefCurEval$parameter$hi[.w]
    if (.hi <= 0) return(TRUE)
  }
  FALSE
}
#' Should the variable be protected from being zero?
#'
#' @param variable name of the variable
#' @param ui rxode2 ui
#' @return boolean of if the variable needs protection
#' @author Matthew L. Fidler
#' @noRd
.rxShouldProtectZeros <- function(variable, ui) {
  # should protect zeros if requested, not in an if/else block
  # and if the variable is known to be something non-zero
  if (!rxode2::rxGetControl(ui, "protectZeros", TRUE)) return(FALSE)
  if (rxode2::rxGetControl(ui, ".ifelse", FALSE)) return(FALSE)
  if (.rxIsKnownNonZeroVariable(variable, ui)) return(FALSE)
  return(TRUE)
}
