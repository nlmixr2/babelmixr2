# Constants inlined as literals; PharmML has math:Constant for pi and e but
# emitting a Real keeps the walker uniform and round-trips exactly.
.rxPmlCnt <- c(
  pi        = "3.141592653589793",
  M_E       = "2.718281828459045",
  M_LN2     = "0.6931471805599453",
  M_LN10    = "2.302585092994046",
  M_PI      = "3.141592653589793",
  M_PI_2    = "1.570796326794897",
  M_PI_4    = "0.7853981633974483",
  M_SQRT2   = "1.414213562373095",
  M_LOG2E   = "1.442695040888963",
  M_LOG10E  = "0.4342944819032518"
)

# rxode2 constructs with no PharmML equivalent.
.rxPmlBad <- c("NA", "NaN", "Inf", "newind", "NEWIND")

.rxPmlBadF <- c("linCmt", "mtime", "tad", "digamma", "trigamma", "tetragamma",
                "pentagamma", "psigamma", "choose", "lchoose")

#' Translate an rxode2 expression to a PharmML math tree
#'
#' @param x A `call`, `name` or scalar, as produced by `str2lang()`.
#'
#' @param ui Optional rxode2 UI, used to resolve which block a symbol lives in
#'   so that `blkIdRef` can be emitted.  When `NULL` (the default) bare
#'   `symbIdRef` references are emitted, which is what unit tests want.
#'
#' @param indent Indent depth.
#'
#' @return character(1) holding the emitted PharmML subtree.
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.rxToPharmml <- function(x, ui = NULL, indent = 0L) {
  if (is.numeric(x) || is.integer(x)) return(.rxToPharmmlScalar(x, indent))
  if (is.name(x)) return(.rxToPharmmlSymbol(x, ui, indent))
  if (is.call(x)) return(.rxToPharmmlCall(x, ui, indent))
  stop("cannot translate '", deparse1(x), "' to PharmML", call. = FALSE)
}

#' @noRd
.rxToPharmmlScalar <- function(x, indent = 0L) {
  if (is.integer(x)) return(.pmlText("ct:Int", x, indent))
  .pmlText("ct:Real", x, indent)
}

#' @noRd
.rxToPharmmlSymbol <- function(x, ui = NULL, indent = 0L) {
  .n <- as.character(x)
  if (.n %in% .rxPmlBad) {
    stop("'", .n, "' has no PharmML equivalent", call. = FALSE)
  }
  if (.n %in% names(.rxPmlCnt)) {
    return(.pmlText("ct:Real", .rxPmlCnt[[.n]], indent))
  }
  if (.n == "time" || .n == "t") {
    return(.pmlNode("ct:SymbRef", attrs = c(symbIdRef = "t"), indent = indent))
  }
  .pmlNode("ct:SymbRef",
           attrs = .rxToPharmmlSymbAttrs(.n, ui),
           indent = indent)
}

#' Build the attribute set for a SymbRef, adding blkIdRef when resolvable
#'
#' @param n symbol name
#' @param ui rxode2 UI or NULL
#' @return named character vector
#' @noRd
.rxToPharmmlSymbAttrs <- function(n, ui = NULL) {
  .blk <- .rxToPharmmlBlockOf(n, ui)
  if (is.na(.blk)) return(c(symbIdRef = n))
  c(blkIdRef = .blk, symbIdRef = n)
}

#' Which PharmML block a symbol belongs to
#'
#' @param n symbol name
#' @param ui rxode2 UI or NULL
#' @return block id, or NA_character_ when unresolvable
#' @noRd
.rxToPharmmlBlockOf <- function(n, ui = NULL) {
  # Stage 3 populates this from the ui; with no ui there is nothing to resolve.
  NA_character_
}

#' @noRd
.rxToPharmmlCall <- function(x, ui = NULL, indent = 0L) {
  .fn <- as.character(x[[1]])
  if (.fn %in% .rxPmlBadF) {
    stop("'", .fn, "()' has no PharmML equivalent", call. = FALSE)
  }
  stop("cannot translate '", deparse1(x), "' to PharmML", call. = FALSE)
}
