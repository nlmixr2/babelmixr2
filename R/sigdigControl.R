# Shared `sigdig` -> tolerance mapping.
#
# These mirror the (unexported) helpers in `nlmixr2est`'s `R/sharedControl.R`
# so that a babelmixr2 estimation method keys its solver and optimizer
# tolerances to `sigdig` exactly the way the nlmixr2est methods do.  Keep them
# in sync with nlmixr2est; when nlmixr2est exports them these can become plain
# re-exports.

#' Optimizer convergence tolerance derived from `sigdig`
#'
#' `10^(-sigdig)` -- the same exponent `sigdig` sets for the ODE `rtol`, so the
#' optimizer converges to exactly the precision the solve can support (no
#' chasing below the solver noise floor).  Mirrors `nlmixr2est:::.sigdigOptTol()`.
#'
#' @param sigdig optimization significant digits
#' @return the tolerance
#' @author Matthew L. Fidler
#' @noRd
.sigdigOptTol <- function(sigdig) 10^(-sigdig)

#' Scale a tuned default tolerance by `sigdig` around `sigdig = 4`
#'
#' Keeps the method's tuned default at `sigdig = 4` and tightens/loosens it by
#' one order of magnitude per significant digit (`default * 10^(4 - sigdig)`).
#' Used where a method's optimizer default should be preserved at the default
#' `sigdig` rather than replaced by the FOCEi formula.  Mirrors
#' `nlmixr2est:::.sigdigScale()`.
#'
#' @param default the tolerance at `sigdig = 4`
#' @param sigdig optimization significant digits
#' @return the scaled tolerance
#' @author Matthew L. Fidler
#' @noRd
.sigdigScale <- function(default, sigdig) default * 10^(4 - sigdig)

#' Scale ODE solver tolerances from the optimization `sigdig`
#'
#' The optimization `sigdig` sets the optimizer tolerances directly
#' (`10^-sigdig`, see [.sigdigOptTol()]); here it sets the ODE solver tolerances
#' with ONE formula -- the same for every solver (stiff, non-stiff,
#' auto-switch):
#'
#'   `rtol = 10^-sigdig`, `atol = 10^(-sigdig-3)`
#'
#' The `rtol` exponent IS `sigdig`, and `atol` sits three orders below.
#' Sensitivity solves match the main solve (the gradient/covariance are built
#' from them); steady-state solves run one order looser.  `tighten` shifts every
#' exponent down by that many orders for a method that needs a tighter solve
#' than the optimizer target.  Mirrors `nlmixr2est:::.rxControlScaleSigdig()`.
#'
#' @param rxControl an [rxode2::rxControl()] object (modified and returned)
#' @param sigdig optimization significant digits; `NULL` leaves `rxControl` as-is
#' @param skip character vector of tolerance field names the user set explicitly
#'   (e.g. from a `rxControl = list(atol = ...)`); these are left untouched so an
#'   explicit `atol`/`rtol` overrides the `sigdig`-derived value
#' @param tighten extra orders of magnitude to tighten every tolerance (default 0)
#' @return `rxControl` with ODE solver tolerances set from `sigdig`
#' @author Matthew L. Fidler
#' @noRd
.rxControlScaleSigdig <- function(rxControl, sigdig, skip = character(0), tighten = 0) {
  if (is.null(sigdig) || is.null(rxControl)) return(rxControl)
  .rtol <- 10^(-(sigdig + tighten))
  .atol <- 10^(-(sigdig + 3 + tighten))
  .set <- function(field, value) {
    if (!(field %in% skip)) rxControl[[field]] <<- rep_len(value, length(rxControl[[field]]))
  }
  .set("rtol", .rtol)
  .set("atol", .atol)
  .set("rtolSens", .rtol)
  .set("atolSens", .atol)
  .set("ssRtol", 10 * .rtol)
  .set("ssAtol", 10 * .atol)
  if (!is.null(rxControl$ssRtolSens)) .set("ssRtolSens", 10 * .rtol)
  if (!is.null(rxControl$ssAtolSens)) .set("ssAtolSens", 10 * .atol)
  rxControl
}

#' Build the `rxControl` for a babelmixr2 estimation control from `sigdig`
#'
#' The shared body every babelmixr2 method's control used to repeat: build (or
#' validate) the solving options and, when `sigdig` is given, key the solver
#' tolerances to it with [.rxControlScaleSigdig()].  A user-supplied
#' `rxControl = list(atol = ...)` wins over the `sigdig`-derived value for the
#' fields it names.
#'
#' @param rxControl `NULL`, an [rxode2::rxControl()] object, or a named list of
#'   solving options
#' @param sigdig optimization significant digits (or `NULL`)
#' @param genRxControl the incoming `genRxControl` flag (from `...`)
#' @return list with `rxControl` (an `rxControl` object) and `genRxControl`
#' @author Matthew L. Fidler
#' @noRd
.babelmixr2RxControlFromSigdig <- function(rxControl, sigdig, genRxControl=FALSE) {
  if (is.null(rxControl)) {
    if (!is.null(sigdig)) {
      rxControl <- .rxControlScaleSigdig(rxode2::rxControl(sigdig=sigdig), sigdig)
    } else {
      rxControl <- rxode2::rxControl(atol=1e-4, rtol=1e-4)
    }
    genRxControl <- TRUE
  } else if (inherits(rxControl, "rxControl")) {
  } else if (is.list(rxControl)) {
    rxControl <- .rxControlScaleSigdig(do.call(rxode2::rxControl, rxControl), sigdig,
                                       skip=names(rxControl))
  } else {
    stop("solving options 'rxControl' needs to be generated from 'rxode2::rxControl'",
         call.=FALSE)
  }
  list(rxControl=rxControl, genRxControl=genRxControl)
}

#' Resolve `sigdigTable` from `sigdig`
#'
#' `sigdigTable` follows `sigdig` (rounded) unless the user gives it, falling
#' back to 3 when there is no `sigdig`.  Matches `nlmixr2est::foceiControl()`.
#'
#' @param sigdigTable the user's `sigdigTable` (possibly `NULL`)
#' @param sigdig optimization significant digits (or `NULL`)
#' @return the resolved `sigdigTable`
#' @author Matthew L. Fidler
#' @noRd
.babelmixr2SigdigTable <- function(sigdigTable, sigdig) {
  if (!is.null(sigdig)) {
    checkmate::assertNumeric(sigdig, lower=1, finite=TRUE, any.missing=TRUE, len=1)
    if (is.null(sigdigTable)) {
      sigdigTable <- round(sigdig)
    }
  }
  if (is.null(sigdigTable)) {
    sigdigTable <- 3
  }
  checkmate::assertIntegerish(sigdigTable, lower=1, len=1, any.missing=FALSE)
  sigdigTable
}
