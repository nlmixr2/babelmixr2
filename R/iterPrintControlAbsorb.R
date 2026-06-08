#' Wrap scalar or list arguments into an iterPrintControl object
#'
#' Internal helper used by `pseudoOptimControl()` and `fmeMcmcControl()`
#' to absorb the scalar `print` / `printNcol` / `useColor` arguments
#' into a single [nlmixr2est::iterPrintControl()] sub-list.  If the
#' user already passed a pre-built `iterPrintControl()` object via the
#' `print` argument (or, on round-trip, via an `iterPrintControl =`
#' slot in `...`), return it directly.
#'
#' Mirrors `nlmixr2est:::.absorbIterPrintControl()`; kept local so we
#' do not depend on a non-exported helper.
#'
#' @param print Either an integer print-frequency or an
#'   `iterPrintControl` object.
#' @param printNcol,useColor Scalar `*Control()` arguments forwarded
#'   to [nlmixr2est::iterPrintControl()] only when `print` is a scalar.
#' @param iterPrintControl Optional pre-built
#'   [nlmixr2est::iterPrintControl()] object.  Wins over `print` and
#'   the other scalars when supplied — used by the round-trip case
#'   where a returned control list is passed back through
#'   `do.call(*Control, .ctl)`.
#' @return An `iterPrintControl` list.
#' @noRd
.absorbIterPrintControl <- function(print = 1L,
                                    printNcol = NULL,
                                    useColor = NULL,
                                    iterPrintControl = NULL) {
  if (!is.null(iterPrintControl)) {
    if (!inherits(iterPrintControl, "iterPrintControl")) {
      stop("`iterPrintControl` must be the result of iterPrintControl()",
           call. = FALSE)
    }
    return(iterPrintControl)
  }
  if (inherits(print, "iterPrintControl")) {
    .conflicts <- character(0)
    if (!is.null(printNcol)) .conflicts <- c(.conflicts, "printNcol")
    if (!is.null(useColor))  .conflicts <- c(.conflicts, "useColor")
    if (length(.conflicts)) {
      warning("ignoring `", paste(.conflicts, collapse = "`, `"),
              "` because `print` was passed as an iterPrintControl object",
              call. = FALSE)
    }
    return(print)
  }
  nlmixr2est::iterPrintControl(every = print,
                               ncol = printNcol,
                               useColor = useColor)
}
