#' Does this model have `if`/`else` statements?
#'
#' @param ui rxode2 ui
#' @return boolean saying if the model has `if` or `ifelse()`
#'   statements
#' @noRd
#' @author Matthew L. Fidler
.bblHasIf <- function(ui) {
  .hasIf <- function(x) {
    if (!is.call(x)) return(FALSE)
    if (identical(x[[1]], quote(`if`)) || identical(x[[1]], quote(`ifelse`))) {
      return(TRUE)
    }
    any(vapply(as.list(x)[-1], .hasIf, logical(1), USE.NAMES=FALSE))
  }
  any(vapply(ui$lstExpr, .hasIf, logical(1), USE.NAMES=FALSE))
}

#' Prune the `if`/`else` branches of the model to estimate
#'
#' This uses `rxode2`'s branch pruning to write each `if`/`else`
#' branch as an arithmetic expression, so software that cannot use
#' nested `if`/`else` statements (like NONMEM) can still be used.
#'
#' @param env nlmixr2 estimation environment; `env$ui` is replaced by
#'   the pruned model
#' @param software name of the software for the message
#' @return nothing, called for its side effects
#' @noRd
#' @author Matthew L. Fidler
.bblPruneIf <- function(env, software) {
  .ui <- rxode2::rxUiDecompress(env$ui)
  if (!.bblHasIf(.ui)) return(invisible())
  .env <- new.env(parent=emptyenv())
  .env$.if <- NULL
  .env$.def1 <- NULL
  .pruned <- rxode2::.rxPrune(as.call(c(quote(`{`), .ui$lstExpr)), envir=.env,
                              strAssign=rxode2::rxModelVars(.ui)$strAssign)
  .pruned <- as.list(str2lang(paste0("{", .pruned, "}")))[-1]
  .new <- rxode2::rxUiDecompress(.ui)
  suppressMessages(rxode2::model(.new) <- .pruned)
  .new <- rxode2::rxUiDecompress(.new)
  # keep what nlmixr2 attached to the ui (like the model name used for
  # the output files)
  for (.v in c("modelName", "boundedTransforms")) {
    if (exists(.v, envir=.ui, inherits=FALSE)) {
      assign(.v, get(.v, envir=.ui, inherits=FALSE), envir=.new)
    }
  }
  env$ui <- .new
  .minfo(paste0("pruned if/else branches for ", software))
  invisible()
}

#' Get the prune option from a (possibly incomplete) control
#'
#' @param control control object or list (may be `NULL`)
#' @return boolean, should the model be pruned
#' @noRd
#' @author Matthew L. Fidler
.bblPruneControl <- function(control) {
  if (is.list(control) && !is.null(control$prune)) {
    return(isTRUE(control$prune[1]))
  }
  FALSE
}
