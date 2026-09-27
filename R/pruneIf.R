#' Does this model have `if`/`else` statements?
#'
#' @param ui rxode2 ui
#' @return boolean saying if the model has `if` or `ifelse()`
#'   statements
#' @noRd
#' @author Matthew L. Fidler
.bblHasIf <- function(ui) {
  .hasIf <- function(x) {
    if (!is.call(x)) {
      return(FALSE)
    }
    if (identical(x[[1]], quote(`if`)) || identical(x[[1]], quote(`ifelse`))) {
      return(TRUE)
    }
    any(vapply(as.list(x)[-1], .hasIf, logical(1), USE.NAMES = FALSE))
  }
  any(vapply(ui$lstExpr, .hasIf, logical(1), USE.NAMES = FALSE))
}

#' Does this model have `if`/`else` statements the software cannot write?
#'
#' `ifelse()` always needs the branches pruned.  When the software
#' cannot write nested `if`/`else` statements (like NONMEM, where only
#' simple `IF` blocks are written), an `else`, `else if` or a nested
#' `if` also needs the branches pruned.
#'
#' @param ui rxode2 ui
#' @param nested can the software write nested `if`/`else` statements
#'   (like Monolix)?
#' @return boolean saying if the model needs its branches pruned
#' @noRd
#' @author Matthew L. Fidler
.bblNeedsPrune <- function(ui, nested = FALSE) {
  .needs <- function(x, inIf = FALSE) {
    if (!is.call(x)) {
      return(FALSE)
    }
    if (identical(x[[1]], quote(`ifelse`))) {
      return(TRUE)
    }
    if (identical(x[[1]], quote(`if`))) {
      if (!nested && (inIf || length(x) > 3L)) {
        return(TRUE)
      }
      return(any(vapply(
        as.list(x)[-1],
        .needs,
        logical(1),
        inIf = TRUE,
        USE.NAMES = FALSE
      )))
    }
    any(vapply(
      as.list(x)[-1],
      .needs,
      logical(1),
      inIf = inIf,
      USE.NAMES = FALSE
    ))
  }
  any(vapply(ui$lstExpr, .needs, logical(1), USE.NAMES = FALSE))
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
  if (!.bblHasIf(.ui)) {
    return(invisible())
  }
  .env <- new.env(parent = emptyenv())
  .env$.if <- NULL
  .env$.def1 <- NULL
  .pruned <- rxode2::.rxPrune(
    as.call(c(quote(`{`), .ui$lstExpr)),
    envir = .env,
    strAssign = rxode2::rxModelVars(.ui)$strAssign
  )
  .pruned <- as.list(str2lang(paste0("{", .pruned, "}")))[-1]
  .new <- rxode2::rxUiDecompress(.ui)
  suppressMessages(rxode2::model(.new) <- .pruned)
  .new <- rxode2::rxUiDecompress(.new)
  # keep what nlmixr2 attached to the ui (like the model name used for
  # the output files)
  for (.v in c("modelName", "boundedTransforms")) {
    if (exists(.v, envir = .ui, inherits = FALSE)) {
      assign(.v, get(.v, envir = .ui, inherits = FALSE), envir = .new)
    }
  }
  env$ui <- .new
  .minfo(paste0("pruned if/else branches for ", software))
  invisible()
}

#' Should the model be pruned for the estimation software?
#'
#' @param env nlmixr2 estimation environment with `env$ui` and
#'   `env$control`
#' @param nested can the software write nested `if`/`else` statements?
#' @return boolean, should the model be pruned; `"auto"` (the default)
#'   prunes only when the `if`/`else` statements cannot be written
#'   directly
#' @noRd
#' @author Matthew L. Fidler
.bblPruneControl <- function(env, nested = FALSE) {
  .prune <- "auto"
  .control <- env$control
  if (is.list(.control) && !is.null(.control$prune)) {
    .prune <- .control$prune[1]
  }
  if (isTRUE(.prune)) {
    return(TRUE)
  }
  if (isFALSE(.prune)) {
    return(FALSE)
  }
  .bblNeedsPrune(rxode2::rxUiDecompress(env$ui), nested = nested)
}
