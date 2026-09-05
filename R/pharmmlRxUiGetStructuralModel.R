#' Symbols that the ParameterModel owns
#'
#' These are the mu-referenced individual parameters, i.e. the left-hand sides
#' that `$getSplitMuModel` dropped out of the model body.
#'
#' @param ui rxode2 UI
#' @return character vector of individual-parameter names
#' @noRd
.pharmmlParameterSymbols <- function(ui) {
  .split <- ui$getSplitMuModel
  setNames(c(.split$pureMuRef, .split$taintMuRef), NULL)
}

#' Build the symbol -> block lookup used to qualify SymbRefs
#'
#' PharmML requires a `blkIdRef` when a symbol is defined in a different block
#' from the one referencing it.  Individual parameters live in `pm1` and
#' covariates in `cm1`; everything else (states, local assignments, the
#' independent variable) is local to the structural model and takes no
#' `blkIdRef`.
#'
#' @param ui rxode2 UI
#' @return named character vector mapping symbol name to block id
#' @noRd
.pharmmlBlockMap <- function(ui) {
  .ret <- character(0)
  .pars <- .pharmmlParameterSymbols(ui)
  if (length(.pars) > 0L) {
    .ret <- c(.ret, setNames(rep(.pmlBlk[["parameter"]], length(.pars)), .pars))
  }
  .covs <- ui$allCovs
  if (length(.covs) > 0L) {
    .ret <- c(.ret, setNames(rep(.pmlBlk[["covariate"]], length(.covs)), .covs))
  }
  .ret
}

#' Statements of the model body, with mu-reference and error lines removed
#'
#' @param ui rxode2 UI
#' @return list of R language objects
#' @noRd
.pharmmlModelStatements <- function(ui) {
  .lst <- ui$getSplitMuModel$modelWithDrop
  Filter(function(.e) {
    # mu-referenced parameter definitions are replaced by the bare symbol
    # `_drop`; residual-error lines are `~` calls and belong to the
    # ObservationModel, not here.
    if (is.name(.e)) return(FALSE)
    if (is.call(.e) && identical(as.character(.e[[1]]), "~")) return(FALSE)
    TRUE
  }, .lst)
}

#' Is this statement an explicit initial condition, e.g. `center(0) <- 100`?
#'
#' @param e R language object
#' @return TRUE when the statement sets an initial condition
#' @noRd
.pharmmlIsInitialCondition <- function(e) {
  if (!is.call(e)) return(FALSE)
  if (!identical(as.character(e[[1]]), "<-")) return(FALSE)
  .lhs <- e[[2]]
  if (!is.call(.lhs)) return(FALSE)
  # a call whose function is a plain name and whose single argument is 0
  is.name(.lhs[[1]]) && length(.lhs) == 2L && identical(.lhs[[2]], 0)
}

#' Is this statement a derivative, e.g. `d/dt(center) <- ...`?
#'
#' @param e R language object
#' @return TRUE when the statement defines a derivative
#' @noRd
.pharmmlIsDerivative <- function(e) {
  if (!is.call(e)) return(FALSE)
  if (!identical(as.character(e[[1]]), "<-")) return(FALSE)
  .lhs <- e[[2]]
  # `d/dt(x)` parses as `/`(d, dt(x)) -- it is a division call, not a single
  # `d/dt` name.
  if (!is.call(.lhs) || !identical(as.character(.lhs[[1]]), "/")) return(FALSE)
  if (length(.lhs) != 3L) return(FALSE)
  if (!is.name(.lhs[[2]]) || !identical(as.character(.lhs[[2]]), "d")) return(FALSE)
  .rhs <- .lhs[[3]]
  is.call(.rhs) && identical(as.character(.rhs[[1]]), "dt") && length(.rhs) == 2L
}

#' The state name a `d/dt()` statement assigns to
#'
#' @param e derivative statement
#' @return state name
#' @noRd
.pharmmlDerivativeState <- function(e) {
  as.character(e[[2]][[3]][[2]])
}

#' Collect explicit initial conditions, keyed by state name
#'
#' @param ui rxode2 UI
#' @return named list of R language objects holding the initial value
#' @noRd
.pharmmlInitialConditions <- function(ui) {
  .ret <- list()
  for (.e in .pharmmlModelStatements(ui)) {
    if (.pharmmlIsInitialCondition(.e)) {
      .ret[[as.character(.e[[2]][[1]])]] <- .e[[3]]
    }
  }
  .ret
}

#' Emit a ct:DerivativeVariable
#'
#' @param state state name
#' @param rhs right-hand side, as an R language object
#' @param init initial value, as an R language object
#' @param ui rxode2 UI
#' @return character(1)
#' @noRd
.pharmmlDerivativeVariable <- function(state, rhs, init, ui) {
  .pmlNode(
    "ct:DerivativeVariable",
    attrs = c(symbId = state, symbolType = "real"),
    children = c(
      .pharmmlAssign(.rxToPharmml(rhs, ui)),
      .pmlNode("ct:IndependentVariable",
               children = .pmlNode("ct:SymbRef", attrs = c(symbIdRef = "t"))),
      .pmlNode("ct:InitialCondition",
               children = c(
                 .pmlNode("ct:InitialValue",
                          children = .pharmmlAssign(.rxToPharmml(init, ui))),
                 .pmlNode("ct:InitialTime",
                          children = .pharmmlAssign(.pmlText("ct:Real", 0)))))))
}

#' PharmML StructuralModel block
#'
#' @param ui rxode2 UI
#'
#' @param indent Indent depth
#'
#' @return character(1) holding the `StructuralModel` block
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.pharmmlStructuralModel <- function(ui, indent = 0L) {
  if (.pharmmlIsLinCmt(ui)) {
    # A solved model is expressed as PK macros rather than derivatives, which
    # keeps the structure the model was written in.
    return(.pmlNode("mdef:StructuralModel",
                    attrs = c(blkId = .pmlBlk[["structural"]]),
                    children = .pharmmlPkMacros(ui),
                    indent = indent))
  }
  .init <- .pharmmlInitialConditions(ui)
  .children <- character(0)

  for (.e in .pharmmlModelStatements(ui)) {
    if (.pharmmlIsInitialCondition(.e)) next # folded into the derivative
    if (.pharmmlIsDerivative(.e)) {
      .state <- .pharmmlDerivativeState(.e)
      .i <- .init[[.state]]
      if (is.null(.i)) .i <- 0
      .children <- c(.children,
                     .pharmmlDerivativeVariable(.state, .e[[3]], .i, ui))
      next
    }
    if (is.call(.e) && identical(as.character(.e[[1]]), "<-") && is.name(.e[[2]])) {
      .children <- c(.children,
                     .pmlNode("ct:Variable",
                              attrs = c(symbId = as.character(.e[[2]]),
                                        symbolType = "real"),
                              children = .pharmmlAssign(.rxToPharmml(.e[[3]], ui))))
      next
    }
    stop("cannot translate the model statement '", deparse1(.e),
         "' to PharmML", call. = FALSE)
  }

  .pmlNode("mdef:StructuralModel",
           attrs = c(blkId = .pmlBlk[["structural"]]),
           children = .children,
           indent = indent)
}

#' @export
rxUiGet.pharmmlStructuralModel <- function(x, ...) {
  .pharmmlStructuralModel(x[[1]])
}
attr(rxUiGet.pharmmlStructuralModel, "rstudio") <- "character"
