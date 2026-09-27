#' @export
rxUiGet.nonmemMod <- function(x, ...) {
  .ui <- x[[1]]
  .state <- rxode2::rxModelVars(.ui)$state
  if (length(.state) == 0) return("")
  # closed-form ADVANs have their own compartments
  if (!is.null(.nonmemLinCmtAdvan(.ui))) return("")
  paste(c(paste0("$MODEL NCOMPARTMENTS=", length(.state)),
          vapply(.state,
               function(s) {
                 paste0("     COMP(", .nmGetVar(s, .ui),
                        ifelse(s == .state[1], ", DEFDOSE", ""), ") ; ",
                        s)
               }, character(1), USE.NAMES=FALSE)),
        collapse="\n")
}
attr(rxUiGet.nonmemMod, "rstudio") <- "nonmemMod"

#' $MODEL section followed by its spacing
#'
#' A closed-form ADVAN has no $MODEL, so it has no spacing either
#'
#' @inheritParams rxUiGet.nonmemMod
#' @return $MODEL text for the control stream
#' @noRd
#' @author Matthew L. Fidler
.nonmemModSection <- function(x, ...) {
  .mod <- rxUiGet.nonmemMod(x, ...)
  if (.mod == "" && !is.null(.nonmemLinCmtAdvan(x[[1]]))) return("")
  paste0(.mod, "\n\n")
}

.nonmemResetUi <- function(ui, extra="") {
  rxode2::rxAssignControlValue(ui, ".nmGetDivideZeroDf",
                               data.frame(expr=character(0),
                                          nm=character(0)))
  rxode2::rxAssignControlValue(ui, ".nmVarNum", 1)
  rxode2::rxAssignControlValue(ui, ".nmGetVarDf",
                               data.frame(var=character(0),
                                          nm=character(0)))

  rxode2::rxAssignControlValue(ui, ".nmVarDZNum", 1)
  rxode2::rxAssignControlValue(ui, ".nmGetDivideZeroDf",
                               data.frame(expr=character(0),
                                          nm=character(0)))
  rxode2::rxAssignControlValue(ui, ".nmPrefixLines", NULL)
  rxode2::rxAssignControlValue(ui, ".nmVarExtra", extra)
}

rxUiGetNonememModelEnv <- new.env(parent=emptyenv())
rxUiGetNonememModelEnv$rxS <- NULL

#' @export
rxUiGet.nonmemModel <- function(x, ...) {
  .ui <- x[[1]]
  rxUiGetNonememModelEnv$rxS <- .ui$loadPrune
  .nonmemResetUi(.ui)
  if (!is.null(.nonmemTnpri(.ui))) {
    return(.nonmemTnpriModel(x, ...))
  }
  .ret <- paste0(
    "$PROBLEM ", .ui$nonmemModelName, " translated from babelmixr2\n; comments show mu referenced model in ui$getSplitMuModel\n\n",
    "$DATA ", .ui$nonmemCsv, " IGNORE=@\n\n",
    rxUiGet.nonmemInput(x, ...), "\n",
    rxUiGet.nonmemSub(x, ...), "\n\n",
    rxUiGet.nonmemPrior(x, ...),
    .nonmemModSection(x, ...),
    rxUiGet.nonmemPkDesErr0(x, ...),
    rxUiGet.nonmemErrF(x, ...),"\n",
    rxUiGet.nonmemTheta(x, ...),"\n\n",
    rxUiGet.nonmemOmega(x, ...),"\n",
    "$SIGMA 1 FIX\n\n",
    rxUiGet.nonmemPriorRecords(x, ...),
    rxUiGet.nonmemEst(x, ...),"\n",
    rxUiGet.nonmemCov(x, ...), "\n\n",
    rxUiGet.nonmemTable(x, ...))
  .ret <- gsub("^ *$", "", .ret)
  .ret
}
attr(rxUiGet.nonmemModel, "rstudio") <- "nonmemModel"

#' The two problem control stream of a TNPRI fit
#'
#' Problem 1 reads the prior from the model specification file of the
#' earlier fit and holds the model code, which NONMEM uses for every
#' problem of the run; problem 2 estimates the model with the prior.
#'
#' @inheritParams rxUiGet.nonmemMod
#' @return control stream
#' @noRd
#' @author Matthew L. Fidler
.nonmemTnpriModel <- function(x, ...) {
  .ui <- x[[1]]
  .tnpri <- .nonmemTnpri(.ui)
  .msf <- .nonmemTnpriMsfName(.ui)
  .input <- rxUiGet.nonmemInput(x, ...)
  .ret <- paste0(
    "$PROBLEM ", .ui$nonmemModelName, " TNPRI prior from ", .msf, "\n",
    "; the prior is the estimates and covariance of an earlier fit of this model;\n",
    "; the model code in this problem is used by every problem\n\n",
    "$DATA ", .ui$nonmemCsv, " IGNORE=@\n\n",
    .input, "\n",
    rxUiGet.nonmemSub(x, ...), "\n\n",
    "$PRIOR TNPRI ", .nonmemTnpriOptions(.tnpri), "\n\n",
    "$MSFI ", .msf, " ONLYREAD\n\n",
    .nonmemModSection(x, ...),
    rxUiGet.nonmemPkDesErr0(x, ...),
    rxUiGet.nonmemErrF(x, ...), "\n",
    "$PROBLEM ", .ui$nonmemModelName, " translated from babelmixr2\n",
    "; comments show mu referenced model in ui$getSplitMuModel\n\n",
    "$DATA ", .ui$nonmemCsv, " IGNORE=@ REWIND\n\n",
    .input, "\n",
    rxUiGet.nonmemTheta(x, ...), "\n\n",
    rxUiGet.nonmemOmega(x, ...), "\n",
    "$SIGMA 1 FIX\n\n",
    rxUiGet.nonmemEst(x, ...), "\n",
    rxUiGet.nonmemCov(x, ...), "\n\n",
    rxUiGet.nonmemTable(x, ...))
  gsub("^ *$", "", .ret)
}
