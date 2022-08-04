#' @export
rxUiGet.nonmemMod <- function(x, ...) {
  .ui <- x[[1]]
  .state <- rxode2::rxModelVars(.ui)$state
  paste(c(paste0("$MODEL NCOMPARTMENTS=", length(.state)),
          vapply(.state,
               function(s) {
                 paste0("     COMP(", .nmGetVar(s, .ui),
                        ifelse(s == .state[1], ", DEFDOSE", ""), ") ; ",
                        s)
               }, character(1), USE.NAMES=FALSE)),
        collapse="\n")
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

#' @export
rxUiGet.nonmemModel <- function(x, ...) {
  .ui <- x[[1]]
  .nonmemResetUi(.ui)
  .ret <- paste0(
    "$PROBLEM ", .ui$nonmemNodelName, " translated from babelmixr2\n\n",
    "$DATA ", .ui$nonmemCsv, " IGNORE=@\n\n",
    rxUiGet.nonmemInput(x, ...), "\n",
    rxUiGet.nonmemSub(x, ...), "\n\n",
    rxUiGet.nonmemMod(x, ...), "\n\n",
    rxUiGet.nonmemPkDesErr0(x, ...),
    rxUiGet.nonmemErrF(x, ...),"\n",
    rxUiGet.nonmemTheta(x, ...),"\n\n",
    rxUiGet.nonmemOmega(x, ...),"\n",
    "$SIGMA 1 FIX\n\n",
    rxUiGet.nonmemEst(x, ...),"\n",
    rxUiGet.nonmemCov(x, ...), "\n\n",
    rxUiGet.nonmemTable(x, ...))
  .ret <- gsub("^ *$", "", .ret)
  .ret
}
