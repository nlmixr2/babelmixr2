#' @export
rxUiGet.nonmemMod <- function(x, ...) {
  .ui <- x[[1]]
  .state <- rxode2::rxModelVars(.ui)$state
  paste(c(paste0("$MODEL NCOMPARTMENTS=", length(.state)),
          vapply(.state,
               function(s){
                 paste0("     COMP(", .nmGetVar(s, .ui),
                        ifelse(s == .state[1], ", DEFDOSE", ""), ") ; ",
                        s)
               }, character(1), USE.NAMES=FALSE)),
        collapse="\n")
}

#' @export
rxUiGet.nonmemModel <- function(x, ...) {
  .ui <- x[[1]]
  .ret <- paste0(
    "$PROBLEM ", .ui$modelName, " translated from babelmixr2\n\n",
    "$DATA ", .ui$modelName, ".csv IGNORE=@\n\n",
    rxUiGet.nonmemInput(x, ...), "\n",
    rxUiGet.nonmemSub(x, ...), "\n\n",
    rxUiGet.nonmemMod(x, ...), "\n\n",
    rxUiGet.nonmemPkDes(x, ...),"\n\n",
    rxUiGet.nonmemErr(x, ...),"\n",
    rxUiGet.nonmemTheta(x, ...),"\n\n",
    rxUiGet.nonmemOmega(x, ...),"\n",
    "$SIGMA 1 FIX\n\n",
    rxUiGet.nonmemEst(x, ...),"\n",
    rxUiGet.nonmemCov(x, ...), "\n\n",
    rxUiGet.nonmemTable(x, ...))
  .ret <- gsub("^ *$", "", .ret)
  .ret
}
