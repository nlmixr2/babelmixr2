
.isPureMuRef <- function(expr, muRefCurEval) {
  # probitInv
  # expit
  # logit
  # probit
  # exp
  # ""
  if (identical(expr[[1]], quote(`=`)) ||
        identical(expr[[1]], quote(`<-`))) {
    if (is.call(expr[[3]])) {
      .call <- expr[[3]]
      .char <- as.character(.call[[2]])
      .callName <- as.character(.call[[1]])
      .w <- which(muRefCurEval$parameter == .char)
      if (length(.w) == 1L) {
        if (muRefCurEval$curEval[.w] == .callName) {
          # This is an additive mu reference expression
          # c(tcl="cl")
          .low <- muRefCurEval$low[.w]
          .hi <- muRefCurEval$hi[.w]
          if (is.na(.low) && length(.call) == 2) {
            return(setNames(as.character(.call[[2]]), .char))
          } else if (!is.na(.low) && is.na(.hi) && length(.call) == 3 && .call[[3]] == .low) {
            return(setNames(as.character(.call[[2]]), .char))
          } else if (!is.na(.low) && !is.na(.hi) && length(.call) == 1 && .call[[3]] == .low && .call[[4]] == .hi) {
            return(setNames(as.character(.call[[2]]), .char))
          }
        }
      }
    } else {
      .char <- as.character(expr[[3]])
      .w <- which(muRefCurEval$parameter == .char)
      if (length(.w) == 1L) {
        if (muRefCurEval$curEval[.w] == "") {
          # This is an additive mu reference expression
          # c(cl="tcl")
          return(setNames(as.character(expr[[2]]), .char))
        }
      }
    }
  }
  return(NULL)
}

#' @export
rxUiGet.getSplitMuModel <- function(x, ...) {
  .ui <- x[[1]]
  nlmixr2est::.saemDropMuRefFromModel(.ui)
}
