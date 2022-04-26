
.getPureMuRef <- function(expr, muRefCurEval) {
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
      .char <- as.character(expr[[2]])
      .callName <- as.character(.call[[1]])
      .w <- which(muRefCurEval$parameter == as.character(.call[[2]]))
      if (length(.w) == 1L) {
        if (muRefCurEval$curEval[.w] == .callName) {
          # This is an additive mu reference expression
          # c(tcl="cl")
          .low <- muRefCurEval$low[.w]
          .hi <- muRefCurEval$hi[.w]
          if (length(.call) == 2) {
            return(setNames(.char, as.character(.call[[2]])))
          } else if (length(.call) == 3 && .call[[3]] == ifelse(is.na(.low), 0, .low)) {
            return(setNames(.char, as.character(.call[[2]])))
          } else if (length(.call) == 4 &&
                       .call[[3]] == ifelse(is.na(.low), 0, .low) &&
                       .call[[4]] == ifelse(is.na(.hi), 1, .hi) &&
                       TRUE) {
            return(setNames(.char, as.character(.call[[2]])))
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
