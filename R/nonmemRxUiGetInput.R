#' @export
rxUiGet.nonmemInput <- function(x, ...) {
  .ui <- x[[1]]
  .ret <- c("$INPUT ID TIME EVID AMT",
            ifelse(rxode2::rxGetControl(.ui, ".hasIi", FALSE), "II", ""),
            "DV CMT",
            ifelse(length(.ui$predDf$cond) > 1, "DVID", ""),
            ifelse(rxode2::rxGetControl(.ui, ".hasSs", FALSE), "SS", ""),
            ifelse(rxode2::rxGetControl(.ui, ".hasRate", FALSE), "RATE", ""),
            vapply(.ui$allCovs, .nmGetVar, character(1), ui=.ui,
                   USE.NAMES=FALSE),
            ifelse(rxode2::rxGetControl(.ui, ".hasCens", FALSE), "CENS", ""),
            ifelse(rxode2::rxGetControl(.ui, ".hasLimit", FALSE), "LIMIT", ""),
            "DROP")
  .ret <- .ret[.ret != ""]
  paste0(paste(.ret, collapse=" "), "\n")
}
