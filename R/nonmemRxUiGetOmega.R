.nonmemHandleOneOmega <- function(om0, ui) {
  .dim <- dim(om0)
  .dimn1 <- dimnames(om0)[[1]][1]
  .w <- which(ui$iniDf$name == .dimn1)
  .fix <- ui$iniDf$fix[.w]
  .block <- ""
  .sigdig <- rxode2::rxGetControl(ui, "iniSigDig", 5)
  if (.dim[1] > 1L) {
    .ret <- paste0("$OMEGA BLOCK(", .dim[1], ") ; ",
                   paste(dimnames(om0)[[1]], collapse=" "), "\n")
    .vec <- om0[lower.tri(om0,TRUE)]
    .i <- .j <- 1
    .ret <- paste0(.ret, "  ")
    for (k in .vec) {
      .ret <- paste0(.ret, " ", signif(k, .sigdig))
      if (.i == .j) {
        .i <- 1
        .j <- .j + 1
        .ret <- paste0(.ret, "\n  ")
      }
    }
    if (.fix) .ret <- paste(.ret, " FIX")
    return(paste0(.ret, "\n"))
  } else {
    paste0("$OMEGA ", signif(om0[1, 1], .sigdig), ifelse(.fix, " FIX", ""),
           " ; ", .dimn1, "\n")
  }
}

#' @export
rxUiGet.nonmemOmega <- function(x, ...) {
  .ui <- x[[1]]
  .lst <- lotri::lotriMatInv(lotri::lotriEst(lotri::as.lotri(.ui$iniDf),drop=TRUE))
  paste(vapply(seq_along(.lst),
         function(m) {
           .nonmemHandleOneOmega(.lst[[m]], .ui)
         }, character(1), USE.NAMES=FALSE), collapse="")
}
