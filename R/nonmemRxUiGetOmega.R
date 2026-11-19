.nonmemHandleOneOmega <- function(om0, ui, rec="$OMEGA", fix=NULL) {
  .dim <- dim(om0)
  .dimn1 <- dimnames(om0)[[1]][1]
  if (is.null(fix)) {
    .w <- which(ui$iniDf$name == .dimn1)
    .fix <- ui$iniDf$fix[.w]
  } else {
    .fix <- fix
  }
  .block <- ""
  .sigdig <- rxode2::rxGetControl(ui, "iniSigDig", 5)
  if (.dim[1] > 1L) {
    .ret <- paste0(rec, " BLOCK(", .dim[1], ") ; ",
                   paste(dimnames(om0)[[1]], collapse=" "), "\n")
    # NONMEM reads a block row by row down the lower triangle
    .ret <- paste0(.ret,
                   paste(vapply(seq_len(.dim[1]), function(i) {
                     paste0("  ", paste0(" ", signif(om0[i, seq_len(i)], .sigdig),
                                         collapse=""))
                   }, character(1), USE.NAMES=FALSE), collapse="\n"))
    if (.fix) .ret <- paste(.ret, " FIX")
    return(paste0(.ret, "\n"))
  } else {
    paste0(rec, " ", signif(om0[1, 1], .sigdig), ifelse(.fix, " FIX", ""),
           " ; ", .dimn1, "\n")
  }
}

#' @export
rxUiGet.nonmemOmega <- function(x, ...) {
  .ui <- x[[1]]
  .lst <- .nonmemOmegaBlocks(.ui)
  paste(vapply(seq_along(.lst),
         function(m) {
           .nonmemHandleOneOmega(.lst[[m]], .ui)
         }, character(1), USE.NAMES=FALSE), collapse="")
}
attr(rxUiGet.nonmemOmega, "rstudio") <- lotri::lotri(a~1)
