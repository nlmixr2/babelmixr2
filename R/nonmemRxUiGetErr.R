.getErr <- function(err, indent=TRUE) {
  .err <- readLines(file.path(system.file(package="babelmixr2"), err))
  .ret <- paste0(paste(gsub("^      ",
                            ifelse(indent, "    ", "  "), .err),
                       collapse="\n"))
  .ret <- strsplit(.ret, "\n+")[[1]]
  paste(.ret[regexpr("^ *$", .ret) == -1], collapse="\n")
}

.nonmemErr0 <- function(ui, pred1, indent=TRUE) {
  if (length(pred1$cond) == 1L) {
    .iniDf <- ui$iniDf
    if (pred1$transform == "untransformed") return(.getErr("err.txt",  indent=indent))
    .transform <- paste(pred1$transform)
    .low <- paste(pred1$trLow)
    if (regexpr("^[0-9]+$", .low) != -1) .low <- paste0(.low, ".0")
    .hi <- paste(pred1$trHi)
    if (regexpr("^[0-9]+$", .hi) != -1) .hi <- paste0(.hi, ".0")
    if (.transform == "boxCox") {
      .w <- which(.iniDf$condition == pred1$cond & .iniDf$err == "boxCox")
      .t <- (.nonmemGetThetaNum(.iniDf$name[.w],ui))
      return(gsub("LAMBDA", .t, .getErr("err-boxCox.txt",  indent=indent)))
    } else if (.transform == "yeoJohnson") {
      # The rest are not
      .w <- which(.iniDf$condition == pred1$cond &.iniDf$err == "yeoJohnson")
      .t <- (.nonmemGetThetaNum(.iniDf$name[.w],ui))
      return(gsub("LAMBDA", .t, .getErr("err-yeoJohnson.txt",  indent=indent)))
    } else if (.transform == "lnorm") {
      return(.getErr("err-lnorm.txt", indent=indent))
    } else if (.transform == "logit") {
      return(gsub("LOW", .low,
                  gsub("HIGH", .hi,
                       .getErr("err-logit.txt",  indent=indent))))
    } else if (.transform == "logit + yeoJohnson") {
      .w <- which(.iniDf$condition == pred1$cond &.iniDf$err == "yeoJohnson")
      .t <- (.nonmemGetThetaNum(.iniDf$name[.w],ui))
      return(gsub("LAMBDA", .t,
                  gsub("LOW", .low,
                       gsub("HIGH", .hi,
                            .getErr("err-logitYeoJohnson.txt",  indent=indent)))))
    }
  }
  return(.getErr("err.txt"))
}

.repEndpoint <- function(var, dvid) {
  .pf <- paste0("RX_PF", dvid)
  .p <- paste0("RX_P", dvid)
  .i <- paste0("RX_IP", dvid)
  .w <- paste0("W", dvid)
  .ret <- gsub("RX_PRED_F_", .pf, var)
  .ret <- gsub("RX_PRED_", .p, .ret)
  .ret <- gsub("IPRED", .i, .ret)
  .ret
}

#'@export
rxUiGet.nonmemErrF <- function(x, ...) {
  .ui <- x[[1]]
  rxode2::rxAssignControlValue(.ui, ".ifelse", TRUE)
  on.exit(rxode2::rxAssignControlValue(.ui, ".ifelse", FALSE))
  .predDf <- .ui$predDf
  .cmtCnt <- rxode2::rxGetControl(.ui, ".cmtCnt", rep(0L, length(.predDf$cond)))
  .ipred <- vapply(seq_along(.predDf$cond),
                   function(i){
                     .pred1 <- .predDf[i, ]
                     .ret <- .nonmemErr0(.ui, .pred1, indent=FALSE)
                     .ret <- .repEndpoint(.ret, .pred1$dvid)
                     .w <- str2lang(paste0("W", .pred1$dvid))
                     .var <- paste0("\n  RX_PRED_ = IPRED\n",
                                    .rxToNonmem(bquote(.(.w) ~
                                                         sqrt(.(rxode2::.rxGetVarianceForErrorType(.ui, .pred1)))),
                                                .ui))
                     .var <- .repEndpoint(.var, .pred1$dvid)
                     # depending on the method the prop can be with regards to the F or the transformed F
                     # So, here we add RX_PRED_ to be the transformed to support both
                     .ret <- paste0(paste0("  ; endpoint nobs=", .cmtCnt[i], "\n"),
                                    .ret,
                                    .var,
                                    paste0("\n  IF (", .w, " .EQ. 0.0) ", .w, " = 1"))
                     .ret
                   }, character(1), USE.NAMES=FALSE)
  .err <- paste(.ipred, collapse="\n")
  .cens <- rxode2::rxGetControl(.ui, ".hasCens", FALSE)
  .limit <- rxode2::rxGetControl(.ui, ".hasLimit", FALSE)
  if (length(.predDf$cond) == 1L) {
    .ipred <- "  IPRED = RX_IP2\n  W     = W1\n"
  } else {
    .ipred <- vapply(seq_along(.predDf$cond), function(i) {
      paste0("  IF (DVID .EQ. ", i, ") THEN\n    IPRED = RX_IP", i, "\n    W     = W", i, "\n  END IF\n")
    }, character(1), USE.NAMES=FALSE)
  }
  if (.cens && .limit) {
    .y <- .getErr("err-cens-limit.txt", FALSE)
  } else if (.cens) {
    .y <- .getErr("err-cens.txt", FALSE)
  } else if (.limit) {
    stop("ylo/yup not implemented (would require laplacian); drop LIMIT or add CENS column",
         call.=FALSE)
  } else {
    .y <- "  Y     = IPRED + W*EPS(1)"
  }
  paste0("\n  ; Write out expressions for ipred and w\n",
         gsub("\n *\n+",
              "\n",
              paste(c(.err,
                      .ipred,
                      .y), collapse="\n")), "\n")
}
