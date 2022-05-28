.nonmemThetaPad <- function(what, left=FALSE) {
  .s <- seq_along(what)
  .max <- max(vapply(.s, function(i) nchar(what[i]), integer(1), USE.NAMES=FALSE))
  vapply(.s, function(i){
    .nc <- nchar(what[i])
    if (.nc == .max) return(what[i])
    .pad <- paste(rep(" ", .max - .nc), collapse='')
    if (left) return(paste0(.pad, what[i]))
    paste0(what[i], .pad)
  }, character(1), USE.NAMES=FALSE)
}

#' @export
rxUiGet.nonmemTheta <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .theta <- .iniDf[!is.na(.iniDf$ntheta), ]
  .s <- seq_along(.theta$est)
  .sigdig <- rxode2::rxGetControl(x[[1]], "iniSigDig", 5)
  .ret <- data.frame(t0=vapply(.s, function(i) ifelse(i == 1, "$THETA", ""),
                               character(1), USE.NAMES=FALSE),
                     t1=vapply(.s,
                               function(i) {
                                 .fix <- .theta$fix[i]
                                 # (est, fixed)
                                 if (.fix) return(paste0(signif(.theta$est[i], .sigdig), ", "))
                                 # (lower, )
                                 .lower <- .theta$lower[i]
                                 if (is.finite(.lower)) return(paste0(signif(.lower, .sigdig), ", "))
                                 # (-INF, )
                                 .upper <- .theta$upper[i]
                                 if (is.finite(.upper)) return("-INF, ")
                                 # (est, )
                                 if (.theta$est[i] == 0.0) return(paste0("0.",
                                                                         paste(rep("0", .sigdig),collapse=""),
                                                                         "1"))
                                 paste(signif(.theta$est[i], .sigdig))
                               }, character(1), USE.NAMES=TRUE),
                     t2=vapply(.s,
                               function(i) {
                                 .fix <- .theta$fix[i]
                                 # (est, fixed)
                                 if (.fix) return("FIXED")
                                 if (!is.finite(.theta$lower[i]) &&
                                       !is.finite(.theta$upper[i])) return("")
                                 if (.theta$est[i] == 0.0) return(paste0("0.",
                                                                         paste(rep("0", .sigdig), collapse=""),
                                                                         "1, "))
                                 paste0(signif(.theta$est[i], .sigdig),
                                        ifelse(is.finite(.theta$upper[i]), ", ", ""))
                               }, character(1), USE.NAMES=FALSE),
                     t3=vapply(.s,
                               function(i) {
                                 .fix <- .theta$fix[i]
                                 # (est, fixed)
                                 if (.fix) return("")
                                 if (!is.finite(.theta$lower[i]) &&
                                       !is.finite(.theta$upper[i])) return("")
                                 if (is.finite(.theta$upper[i]))
                                   return(paste(signif(.theta$upper[i], .sigdig)))
                                 ""
                               }, character(1), USE.NAMES=FALSE),
                     t4=vapply(.s,
                               function(i) paste(.theta$ntheta[i]),
                               character(1), USE.NAMES=FALSE),
                     t5=vapply(.s,
                               function(i) paste(.theta$name[i]),
                               character(1), USE.NAMES=FALSE))

  .ret$t0 <- .nonmemThetaPad(.ret$t0)
  .ret$t1 <- .nonmemThetaPad(.ret$t1)
  .ret$t2 <- .nonmemThetaPad(.ret$t2)
  .ret$t3 <- .nonmemThetaPad(.ret$t3)
  .ret$t4 <- .nonmemThetaPad(.ret$t4, TRUE)
  .ret$t5 <- .nonmemThetaPad(.ret$t5)
  paste(with(.ret,
             paste0(t0, " (", t1, t2, t3, ") ; ", t4, " - ", t5)),
        collapse="\n")
}
