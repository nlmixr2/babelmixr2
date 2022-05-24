nonmemControl <- function(est=c("focei", "posthoc"),
                          cov=c("r,s", "r", "s", ""),
                          maxeval=100000,
                          sigdig=3,
                          print=1,
                          noabort=TRUE) {
  checkmate::assertIntegerish(maxeval, lower=100, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(sigdig, lower=1, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(print, lower=1, len=1, any.missing=FALSE)
  checkmate::assertLogical(abort, len=1, any.missing=FALSE)
  .ret <- list(est=match.arg(est),
               cov=match.arg(cov),
               maxeval=maxeval,
               print=print,
               noabort=noabort)
  class(.ret) <- "nonmemControl"
  .ret
}
