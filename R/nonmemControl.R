nonmemControl <- function(est=c("focei", "posthoc"),
                          advanOde=c("advan13", "advan8", "advan6"),
                          cov=c("r,s", "r", "s", ""),
                          maxeval=100000,
                          tol=6,
                          sigl=6,
                          sigdig=2,
                          print=1,
                          noabort=TRUE) {
  # nonmem manual slides suggest tol=6, sigl=6 sigdig=2
  checkmate::assertIntegerish(maxeval, lower=100, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(sigdig, lower=1, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(print, lower=1, len=1, any.missing=FALSE)
  checkmate::assertLogical(abort, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(tol, lower=1, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(sigl, lower=1, upper=14, len=1, any.missing=FALSE)
  .ret <- list(est=match.arg(est),
               cov=match.arg(cov),
               advanOde=match.arg(advanOde),
               maxeval=maxeval,
               print=print,
               noabort=noabort,
               tol=tol,
               sigl=sigl,
               sigdig=sigdig)
  class(.ret) <- "nonmemControl"
  .ret
}
