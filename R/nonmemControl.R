#' NONMEM estimation control
#'
#' @param est NONMEM estimation method
#' @param advanOde The ODE solving method for NONMEM
#' @param cov The NONMEM covariance method
#' @param maxeval NONMEM's maxeval (for non posthoc methods)
#' @param tol NONMEM tolerance for ODE solving advan
#' @param sigl NONMEM sigl estimation option
#' @param sigdig the significant digits for NONMEM
#' @param print The print number for NONMEM
#' @param extension NONMEM file extensions
#' @param protectZeros Add methods to protect divide by zero
#' @param noabort Add the `NOABORT` option for `$EST`
#'
#' @return babelmixr control option for generating NONMEM control stream and reading it back into
#'   `babelmixr2`/`nlmixr2`
#'
#' @author Matthew L. Fidler
#'
#' @export
#'
#' @examples
#' nonmemControl()
nonmemControl <- function(est=c("focei", "posthoc"),
                          advanOde=c("advan13", "advan8", "advan6"),
                          cov=c("r,s", "r", "s", ""),
                          maxeval=100000,
                          tol=6,
                          sigl=12,
                          sigdig=3,
                          print=1,
                          extension=getOption("babelmixr2.nmModelExtension", ".nmctl"),
                          outputExtension=getOption("babelmixr2.nmOutputExtension", ".lst"),
                          runCommand=getOption("babelmixr2.nonmem", ""),
                          iniSigDig=5,
                          protectZeros=TRUE,
                          muRef=FALSE,
                          noabort=TRUE) {
  # nonmem manual slides suggest tol=6, sigl=6 sigdig=2
  checkmate::assertIntegerish(maxeval, lower=100, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(sigdig, lower=1, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(print, lower=1, len=1, any.missing=FALSE)
  checkmate::assertLogical(noabort, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(tol, lower=1, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(sigl, lower=1, upper=14, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(iniSigDig, lower=1, len=1, any.missing=FALSE)
  checkmate::assertLogical(protectZeros, len=1, any.missing=FALSE)
  checkmate::assertLogical(muRef, len=1, any.missing=FALSE)
  if (runCommand != "") checkmate::assertCharacter(runCommand, pattern="%s", min.len=1, max.len=1)
  .ret <- list(est=match.arg(est),
               cov=match.arg(cov),
               advanOde=match.arg(advanOde),
               maxeval=maxeval,
               print=print,
               noabort=noabort,
               iniSigDig=iniSigDig,
               tol=tol,
               sigl=sigl,
               muRef=muRef,
               sigdig=sigdig,
               runCommand=runCommand,
               outputExtension=outputExtension)
  class(.ret) <- "nonmemControl"
  .ret
}

#' @export
getValidNlmixrCtl.nonmem <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- nonmemControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("nonmemControl", .ctl)
  if (!inherits(.ctl, "nonmemControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- nonmemControl()
  } else {
    .ctl <- do.call(nonmemControl, .ctl)
  }
  .ctl
}

#' @export
nmObjHandleControlObject.nonmemControl <- function(control, env) {
  assign("nonmemControl", control, envir=env)
}

#' @export
nmObjGetControl.nonmem <- function(x, ...) {
  .env <- x[[1]]
  if (exists("nonmemControl", .env)) {
    .control <- get("nonmemControl", .env)
    if (inherits(.control, "nonmemControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "nonmemControl")) return(.control)
  }
  stop("cannot find nonmem related control object", call.=FALSE)
}

