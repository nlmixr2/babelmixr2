## Declared non-Gaussian random effect distributions in the NONMEM
## translation.
##
## By the time this sees a model, nlmixr2est's pre-processing hook has run
## `rxode2::rxEtaDistExpand()`, so the declaration is already ordinary
## model code: a latent standard normal with a fixed unit `$OMEGA`, the
## copula's correlation in ordinary thetas, and a `phiU()` + inverse CDF
## line per declared random effect.  `phiU()` is NONMEM's `PHI(x)+DEL`
## idiom and `tanh()` is `DTANH()`, both added to the translation table --
## so for every family whose quantile function is ELEMENTARY the whole
## thing translates as plain NONMEM arithmetic, with nothing special to
## write.
##
## What does not translate is the three families whose quantile function
## is not elementary -- gamma, beta and Student t.  NONMEM reaches those
## through a separate mechanism entirely: `$ABBR FUNCTION GAMMACDFINV(VQ,10)`
## plus a `VQ` argument vector filled slot by slot (Bauer, NONMEM 7.5.1,
## gamma_indpar.pdf), which is a sequence of statements rather than an
## expression and so cannot come out of the expression translator.  Those
## are refused here, naming the function and what NONMEM would need,
## rather than translated into something that means a different model.

#' rxode2 inverse CDFs that NONMEM has no expression form for
#'
#' Each maps to the `$ABBR FUNCTION` routine NONMEM would need, so the
#' error can say what is actually missing.
#'
#' @noRd
#' @author Matthew L. Fidler
.nonmemEtaDistNoExpr <- c(
  "gammapInv"   = "GAMMACDFINV",
  "gammaqInv"   = "GAMMACDFINV",
  "ibetaInv"    = "a beta inverse CDF, which NONMEM does not ship",
  "studentTInv" = "a Student t inverse CDF, which NONMEM does not ship")

#' Refuse a declared distribution NONMEM cannot be given as an expression
#'
#' Called before the control stream is written, so the message names the
#' random effect the user declared rather than an internal function they
#' never wrote.
#'
#' @param ui rxode2 ui, after `rxEtaDistExpand()`
#' @return nothing, called for the error
#' @noRd
#' @author Matthew L. Fidler
.nonmemAssertEtaDist <- function(ui) {
  .info <- try(get("etaDistInfo", envir=rxode2::rxUiDecompress(ui)), silent=TRUE)
  if (inherits(.info, "try-error") || is.null(.info)) return(invisible())
  .txt <- paste(vapply(ui$lstExpr, deparse1, character(1), USE.NAMES=FALSE),
                collapse="\n")
  .bad <- names(.nonmemEtaDistNoExpr)[
    vapply(names(.nonmemEtaDistNoExpr),
           function(.f) regexpr(paste0("\\b", .f, "\\("), .txt) != -1,
           logical(1), USE.NAMES=FALSE)]
  if (length(.bad) == 0L) return(invisible())
  .nms <- .info$etaDist$name
  stop("the distribution(s) declared on '", paste(.nms, collapse="', '"),
       "' need ", paste(unique(.nonmemEtaDistNoExpr[.bad]), collapse=" and "),
       ", which NONMEM reaches through a '$ABBR FUNCTION' argument vector ",
       "rather than an expression; babelmixr2 does not write that yet.  ",
       "A family with an elementary quantile function -- see ",
       "'lotri::lotriEtaDists()' -- translates as ordinary NONMEM arithmetic",
       call.=FALSE)
}
