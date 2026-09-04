# Constants inlined as literals; PharmML has math:Constant for pi and e but
# emitting a Real keeps the walker uniform and round-trips exactly.
.rxPmlCnt <- c(
  pi        = "3.141592653589793",
  M_E       = "2.718281828459045",
  M_LN2     = "0.6931471805599453",
  M_LN10    = "2.302585092994046",
  M_PI      = "3.141592653589793",
  M_PI_2    = "1.570796326794897",
  M_PI_4    = "0.7853981633974483",
  M_SQRT2   = "1.414213562373095",
  M_LOG2E   = "1.442695040888963",
  M_LOG10E  = "0.4342944819032518"
)

# rxode2 constructs with no PharmML equivalent.
.rxPmlBad <- c("NA", "NaN", "Inf", "newind", "NEWIND")

.rxPmlBadF <- c("linCmt", "mtime", "tad", "digamma", "trigamma", "tetragamma",
                "pentagamma", "psigamma", "choose", "lchoose")

#' Translate an rxode2 expression to a PharmML math tree
#'
#' @param x A `call`, `name` or scalar, as produced by `str2lang()`.
#'
#' @param ui Optional rxode2 UI, used to resolve which block a symbol lives in
#'   so that `blkIdRef` can be emitted.  When `NULL` (the default) bare
#'   `symbIdRef` references are emitted, which is what unit tests want.
#'
#' @param indent Indent depth.
#'
#' @return character(1) holding the emitted PharmML subtree.
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.rxToPharmml <- function(x, ui = NULL, indent = 0L) {
  if (is.numeric(x) || is.integer(x)) return(.rxToPharmmlScalar(x, indent))
  if (is.name(x)) return(.rxToPharmmlSymbol(x, ui, indent))
  if (is.call(x)) return(.rxToPharmmlCall(x, ui, indent))
  stop("cannot translate '", deparse1(x), "' to PharmML", call. = FALSE)
}

#' @noRd
.rxToPharmmlScalar <- function(x, indent = 0L) {
  if (is.integer(x)) return(.pmlText("ct:Int", x, indent))
  .pmlText("ct:Real", x, indent)
}

#' @noRd
.rxToPharmmlSymbol <- function(x, ui = NULL, indent = 0L) {
  .n <- as.character(x)
  if (.n %in% .rxPmlBad) {
    stop("'", .n, "' has no PharmML equivalent", call. = FALSE)
  }
  if (.n %in% names(.rxPmlCnt)) {
    return(.pmlText("ct:Real", .rxPmlCnt[[.n]], indent))
  }
  if (.n == "time" || .n == "t") {
    return(.pmlNode("ct:SymbRef", attrs = c(symbIdRef = "t"), indent = indent))
  }
  .pmlNode("ct:SymbRef",
           attrs = .rxToPharmmlSymbAttrs(.n, ui),
           indent = indent)
}

#' Build the attribute set for a SymbRef, adding blkIdRef when resolvable
#'
#' @param n symbol name
#' @param ui rxode2 UI or NULL
#' @return named character vector
#' @noRd
.rxToPharmmlSymbAttrs <- function(n, ui = NULL) {
  .blk <- .rxToPharmmlBlockOf(n, ui)
  if (is.na(.blk)) return(c(symbIdRef = n))
  c(blkIdRef = .blk, symbIdRef = n)
}

#' Which PharmML block a symbol belongs to
#'
#' @param n symbol name
#' @param ui rxode2 UI or NULL
#' @return block id, or NA_character_ when unresolvable
#' @noRd
.rxToPharmmlBlockOf <- function(n, ui = NULL) {
  # Stage 3 populates this from the ui; with no ui there is nothing to resolve.
  NA_character_
}

# R binary operators -> math:Binop/@op
.rxPmlBinop <- c("+" = "plus", "-" = "minus", "*" = "times",
                 "/" = "divide", "^" = "power")

# Two-argument R functions -> math:Binop/@op
.rxPmlBinopF <- c(atan2 = "atan2", max = "max", min = "min")

# One-argument R functions mapping straight onto math:Uniop/@op.  PharmML's
# Uniop vocabulary is far richer than Monolix's, so most rxode2 functions land
# here rather than needing a rewrite.
.rxPmlUniop <- c(
  "-" = "minus",
  abs = "abs", exp = "exp", log = "log", log2 = "log2", log10 = "log10",
  sqrt = "sqrt", floor = "floor", ceiling = "ceiling", sign = "sign",
  factorial = "factorial", lfactorial = "factln",
  gammafn = "gamma", lgammafn = "gammaln", lgamma = "gammaln",
  loggamma = "gammaln",
  logit = "logit", expit = "logistic", probit = "probit",
  pnorm = "normcdf", phi = "normcdf",
  sin = "sin", cos = "cos", tan = "tan",
  asin = "arcsin", acos = "arccos", atan = "arctan",
  sinh = "sinh", cosh = "cosh", tanh = "tanh",
  asinh = "arcsinh", acosh = "arccosh", atanh = "arctanh"
)

# R logical binary operators -> math:LogicBinop/@op
.rxPmlLogicBinop <- c("<" = "lt", "<=" = "leq", ">" = "gt", ">=" = "geq",
                      "==" = "eq", "!=" = "neq", "&" = "and", "&&" = "and",
                      "|" = "or", "||" = "or")

# One-argument R functions with no direct Uniop; rewritten as an equivalent
# expression before translation.
.rxPmlRewrite <- list(
  log1p = function(a) bquote(log(1 + .(a))),
  expm1 = function(a) bquote(exp(.(a)) - 1),
  cospi = function(a) bquote(cos(pi * .(a))),
  sinpi = function(a) bquote(sin(pi * .(a))),
  tanpi = function(a) bquote(tan(pi * .(a))),
  log1pexp = function(a) bquote(log(1 + exp(.(a)))),
  lgamma1p = function(a) bquote(lgamma(1 + .(a)))
)

#' @noRd
.rxToPharmmlCall <- function(x, ui = NULL, indent = 0L) {
  .fn <- as.character(x[[1]])
  .nargs <- length(x) - 1L

  if (.fn %in% .rxPmlBadF) {
    stop("'", .fn, "()' has no PharmML equivalent", call. = FALSE)
  }

  # ( a ) is transparent
  if (.fn == "(") return(.rxToPharmml(x[[2]], ui, indent))

  if (.fn == "!" && .nargs == 1L) {
    return(.pmlNode("math:LogicUniop", attrs = c(op = "not"),
                    children = .rxToPharmml(x[[2]], ui),
                    indent = indent))
  }

  if (.fn %in% names(.rxPmlLogicBinop) && .nargs == 2L) {
    return(.pmlNode("math:LogicBinop",
                    attrs = c(op = .rxPmlLogicBinop[[.fn]]),
                    children = c(.rxToPharmml(x[[2]], ui),
                                 .rxToPharmml(x[[3]], ui)),
                    indent = indent))
  }

  if (.fn == "ifelse" && .nargs == 3L) {
    return(.rxToPharmmlPiecewise(x[[2]], x[[3]], x[[4]], ui, indent))
  }

  if (.fn == "if") {
    return(.rxToPharmmlPiecewise(x[[2]], x[[3]],
                                 if (length(x) > 3L) x[[4]] else NULL,
                                 ui, indent))
  }

  # unary minus/plus
  if (.fn %in% c("-", "+") && .nargs == 1L) {
    if (.fn == "+") return(.rxToPharmml(x[[2]], ui, indent))
    return(.pmlNode("math:Uniop", attrs = c(op = "minus"),
                    children = .rxToPharmml(x[[2]], ui),
                    indent = indent))
  }

  if (.fn %in% names(.rxPmlBinop) && .nargs == 2L) {
    return(.pmlNode("math:Binop",
                    attrs = c(op = .rxPmlBinop[[.fn]]),
                    children = c(.rxToPharmml(x[[2]], ui),
                                 .rxToPharmml(x[[3]], ui)),
                    indent = indent))
  }

  if (.fn %in% names(.rxPmlBinopF) && .nargs == 2L) {
    return(.pmlNode("math:Binop",
                    attrs = c(op = .rxPmlBinopF[[.fn]]),
                    children = c(.rxToPharmml(x[[2]], ui),
                                 .rxToPharmml(x[[3]], ui)),
                    indent = indent))
  }

  if (.fn %in% names(.rxPmlRewrite) && .nargs == 1L) {
    return(.rxToPharmml(.rxPmlRewrite[[.fn]](x[[2]]), ui, indent))
  }

  if (.fn %in% names(.rxPmlUniop) && .nargs == 1L) {
    return(.pmlNode("math:Uniop",
                    attrs = c(op = .rxPmlUniop[[.fn]]),
                    children = .rxToPharmml(x[[2]], ui),
                    indent = indent))
  }

  stop("cannot translate '", deparse1(x), "' to PharmML", call. = FALSE)
}

#' Emit a math:Piece
#'
#' @param value already-emitted value expression
#' @param cond already-emitted condition expression, or NULL for Otherwise
#' @return character(1)
#' @noRd
.rxToPharmmlPiece <- function(value, cond = NULL) {
  .condChild <- if (is.null(cond)) .pmlNode("math:Otherwise") else cond
  .pmlNode("math:Piece",
           children = c(value, .pmlNode("math:Condition", children = .condChild)))
}

#' Translate an if/else or ifelse() into a math:Piecewise
#'
#' @param test unevaluated condition
#' @param yes unevaluated true branch
#' @param no unevaluated false branch, or NULL when absent
#' @param ui rxode2 UI or NULL
#' @param indent indent depth
#' @return character(1)
#' @noRd
.rxToPharmmlPiecewise <- function(test, yes, no = NULL, ui = NULL, indent = 0L) {
  .pieces <- .rxToPharmmlPiece(.rxToPharmml(yes, ui), .rxToPharmml(test, ui))
  if (!is.null(no)) {
    .pieces <- c(.pieces, .rxToPharmmlPiece(.rxToPharmml(no, ui)))
  }
  .pmlNode("math:Piecewise", children = .pieces, indent = indent)
}
