.rxNMcnt <- c(
  # "band"
  # "bsmm"
  # "categorical"
  # "categories"
  ########
  # "amtDose"
  # "inftDose"
  #tlast="tDose",
  time="TIME",
  M_E="2.718281828459045090796",
  M_LOG2E="1.442695040888963387005",
  M_LOG10E="0.4342944819032518166679",
  M_LN2="0.6931471805599452862268",
  M_LN10="2.302585092994045901094",
  M_PI="3.141592653589793115998",
  M_PI_2="1.570796326794896557999",
  M_PI_4="0.7853981633974482789995",
  M_1_PI="0.3183098861837906912164",
  M_2_PI="0.6366197723675813824329",
  M_2_SQRTPI="1.128379167095512558561",
  M_SQRT2="1.414213562373095145475",
  M_SQRT1_2="0.707106781186547461715",
  M_SQRT_3="1.732050807568877193177",
  M_SQRT_32="5.656854249492380581898",
  M_LOG10_2="0.3010299956639811980175",
  M_2PI="6.283185307179586231996",
  M_SQRT_PI="1.772453850905515881919",
  M_1_SQRT_2PI="0.3989422804014327028632",
  M_SQRT_2dPI="0.7978845608028654057264",
  M_LN_SQRT_PI="0.5723649429246999709164",
  M_LN_SQRT_2PI="0.918938533204672669541",
  M_LN_SQRT_PId2="0.2257913526447273278031",
  pi="3.141592653589793115998"
  )

.rxNMbad <- c("NA", "NaN", "Inf", "newind", "NEWIND")

.rxNMbadF <- c("digamma", "trigamma", "tetragamma", "pentagamma", "psigamma", "choose", "lchoose", "qnorm")

.rxNMsingle <- list(
  "gammafn" = c("DEXP(GAMLN(", "))"),
  "lgammafn" = c("GAMLN(", ")"),
  "lgamma" = c("GAMLN(", ")"),
  "loggamma" = c("GAMLN(", ")"),
  "cospi" = c("DCOS(3.141592653589793115998*(", "))"),
  "sinpi" = c("DSIN(3.141592653589793115998*(", "))"),
  "tanpi" = c("DTAN(3.141592653589793115998*(", "))"),
  "log1p" = c("DLOG(1+", ")"),
  "expm1" = c("(DEXP(", ")-1)"),
  "lfactorial" = c("GAMLN((", ")+1)"),
  "lgamma1p" = c("GAMLN((", ")+1)"),
  "expm1" = c("(DEXP(", ")-1)"),
  "log10" = c("DLOG10(", ")"),
  "log2" = c("(DLOG(", ")*1.442695040888963387005)"),
  "log1pexp" = c("DLOG(1+DEXP(", "))", "log1pexp"),
  "phi" = c("PHI(", ")"),
  "pnorm" = c("PHI(", ")"),
  "fabs"=c("DABS(", ")"),
  "sqrt"=c("DSQRT(", ")"),
  "exp"=c("DEXP(", ")"),
  "abs"=c("DABS(", ")"),
  "log"=c("DLOG(", ")"),
  "log10"=c("DLOG10(", ")"),
  "normcdf"=c("PHI(", ")"),
  "sin"=c("DSIN(", ")"),
  "cos"=c("DCOS(", ")"),
  "tan"=c("DTAN(", ")"),
  "asin"=c("DASIN(", ")"),
  "acos"=c("DACOS(", ")"),
  "atan"=c("DATAN(", ")"),
  "sinh"=c("DSINH(", ")"),
  "cosh"=c("DCOSH(", ")"),
  "atan2"=c("DATAN2(", ")"),
  "floor"=c("FLOOR(", ")"),
  "ceil"=c("CEILING(", ")"),
  "factorial" = c("DEXP(GAMLN((", ")+1))")
)

.rxNMeq <- c()
#' Handle numbers and symbols
#'
#' @param x Expression
#' @param ui User interface
#' @return Symbol, converted to NONMEM compatible name
#' @author Matthew L. Fidler
#' @noRd
.rxToNonmemHandleNamesOrAtomic <- function(x, ui=NULL) {
  if (is.character(x)) {
    stop("strings in nlmixr<->monolix are not supported", call.=FALSE)
  } else {
    if (is.na(.ret) | (.ret %in% .rxNMbad)) {
      stop("'", .ret, "' cannot be translated to NONMEM", call.=FALSE)
    }
  }
  .v <- .rxMcnt[.ret]
  if (is.na(.v)) {
    if (is.numeric(.ret)) {
      return(.ret)
    } else if (regexpr("(?:-)?(?:(?:0|(?:[1-9][0-9]*))|(?:(?:[0-9]+\\.[0-9]*)|(?:[0-9]*\\.[0-9]+))(?:(?:[Ee](?:[+\\-])?[0-9]+))?|[0-9]+[Ee](?:[\\-+])?[0-9]+)",
                       .ret, perl=TRUE) != -1) {
      return(.ret)
    } else {
      return(gsub("[.]", "__", toupper(.ret)))
    }
  } else {
    return(.v)
  }
}

.rxIsPossibleBinaryOperator <- function(expr) {
  identical(expr, quote(`*`)) ||
    identical(expr, quote(`^`)) ||
    identical(expr, quote(`+`)) ||
    identical(expr, quote(`-`)) ||
    identical(expr, quote(`/`))
}

.rxIsLogicalOperator <- function(expr) {
  identical(expr, quote(`==`)) ||
    identical(expr, quote(`>`)) ||
    identical(expr, quote(`<`)) ||
    identical(expr, quote(`<=`)) ||
    identical(expr, quote(`>=`)) ||
    identical(expr, quote(`!=`)) ||
    identical(expr, quote(`&&`)) ||
    identical(expr, quote(`||`)) ||
    identical(expr, quote(`|`)) ||
    identical(expr, quote(`&`))
}


.rxToNonmem <- function(x, ui) {
  if (is.name(x) || is.atomic(x)) {
    return(.rxToNonmemHandleNamesOrAtomic(x))
  } else if (is.call(x)) {
    if (identical(x[[1]], quote(`(`))) {
      return(paste0("(", .rxToNonmem(x[[2]], ui=ui), ")"))
    } else if (identical(x[[1]], quote(`{`))) {
      .x2 <- x[-1]
      .ret <- paste(lapply(.x2, function(x) {
        .rxToNonmem(x, ui=ui)
      }), collapse = "\n")
      return(.ret)
    } else if (.rxIsPossibleBinaryOperator(x[[1]])) {
    }
  }

}
