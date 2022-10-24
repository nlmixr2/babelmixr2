rex::register_shortcuts("babelmixr2")

.rxNMcnt <- c(
  # "band"
  # "bsmm"
  # "categorical"
  # "categories"

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

.rxNMprotectZero <-
  c("gammafn", "lgammafn", "lgamma", "loggamma", "log10", "log2", "sqrt", "log")

.rxNmProtectZeroP1 <- c("log1p", "lfactorial", "lgamma1p", "factorial")

# "log1pexp" = c("DLOG(1+DEXP(", "))", "log1pexp"), ???

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

.rxNMlogic <- c(
  `==`=".EQ.",
  `>`=".GT.",
  `>=`=".GE.",
  `<`=".LT.",
  `<=`=".LE.",
  `!=`=".NE.",
  `&&`=".AND.",
  `||`=".OR.",
  `&`=".AND.",
  `|`=".OR.")

.rxNMbin <- c(`*`="*",
              `^`="**",
              `**`="**",
              `+`="+",
              `-`="-",
              `/`="/")

#' Help show where code came from (but stop at 500 characters of deparsing and
#' stop on the first line)
#' 
#' @param x rxode2 expression line
#' @inheritParams deparse
#' @return R expression as NONMEM comement
#' @author Bill Denney
#' @noRd
.babelmixr2Deparse <- function(x, width.cutoff=500L) {
  ret <- deparse(x, width.cutoff=width.cutoff)
  if (length(ret) > 1) {
    ret <- paste0(ret[1], "...")
  }
  paste0(" ; ", ret)
}

.nmNumReg <- function(var="ETA", C=FALSE) {
  .nmEtaNum <- rex::rex(var, "1":"9")
  .nmEtaNum2 <- rex::rex(var, "1":"9", "0":"9")
  if (C) {
    .var2 <- paste0(var, "C")
    rex::rex(or(var, .var2, .nmEtaNum, .nmEtaNum2))
  } else {
    rex::rex(or(var, .nmEtaNum, .nmEtaNum2))
  }
}


.nmRes <- rex::rex(start,
                   or("GETETA", "SIMETA", "SIMEPS",
                      "COMSAV", "NWIND", "ETEXT", "IERPRD", "MSEC",
                      "MFIRST", "NETEXT", .nmNumReg("ETA"),.nmNumReg("THETA"), .nmNumReg("EPS"),
                      .nmNumReg("MU_"),.nmNumReg("UM_"),
                      .nmNumReg("A"), .nmNumReg("B"), .nmNumReg("C"), .nmNumReg("D"), .nmNumReg("E"),
                      .nmNumReg("F"), .nmNumReg("P"), .nmNumReg("Q"), .nmNumReg("MC"), .nmNumReg("ME"),
                      .nmNumReg("MG"), .nmNumReg("MT"), .nmNumReg("ROCM"),
                      .nmNumReg("S", TRUE), .nmNumReg("F", TRUE),
                      "FO", .nmNumReg("R", TRUE), .nmNumReg("D", TRUE),
                      .nmNumReg("ALAG", TRUE), .nmNumReg("TSCALE", TRUE),
                      .nmNumReg("XSCALE", TRUE), "A_0FLG", "A_0", "DADT",
                      "CALLFL", "Y", "NEWL2", "ICALL", "EXIT", "CALL",
                      "GETETA", "SIMETA", "SIMEPS", "COMSAV", "NWIND", "ETEXT",
                      "IERPRD", "MSEC", "MFIRST", "NETEXT", "IPRED",
                      "IPRE", "IPR"),
                   end)
#' Gets variable, respecting the many reserved names in NONMEM
#'
#'
#' @param var Variable name in rxode2 syntax
#' @param ui UI for saving and retrieving information
#' @return NONMEM-compatible variable name
#' @author Matthew L. Fidler
#' @noRd
.nmGetVar <- function(var, ui) {
  .reserved <- rxode2::rxGetControl(ui, ".nmGetVarReservedDf",
                                    data.frame(var=character(0),
                                               nm=character(0)))
  .uvar <- gsub(".", "_", toupper(var), fixed=TRUE)
  .w <- which(.reserved$var == var)
  if (length(.w) == 1) {
    var <- .reserved$nm[.w]
  } else if (regexpr(.nmRes, .uvar, perl=TRUE) != -1) {
    .num <- rxode2::rxGetControl(ui, ".nmVarResNum", 1)
    .newVar <- sprintf("RXR%d", .num)
    rxode2::rxAssignControlValue(ui, ".nmVarResNum", .num + 1)
    .reserved <- rbind(.reserved, data.frame(var=var, nm=.newVar))
    rxode2::rxAssignControlValue(ui, ".nmGetVarReservedDf", .reserved)
    var <- .newVar
  }
  .uvar <- gsub(".", "_", toupper(var), fixed=TRUE)
  .var <- rxode2::rxGetControl(ui, ".nmGetVarDf",
                               data.frame(var=character(0),
                                          nm=character(0)))
  .w <- which(.var$var == var)
  if (length(.w) == 1) return(.var$nm[.w])
  .w <- which(.var$nm == .uvar)
  .doRx <- FALSE
  if (length(.w) == 1) {
    .doRx <- TRUE
  }
  .extra <- rxode2::rxGetControl(ui, ".nmVarExtra", "")
  if (.doRx) {
    .num <- rxode2::rxGetControl(ui, ".nmVarNum", 1)
    .newVar <- sprintf("RX%s%03d", .extra, .num)
    rxode2::rxAssignControlValue(ui, ".nmVarNum", .num + 1)
  } else {
    .newVar <- .uvar
  }
  .var <- rbind(.var, data.frame(var=var, nm=.newVar))
  rxode2::rxAssignControlValue(ui, ".nmGetVarDf", .var)
  .newVar
}
#'
#' @param x Expression
#' @param ui User interface
#' @return Symbol, converted to NONMEM compatible name
#' @author Matthew L. Fidler
#' @noRd
.rxToNonmemHandleNamesOrAtomic <- function(x, ui) {
  if (is.character(x)) stop("strings in nlmixr<->monolix are not supported", call.=FALSE)
  .ret <- as.character(x)
  if (tolower(.ret) %in% c("t", "time")) return("TIME")
  if (exists(".thetaMu", ui)) {
    .thetaMu <- ui$.thetaMu
    .w <- which(names(.thetaMu) == .ret)
    if (length(.w) == 1) {
      .ret <- .thetaMu[.w]
      if (!is.na(.ret)) return(.ret)
    }
  }
  .ref <- .nonmemGetThetaNum(x, ui)
  if (!is.na(.ref)) return(.ref)
  .ref <- .nonmemGetEtaNum(x, ui)
  if (!is.na(.ref)) return(.ref)
  if (is.na(.ret) | (.ret %in% .rxNMbad)) {
    stop("'", .ret, "' cannot be translated to NONMEM", call.=FALSE)
  }
  .v <- .rxNMcnt[.ret]
  if (is.na(.v)) {
    if (is.numeric(.ret)) {
      .ret <- gsub("e", "D", as.character(.ret))
      return(.ret)
    } else if (regexpr("^(?:-)?(?:(?:0|(?:[1-9][0-9]*))|(?:(?:[0-9]+\\.[0-9]*)|(?:[0-9]*\\.[0-9]+))(?:(?:[Ee](?:[+\\-])?[0-9]+))?|[0-9]+[Ee](?:[\\-+])?[0-9]+)$",
                       .ret, perl=TRUE) != -1) {
      .ret <- gsub("e", "D", .ret)
      return(.ret)
    } else {
      .cmt <-.rxGetCmtNumber(.ret, ui, error=FALSE)
      if (!is.na(.cmt)) {
        return(paste0("A(", .cmt, ")"))
      }
      return(.nmGetVar(.ret, ui))
    }
  } else {
    return(.v)
  }
}

.rxIsPossibleBinaryOperator <- function(expr) {
  identical(expr, quote(`*`)) ||
    identical(expr, quote(`^`)) ||
    identical(expr, quote(`**`)) ||
    identical(expr, quote(`+`)) ||
    identical(expr, quote(`-`)) ||
    identical(expr, quote(`/`))
}

#'  Handle d/dt() expression
#'
#' @param x d/dt() expression
#' @param ui rxode2 ui
#' @return DADT(#) for NONMEM
#' @author Matthew L. Fidler
#' @noRd
.rxToNonmemHandleDdt <- function(expr, ui) {
  stopifnot(.rxIsDdt(expr))
  .cmt <-.rxGetCmtNumber(expr[[3]][[2]], ui)
  sprintf("DADT(%g)", .cmt)
}

#' Handle d/dt() line and add rxode code as a comment
#'
#' @param x expression 
#' @param ui rxode2 ui object
#' @return d/dt() nonmem line
#' @author Matthew L. Fidler with influence from Bill Denney
#' @noRd
.rxToNonmemHandleDdtLine <- function(x, ui) {
  paste0(.rxToNonmemHandleDdt(x[[2]], ui), " = ",
         .rxToNonmem(x[[3]], ui=ui))
}

#' Protect Zeros for dlog(x) or dsqrt(x)
#'
#' @param x Expression to protect
#' @param ui rxode2 to get information
#' @param one if this is protecting a plus one expression like `lfactorial()`
#' @return expression, with prefix lines calculated
#' @author Matthew L. Fidler
#' @noRd
.rxProtectPlusZero <- function(x, ui, one=FALSE) {
  .ret <- .rxToNonmem(x, ui=ui)
  if (.rxShouldProtectZeros(.ret, ui)) {
    .df <- rxode2::rxGetControl(ui, ".nmGetDivideZeroDf",
                                data.frame(expr=character(0),
                                           nm=character(0)))
    .expr <- paste0(.ret, ifelse(one, "+++1", ""))
    .w <- which(.df$expr == .expr)
    if (length(.w) == 1) {
      # Previously protected this expression
      .ret <- .df$nm[.w]
    } else {
      .prefixLines <- rxode2::rxGetControl(ui, ".nmPrefixLines", NULL)
      .num <- rxode2::rxGetControl(ui, ".nmVarDZNum", 1)
      .extra <- rxode2::rxGetControl(ui, ".nmVarExtra", "")
      .newVar <- sprintf("RXDZ%s%03d", .extra, .num)
      rxode2::rxAssignControlValue(ui, ".nmVarDZNum", .num + 1)
      .sigdig <- rxode2::rxGetControl(ui, "iniSigDig", 5)
      .num <- paste0(ifelse(one, "-1.", "0."), paste(rep("0", .sigdig), collapse=""), "1")
      .prefixLines <- c(.prefixLines,
                        paste0(.rxToNonmemGetIndent(ui),
                               .newVar, "=", .ret),
                        paste0(.rxToNonmemGetIndent(ui),
                               "IF (", .newVar, " .LE. ", .num, ") THEN"))
      .rxToNonmemIndent(ui)
      .prefixLines <- c(.prefixLines,
                        paste0(.rxToNonmemGetIndent(ui),
                               .newVar, "=", .num),
                        paste0(.rxToNonmemGetIndent(ui, FALSE),
                               "END IF"))
      .df <- rbind(.df,
                   data.frame(expr=.expr, nm=.newVar))
      rxode2::rxAssignControlValue(ui, ".nmGetDivideZeroDf", .df)
      rxode2::rxAssignControlValue(ui, ".nmPrefixLines", .prefixLines)
      .ret <- .newVar
    }
  }
  .ret
}

#'  Protect Zero expressions but preserves negative/positive
#'
#' @param x expression to protect
#' @param ui User interface
#' @return expression, but adds prefix lines to protect the expression
#' @author Matthew L. Fidler
#' @noRd
.rxProtectPlusOrMinusZero <- function(x, ui) {
  .denom <- .rxToNonmem(x, ui=ui)
  if (.rxShouldProtectZeros(.denom, ui)) {
    .df <- rxode2::rxGetControl(ui, ".nmGetDivideZeroDf",
                                data.frame(expr=character(0),
                                           nm=character(0)))
    .w <- which(.df$expr == .denom)
    if (length(.w) == 1) {
      # Previously protected this expression
      .denom <- .df$nm[.w]
    } else {
      .prefixLines <- rxode2::rxGetControl(ui, ".nmPrefixLines", NULL)
      .num <- rxode2::rxGetControl(ui, ".nmVarDZNum", 1)
      .extra <- rxode2::rxGetControl(ui, ".nmVarExtra", "")
      .newVar <- sprintf("RXDZ%s%03d", .extra, .num)
      rxode2::rxAssignControlValue(ui, ".nmVarDZNum", .num + 1)
      .sigdig <- rxode2::rxGetControl(ui, "iniSigDig", 5)
      .num <- paste0("0.", paste(rep("0", .sigdig), collapse=""), "1")
      .prefixLines <- c(.prefixLines,
                        paste0(.rxToNonmemGetIndent(ui),
                               .newVar, "=", .denom),
                        paste0(.rxToNonmemGetIndent(ui),
                               "IF (", .newVar, " .GE. 0.0 .AND. ",
                               .newVar, " .LE. ", .num, ") THEN"))
      .rxToNonmemIndent(ui)
      .prefixLines <- c(.prefixLines,
                        paste0(.rxToNonmemGetIndent(ui),
                               .newVar, "=", .num),
                        paste0(.rxToNonmemGetIndent(ui, FALSE), "END IF\n",
                               .rxToNonmemGetIndent(ui), "IF (",
                               .newVar, " .GE. -", .num, " .AND. ",
                               .newVar, " .LT. 0.) THEN"))
      .rxToNonmemIndent(ui)
      .prefixLines <- c(.prefixLines,
                        paste0(.rxToNonmemGetIndent(ui),
                               .newVar, "= -", .num),
                        paste0(.rxToNonmemGetIndent(ui, FALSE),
                               "END IF"))
      .df <- rbind(.df,
                   data.frame(expr=.denom, nm=.newVar))
      rxode2::rxAssignControlValue(ui, ".nmGetDivideZeroDf", .df)
      rxode2::rxAssignControlValue(ui, ".nmPrefixLines", .prefixLines)
      .denom <- .newVar
    }
  }
  .denom
}

#' This handles divide by zero for NONMEM control streams
#'
#' @param x2 numerator
#' @param x3 denominator
#' @param ui rxode2 ui
#' @return divion operator, as a side effect lines are prepended to protect divide by zero errors
#' @author Matthew L. Fidler
#' @noRd
.rxToNonmemHandleDivideZero <- function(x2, x3, x, ui) {
  paste0(
    .rxToNonmem(x2, ui=ui),
    .rxNMbin[as.character(x[[1]])],
    .rxProtectPlusOrMinusZero(x3, ui))
}
#' This is where binary operators are converted to NONMEM operators
#'
#' @param x Binary operator R expresion
#' @param ui rxode2 user interface function
#' @return NONMEM equivalent binary operator
#' @author Matthew L. Fidler
#' @noRd
.rxToNonmemHandleBinaryOperator <- function(x, ui) {
  if (identical(x[[1]], quote(`/`))) {
    .x2 <- x[[2]]
    .x3 <- x[[3]]
    if (.rxIsDdt(x)) {
      return(.rxToNonmemHandleDdt(x, ui))
    } else {
      if (length(.x2) == 2 && length(.x3) == 2) {
        if (identical(.x2[[1]], quote(`df`)) &&
              identical(.x3[[1]], quote(`dy`))) {
          stop('df()/dy() is not supported in NONMEM conversion', call.=FALSE)
        }
      }
      return(.rxToNonmemHandleDivideZero(.x2, .x3, x, ui))
    }
  } else if (identical(x[[1]], quote(`^`)) ||
               identical(x[[1]], quote(`**`))) {
    .needProtect <-TRUE
    if (is.numeric(x[[3]]) && x[[3]] > 0) {
      .needProtect <- FALSE
    }
    .ret <- paste0(
      ifelse(.needProtect,
             .rxProtectPlusOrMinusZero(x[[2]], ui),
             .rxToNonmem(x[[2]], ui)),
      .rxNMbin[as.character(x[[1]])],
      .rxToNonmem(x[[3]], ui=ui)
    )
  } else {
    .ret <- paste0(
      .rxToNonmem(x[[2]], ui=ui),
      .rxNMbin[as.character(x[[1]])],
      .rxToNonmem(x[[3]], ui=ui)
    )
  }
  return(.ret)
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

# Indent more
.rxToNonmemIndent <- function(ui) {
  rxode2::rxAssignControlValue(ui, ".nmIndent",
                               rxode2::rxGetControl(ui, ".nmIndent", 2) + 2)
}

# Indent less
.rxToNonmemUnIndent <- function(ui) {
  rxode2::rxAssignControlValue(ui, ".nmIndent",
                               max(2, rxode2::rxGetControl(ui, ".nmIndent", 2) - 2))
}

# Get spaces for NONMEM indentation (and maybe indent more or less)
.rxToNonmemGetIndent <- function(ui, ind=NA) {
  if (is.na(ind)) {
  } else if (ind) {
    .rxToNonmemIndent(ui)
  } else {
    .rxToNonmemUnIndent(ui)
  }
  .nindent <- rxode2::rxGetControl(ui, ".nmIndent", 2)
  paste(vapply(seq(1, .nindent), function(x) " ", character(1), USE.NAMES=FALSE), collapse="")
}

.rxToNonmemHandleIfExpressions <- function(x, ui) {
  if (rxode2::rxGetControl(ui, ".ifelse", FALSE)) {
    stop("babelmixr2 NONMEM translator will not handle nested if/else models")
  }
    #rxode2::rxAssignControlValue(ui, ".ifelse", TRUE)
  .ret <- paste0(.rxToNonmemGetIndent(ui), "IF (", .rxToNonmem(x[[2]], ui=ui), ") THEN\n")
  .rxToNonmemIndent(ui)
  rxode2::rxAssignControlValue(ui, ".ifelse", TRUE)
  on.exit(rxode2::rxAssignControlValue(ui, ".ifelse", FALSE))
  .ret <- paste0(.ret, .rxToNonmem(x[[3]], ui=ui))
  x <- x[-c(1:3)]
  if (length(x) == 1) x <- x[[1]]
  while(identical(x[[1]], quote(`if`))) {
    stop("babelmixr2 will not allow `else if` or `else` statements in NONMEM models",
         call.=FALSE)
    ## .ret <- paste0(.ret, "\n",
    ##                .rxToNonmemGetIndent(ui, FALSE), "ELSE IF (", .rxToNonmem(x[[2]], ui=ui), ") THEN\n")
    ## .rxToNonmemIndent(ui)
    ## .ret <- paste0(.ret, .rxToNonmem(x[[3]], ui=ui))
    ## x <- x[-c(1:3)]
    ## if (length(x) == 1) x <- x[[1]]
  }
  if (is.null(x)) {
    .ret <- paste0(.ret, "\n",
                   .rxToNonmemGetIndent(ui, FALSE), "END IF\n")
  }  else {
    stop("babelmixr2 will not allow `else if` or `else` statements in NONMEM models",
         call.=FALSE)
    ## .ret <- paste0(.ret, "\n",
    ##                .rxToNonmemGetIndent(ui, FALSE), "ELSE\n")
    ## .rxToNonmemIndent(ui)
    ## .ret <- paste0(.ret, .rxToNonmem(x, ui=ui),
    ##                "\n",
    ##                .rxToNonmemGetIndent(ui, FALSE), "END IF\n")
  }
  return(.ret)
}

.nonmemGetCmtProperties <- function(ui) {
  rxode2::rxGetControl(ui, ".cmtProperties",
                       data.frame(cmt=integer(0),
                                  f=character(0),
                                  dur=character(0),
                                  lag=character(0),
                                  rate=character(0),
                                  init=character(0)))
}

.nonmemSetCmtProperty <- function(ui, state, extra, type="f") {
  .prop <- .nonmemGetCmtProperties(ui)
  .state <- rxode2::rxState(ui)
  .cmt <- which(state == .state)
  .w <- which(.prop$cmt == .cmt)
  if (length(.w) == 0L) {
    .prop <- rbind(.prop,
                   data.frame(cmt=.cmt,
                              f=NA_character_,
                              dur=NA_character_,
                              lag=NA_character_,
                              rate=NA_character_,
                              init=NA_character_))
    .w <- which(.prop$cmt == .cmt)
  }
  if (type == "f") {
    .prop[.w, "f"] <- extra
  } else if (type == "dur") {
    .prop[.w, "dur"] <- extra
  } else if (type == "lag") {
    .prop[.w, "lag"] <- extra
  } else if (type == "rate") {
    .prop[.w, "rate"] <- extra
  } else if (type == "init") {
    .prop[.w, "init"] <- extra
  }
  rxode2::rxAssignControlValue(ui, ".cmtProperties", .prop)
}

.nonmemSetCmtProperty <- function(ui, state, extra, type="f") {
  .prop <- .nonmemGetCmtProperties(ui)
  .state <- rxode2::rxState(ui)
  .cmt <- .rxGetCmtNumber(state, ui)
  .w <- which(.prop$cmt == .cmt)
  if (length(.w) == 0L) {
    .prop <- rbind(.prop,
                   data.frame(cmt=.cmt,
                              f=NA_character_,
                              dur=NA_character_,
                              lag=NA_character_,
                              rate=NA_character_,
                              init=NA_character_))
    .w <- which(.prop$cmt == .cmt)
  }
  if (type == "f") {
    .prop[.w, "f"] <- extra
  } else if (type == "dur") {
    .prop[.w, "dur"] <- extra
  } else if (type == "lag") {
    .prop[.w, "lag"] <- extra
  } else if (type == "rate") {
    .prop[.w, "rate"] <- extra
  } else if (type == "init") {
    .prop[.w, "init"] <- extra
  }
  rxode2::rxAssignControlValue(ui, ".cmtProperties", .prop)
}

.rxToNonmemHandleAssignmentOperator <- function(x, ui) {
  if (length(x[[2]]) == 1L) {
    .rxToNonmemHandleAssignmentOperatorSimpleLHS(x, ui)
  } else {
    .rxToNonmemHandleAssignmentOperatorComplexLHS(x, ui)
  }
}

# When there is a simple left-hand-side assignment (e.g. set a variable)
.rxToNonmemHandleAssignmentOperatorSimpleLHS <- function(x, ui) {
  stopifnot(length(x[[2]]) == 1)
  .var <- .rxToNonmem(x[[2]], ui=ui)
  .val <- .rxToNonmem(x[[3]], ui=ui)
  .prefixLines <- rxode2::rxGetControl(ui, ".nmPrefixLines", NULL)
  .extra <- ""
  if (!is.null(.prefixLines)) {
    .extra <- paste0(paste(.prefixLines, collapse="\n"), "\n")
    rxode2::rxAssignControlValue(ui, ".nmPrefixLines", NULL)
  }
  paste0(
    .extra,
    .rxToNonmemGetIndent(ui), .var, "=", .val,
    .babelmixr2Deparse(x) # Show the source where this code came from
  )
}

# When there is a more complex left-hand-side assignment (e.g. initial condition
# setting or lag time setting)
.rxToNonmemHandleAssignmentOperatorComplexLHS <- function(x, ui) {
  if (.rxIsDdt(x[[2]])) {
    # Currently d/dt() the only 3-long option that is used in practice
    # Specifying the jacobian is possible but an error will the thrown
    # anyway in this implementation and I don't think people use it in
    # nlmixr models.  The information may be thrown away depending on
    # what procedure is run.
    .ret <- .rxToNonmemHandleDdtLine(x, ui)
  } else if (identical(x[[2]][[2]], 0)) {
    # set initial conditions
    return(paste0(.rxToNonmemGetIndent(ui),
                  .rxToNonmemHandleInitialConditions(x, ui)))
  } else if (length(x[[2]]) == 2) {
    return(paste0(.rxToNonmemGetIndent(ui),
                  .rxToNonmemHandleAssignmentPrefix(x, ui)))
  } else {
    stop("the left hand expression '", deparse1(x), "' is not supported",
         call.=FALSE)
  }
  .prefixLines <- rxode2::rxGetControl(ui, ".nmPrefixLines", NULL)
  .extra <- ""
  if (!is.null(.prefixLines)) {
    .extra <- paste0(paste(.prefixLines, collapse="\n"), "\n")
    rxode2::rxAssignControlValue(ui, ".nmPrefixLines", NULL)
  }
  paste0(
    .extra,
    .rxToNonmemGetIndent(ui),
    .ret,
    .babelmixr2Deparse(x)
  )
}

#' Handle compartment number initial conditions
#'
#' @param x rxode2 expression line 
#' @param ui User interface
#' @return String for NONMEM style initial condition
#' @author Matt Fidler and Bill Denney
#' @noRd
.rxToNonmemHandleInitialConditions <- function(x, ui) {
  .state <-  as.character(x[[2]][[1]])
  # Cannot use ifelse in the block
  rxode2::rxAssignControlValue(ui, ".ifelse", TRUE)
  on.exit(rxode2::rxAssignControlValue(ui, ".ifelse", FALSE))
  if (length(x[[3]]) != 1L) {
    stop("the complex initial condition is not supported in the nonmem conversion\n",
         deparse1(x), call.=FALSE)
  }
  .extra <- paste0(.rxToNonmem(x[[3]], ui=ui),
                   .babelmixr2Deparse(x))
  .nonmemSetCmtProperty(ui, .state, .extra, type="init")
  return(paste0(.rxToNonmemGetIndent(ui),
                    ";", .state, "(0) defined in $PK block"))
}

.rxToNonmemHandleAssignmentPrefix <- function(x, ui) {
  .lhs <- x[[2]]
  stopifnot(length(.lhs) == 2)
  .state <- as.character(.lhs[2])
  .prefix <-
    list(
      f="f", # bioavailability
      F="f",
      alag="lag", # lag time
      lag="lag",
      rate="rate", # rate of infusion
      dur="dur" # duration of infusion
    )[[as.character(.lhs[[1]])]]
  if (is.null(.prefix)) {
    stop("unknown rxode2 assignment type:\n", deparse1(x),
         .call.=FALSE)
  }
  if (length(x[[3]]) != 1L) {
    stop("the complex initial condition is not supported in the nonmem conversion\n",
         deparse1(x), call.=FALSE)
  }
  .extra <- paste0(.rxToNonmem(x[[3]], ui=ui),
                   .babelmixr2Deparse(x))
  .nonmemSetCmtProperty(ui, .state, .extra, type=.prefix)
  paste0("; ", .prefix, "(", .state, ") defined in $PK block")
}

.rxToNonmemHandleCall <- function(x, ui) {
  if (identical(x[[1]], quote(`(`))) {
    return(paste0("(", .rxToNonmem(x[[2]], ui=ui), ")"))
  } else if (identical(x[[1]], quote(`{`))) {
    .x2 <- x[-1]
    .ret <- paste(lapply(.x2, function(x) {
      .rxToNonmem(x, ui=ui)
    }), collapse = "\n")
    return(.ret)
  } else if (.rxIsPossibleBinaryOperator(x[[1]])) {
    if (length(x) == 3) {
      return(.rxToNonmemHandleBinaryOperator(x, ui))
    } else {
      ## Unary Operators
      return(paste(
        as.character(x[[1]]),
        .rxToNonmem(x[[2]], ui=ui)
      ))
    }
  } else if (identical(x[[1]], quote(`if`))) {
    return(.rxToNonmemHandleIfExpressions(x, ui))
  } else if (.rxIsLogicalOperator(x[[1]])) {
    return(paste0(.rxToNonmem(x[[2]], ui=ui), .rxNMlogic[as.character(x[[1]])], .rxToNonmem(x[[3]], ui=ui)))
  } else if (identical(x[[1]], quote(`!`)) ) {
    return(paste0(".NOT. (", .rxToNonmem(x[[2]], ui=ui), ")"))
  } else if (.rxIsAssignmentOperator(x[[1]])) {
    return(.rxToNonmemHandleAssignmentOperator(x, ui))
  } else if (identical(x[[1]], quote(`[`))) {
    .type <- toupper(as.character(x[[2]]))
    if (any(.type == c("THETA", "ETA"))) {
      stop("'THETA'/'ETA' not supported by babelmixr2", call.=FALSE);
    }
  } else if (identical(x[[1]], quote(`log1pmx`))) {
      if (length(x == 2)) {
        .a <- .rxToNonmem(x[[2]], ui=ui)
        return(paste0("(DLOG(1+", .a, ")-(", .a, "))"))
      } else {
        stop("'log1pmx' only takes 1 argument", call. = FALSE)
      }
  } else if ((identical(x[[1]], quote(`pnorm`))) |
               (identical(x[[1]], quote(`normcdf`))) |
               (identical(x[[1]], quote(`phi`)))) {
    if (length(x) == 4) {
      .q <- .rxToNonmem(x[[2]], ui=ui)
      .mean <- .rxToNonmem(x[[3]], ui=ui)
      .sd <- .rxToNonmem(x[[4]], ui=ui)
      return(paste0("PHI(((", .q, ")-(", .mean, "))/(", .sd, "))"))
    } else if (length(x) == 3) {
      .q <- .rxToNonmem(x[[2]], ui=ui)
      .mean <- .rxToNonmem(x[[3]], ui=ui)
      return(paste0("PHI(((", .q, ")-(", .mean, ")))"))
    } else if (length(x) == 2) {
      .q <- .rxToNonmem(x[[2]], ui=ui)
      return(paste0("PHI(", .q, ")"))
    } else {
      stop("'pnorm' can only take 1-3 arguments", call. = FALSE)
    }
  } else {
    # handle single function translations
    if (length(x[[1]]) == 1) {
      .x1 <- as.character(x[[1]])
      .xc <- .rxNMsingle[[.x1]]
      if (!is.null(.xc)) {
        if (length(x) == 2) {
          if (.x1 %in% .rxNMprotectZero) {
            .expr <- .rxProtectPlusZero(x[[2]], ui=ui, one=FALSE)
          } else if (.x1 %in% .rxNmProtectZeroP1) {
            .expr <- .rxProtectPlusZero(x[[2]], ui=ui, one=TRUE)
          } else {
            .expr <- .rxToNonmem(x[[2]], ui=ui)
          }
          .ret <- paste0(
            .xc[1], .expr,
            .xc[2])
          if (.ret == "DEXP(1)") {
            return("2.718281828459045090796")
          }
          return(.ret)
        } else {
          stop(sprintf("'%s' only acceps 1 argument", .x1), call. = FALSE)
        }
      }
    }
    # There are no identical functions from nlmixr2 to NONMEM.
    .ret0 <- c(list(as.character(x[[1]])), lapply(x[-1], .rxToNonmem, ui=ui))
    .fun <- paste(.ret0[[1]])
    .ret0 <- .ret0[-1]
    .ret <- paste0("(", paste(unlist(.ret0), collapse = ","), ")")
    if (any(.fun == c("cmt", "dvid"))) {
      return("")
    } else if (any(.fun == c("max", "min"))) {
      .ret0 <- unlist(.ret0)
      .ret <- paste0(toupper(.fun), "(", paste(.ret0, collapse = ","), ")")
    } else if (.fun == "sum") {
      .ret <- paste0("(", paste(paste0("(", unlist(.ret0), ")"), collapse = "+"), ")")
    } else if (.fun == "prod") {
      .ret <- paste0("(", paste(paste0("(", unlist(.ret0), ")"), collapse = "*"), ")")
    } else if (.fun == "probitInv") {
      ##erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1 (probitInv=pnorm)
      if (length(.ret0) == 1) {
        .ret <- paste0("PHI(", unlist(.ret0)[1], ")")
      } else if (length(.ret0) == 2) {
        .ret0 <- unlist(.ret0)
        .p <- paste0("PHI(", .ret0[1], ")")
        ## return (high-low)*p+low;
        .low <- paste(.ret0[2])
        if (regexpr("^-?[0-9]+$", .low) != -1) .low <- paste0(.low, ".0")
        .ret <- paste0(
          "(1.0-(", .low, "))*(", .p,
          ")+(", .low, ")"
        )
      } else if (length(.ret0) == 3) {
        .ret0 <- unlist(.ret0)
        .low <- paste(.ret0[2])
        if (regexpr("^-?[0-9]+$", .low) != -1) .low <- paste0(.low, ".0")
        .hi <- paste(.ret0[3])
        if (regexpr("^-?[0-9]+$", .hi) != -1) .hi <- paste0(.hi, ".0")
        .p <- paste0("PHI(", .ret0[1], ")")
        .ret <- paste0(
          "((", .hi, ")-(", .low, "))*(", .p,
          ")+(", .low, ")"
        )
      } else {
        stop("'probitInv' requires 1-3 arguments",
             call. = FALSE
             )
      }
    } else if (.fun == "probit") {
      ##erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2) (probit=qnorm )
      stop("probit not supported in nonmem (though probitInv is supported)")
    } else if (.fun == "logit") {
      if (length(.ret0) == 1) {
        .ret <- paste0("-DLOG(1/(", unlist(.ret0), ")-1)")
      } else if (length(.ret0) == 2) {
        .ret0 <- unlist(.ret0)
        .low <- paste(.ret0[2])
        if (regexpr("^-?[0-9]+$", .low) != -1) .low <- paste0(.low, ".0")
        .p <- paste0(
          "((", .ret0[1], ")-(", .low, "))/(1.0-",
          "(", .low, "))"
        )
        .ret <- paste0("-DLOG(1/(", .p, ")-1)")
      } else if (length(.ret0) == 3) {
        .ret0 <- unlist(.ret0)
        .low <- paste(.ret0[2])
        if (regexpr("^-?[0-9]+$", .low) != -1) .low <- paste0(.low, ".0")
        .hi <- paste(.ret0[3])
        if (regexpr("^-?[0-9]+$", .hi) != -1) .hi <- paste0(.hi, ".0")
        ## (x-low)/(high-low)
        .p <- paste0(
          "((", .ret0[1], ")-(", .low,
          "))/((", .hi, ")-(", .low, "))"
        )
        .ret <- paste0("-DLOG(1/(", .p, ")-1)")
      } else {
        stop("'logit' requires 1-3 arguments",
             call. = FALSE
             )
      }
    } else if (any(.fun == c("expit", "invLogit", "logitInv"))) {
      if (length(.ret0) == 1) {
        .ret <- paste0("1/(1+DEXP(-(", unlist(.ret0)[1], ")))")
      } else if (length(.ret0) == 2) {
        .ret0 <- unlist(.ret0)
        .low <- paste(.ret0[2])
        if (regexpr("^-?[0-9]+$", .low) != -1) .low <- paste0(.low, ".0")
        .p <- paste0("1/(1+DEXP(-(", .ret0[1], ")))")
        ## return (high-low)*p+low;
        .ret <- paste0(
          "(1.0-(", .low, "))*(", .p,
          ")+(", .low, ")"
        )
      } else if (length(.ret0) == 3) {
        .ret0 <- unlist(.ret0)
        .p <- paste0("1/(1+DEXP(-(", .ret0[1], ")))")
        .low <- paste(.ret0[2])
        if (regexpr("^-?[0-9]+$", .low) != -1) .low <- paste0(.low, ".0")
        .hi <- paste(.ret0[3])
        if (regexpr("^-?[0-9]+$", .hi) != -1) .hi <- paste0(.hi, ".0")
        .ret <- paste0(
          "((", .hi, ")-(", .low, "))*(", .p,
          ")+(", .low, ")"
        )
      } else {
        stop("'expit' requires 1-3 arguments",
             call. = FALSE)
      }
    } else {
      stop(sprintf(gettext("function '%s' is not supported in NONMEM<->nlmixr"), .fun),
           call. = FALSE)
    }
  }
}

.rxIsAssignmentOperator <- function(expr) {
  identical(expr, quote(`=`)) ||
    identical(expr, quote(`<-`)) ||
    identical(expr, quote(`~`))
}

.rxToNonmem <- function(x, ui) {
  if (is.name(x) || is.atomic(x)) {
    .rxToNonmemHandleNamesOrAtomic(x, ui)
  } else if (is.call(x)) {
    .rxToNonmemHandleCall(x, ui)
  } else {
    stop("unrecognized model part")
  }
}

#' Convert RxODE syntax to NONMEM syntax
#'
#' @param x Expression
#' @param ui rxode2 ui
#' @return NONMEM syntax
#' @author Matthew Fidler
#' @export
rxToNonmem <- function(x, ui) {
  ui <- rxode2::assertRxUi(ui)
  ui <- rxode2::rxUiDecompress(ui)
  if (is(substitute(x), "character")) {
    force(x)
  } else if (is(substitute(x), "{")) {
    x <- deparse1(substitute(x))
    if (x[1] == "{") {
      x <- x[-1]
      x <- x[-length(x)]
    }
    x <- paste(x, collapse = "\n")
  } else {
    .xc <- as.character(substitute(x))
    x <- substitute(x)
    if (length(.xc == 1)) {
      .found <- FALSE
      .frames <- seq(1, sys.nframe())
      .frames <- .frames[.frames != 0]
      for (.f in .frames) {
        .env <- parent.frame(.f)
        if (exists(.xc, envir = .env)) {
          .val2 <- try(get(.xc, envir = .env), silent = TRUE)
          if (inherits(.val2, "character")) {
            .val2 <- eval(parse(text = paste0("quote({", .val2, "})")))
            return(.rxToNonmem(.val2, ui=ui))
          } else if (inherits(.val2, "numeric") || inherits(.val2, "integer")) {
            return(sprintf("%s", .val2))
          }
        }
      }
    } else {
      stop("too many lines to parse")
    }
    return(.rxToNonmem(x, ui=ui))
  }
  return(.rxToNonmem(eval(parse(text = paste0("quote({", x, "})"))),
                      ui=ui))
}
