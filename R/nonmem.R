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

#' Handle numbers and symbols
#'
#' @param x Expression
#' @param ui User interface
#' @return Symbol, converted to NONMEM compatible name
#' @author Matthew L. Fidler
#' @noRd
.rxToNonmemHandleNamesOrAtomic <- function(x, ui=NULL) {
  if (is.character(x)) stop("strings in nlmixr<->monolix are not supported", call.=FALSE)
  .ret <- as.character(x)
  if (is.na(.ret) | (.ret %in% .rxNMbad)) {
    stop("'", .ret, "' cannot be translated to NONMEM", call.=FALSE)
  }
  .v <- .rxNMcnt[.ret]
  if (is.na(.v)) {
    if (is.numeric(.ret)) {
      return(.ret)
    } else if (regexpr("(?:-)?(?:(?:0|(?:[1-9][0-9]*))|(?:(?:[0-9]+\\.[0-9]*)|(?:[0-9]*\\.[0-9]+))(?:(?:[Ee](?:[+\\-])?[0-9]+))?|[0-9]+[Ee](?:[\\-+])?[0-9]+)",
                       .ret, perl=TRUE) != -1) {
      return(.ret)
    } else {
      .states <- rxode2::rxModelVars(ui)$state
      .w <- which(.ret == .states)
      if (length(.w) == 1) {
        return(paste0("A(", .w, ")"))
      }
      return(gsub("[.]", "__", toupper(.ret)))
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

.rxToNonmemHandleBinaryOperator <- function(x, ui) {
  if (identical(x[[1]], quote(`/`))) {
    .x2 <- x[[2]]
    .x3 <- x[[3]]
    ## df(%s)/dy(%s)
    if (identical(.x2, quote(`d`)) &&
          identical(.x3[[1]], quote(`dt`))) {
      if (length(.x3[[2]]) == 1) {
        .state <- as.character(.x3[[2]])
      } else {
        .state <- .rxToNonmem(.x3[[2]], ui=ui)
      }
      .states <- rxode2::rxModelVars(ui)$state
      .num <- which(.state == .states)
      return(paste0("DADT(", .num, ")"))
    } else {
      if (length(.x2) == 2 && length(.x3) == 2) {
        if (identical(.x2[[1]], quote(`df`)) &&
              identical(.x3[[1]], quote(`dy`))) {
          stop('df()/dy() is not supported in NONMEM conversion', call.=FALSE)
        }
      }
      .ret <- paste0(
        .rxToNonmem(.x2, ui=ui),
        .rxNMbin[as.character(x[[1]])],
        .rxToNonmem(.x3, ui=ui)
      )
    }
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


.rxToNonmemIndent <- function(ui) {

}

.rxToNonmemUnIndent <- function(ui) {

}

.rxToNonmemGetIndent <- function(ui) {

}

.rxToNonmemHandleIfExpressions <- function(x, ui) {
  .ret <- paste0("IF ", .rxToNonmem(x[[2]], ui=ui), " THEN\n",
                 "  ", .rxToNonmem(x[[3]], ui=ui))
  x <- x[-c(1:3)]
  if (length(x) == 1) x <- x[[1]]
  while(identical(x[[1]], quote(`if`))) {
    .ret <- paste0(.ret, "\nELSE IF ", .rxToNonmem(x[[2]], ui=ui), " THEN\n",
                   "  ", .rxToNonmem(x[[3]], ui=ui))
    x <- x[-c(1:3)]
    if (length(x) == 1) x <- x[[1]]
  }
  if (is.null(x)) {
    .ret <- paste0(.ret, "\nEND IF\n")
  }  else {
    .ret <- paste0(.ret, "\nELSE \n",
                   "  ", .rxToNonmem(x, ui=ui),
                   "\nEND IF\n")
  }
  return(.ret)
}

.nonmemGetCmtProperties <- function(ui) {
  rxode2::rxGetControl(ui, ".cmtProperties",
                       data.frame(cmt=integer(0),
                                  f=character(0),
                                  dur=character(0),
                                  lag=character(0),
                                  rate=character(0)))
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
                              rate=NA_character_))
    .w <- which(.prop$cmt == .cmt)
  }
  if (type == "f") {
    .prop[.w, "f"] <- param
  } else if (type == "dur") {
    .prop[.w, "dur"] <- param
  } else if (type == "lag") {
    .prop[.w, "lag"] <- param
  } else if (type == "rate") {
    .prop[.w, "rate"] <- param
  }
  rxode2::rxAssignControlValue(ui, ".cmtProperties", .adm)
}

.rxToNonmemHandleAssignmentOperator <- function(x, ui) {
  if (any(as.character(x[[2]])[1] == c("alag", "lag", "F", "f", "rate", "dur"))) {
    # I believe these have to be in the $PK block
    if (any(as.character(x[[2]])[1] == c("alag", "lag"))) {
      .state <- as.character(x[[2]][[2]])
      if (length(x[[3]]) == 1L) {
        .extra <- .rxToNonmem(x[[3]], ui=ui)
        .nonmemSetCmtProperty(ui, .state, .extra, type="lag")
      } else {
        stop("the complex lag time is not supported by babelmixr2",
             call.=FALSE)
      }
    }
    if (any(as.character(x[[2]])[1] == c("F", "f"))) {
      .state <- as.character(x[[2]][[2]])
      if (length(x[[3]]) == 1L) {
        .extra <- .rxToNonmem(x[[3]], ui=ui)
        .nonmemSetCmtProperty(ui, .state, .extra, type="f")
      } else {
        stop("the complex F is not supported by babelmixr2",
             call.=FALSE)
      }
    }
    if (as.character(x[[2]])[1] == "rate") {
      .state <- as.character(x[[2]][[2]])
      if (length(x[[3]]) == 1L) {
        .extra <- .rxToNonmem(x[[3]], ui=ui)
        .nonmemSetCmtProperty(ui, .state, .extra, type="rate")
      } else {
        stop("the complex rate is not supported by babelmixr2",
             call.=FALSE)
      }
    }
    if (as.character(x[[2]])[1] == "dur") {
      .state <- as.character(x[[2]][[2]])
      if (length(x[[3]]) == 1L) {
        .extra <- .rxToNonmem(x[[3]], ui=ui)
        .nonmemSetCmtProperty(ui, .state, .extra, type="dur")
      } else {
        stop("the complex dur is not supported by babelmixr2",
             call.=FALSE)
      }
    }
    return(paste0(";", as.character(x[[2]])[1], " defined in $PK block"))
  }
  .var <- .rxToNonmem(x[[2]], ui=ui)
  return(paste(.var, "=", .rxToNonmem(x[[3]], ui=ui)))
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
  }  else if (identical(x[[1]], quote(`if`))) {
    return(.rxToNonmemHandleIfExpressions(x, ui))
  } else if (.rxIsLogicalOperator(x[[1]])) {
    return(paste0(.rxToNonmem(x[[2]], ui=ui), .rxNMlogic[as.character(x[[1]])], .rxToNonmem(x[[3]], ui=ui)))
  } else if (identical(x[[1]], quote(`!`)) ) {
    return(paste0(".NOT. (", .rxToNonmem(x[[2]], ui=ui), ")"))
  } else if (.rxIsAssignmentOperator(x[[1]])) {
    return(.rxToNonmemHandleAssignmentOperator(x))
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
          .ret <- paste0(
            .xc[1], .rxToNonmem(x[[2]], ui=ui),
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
    # There are no identical functions from nlmixr to NONMEM.
    .fun <- paste(.ret0[[1]])
    .ret0 <- .ret0[-1]
    .ret <- paste0("(", paste(unlist(.ret0), collapse = ","), ")")
    if (.ret == "(0)") {
      .state <- rxode2::rxModelVars(ui)
      .cmt <- which(.fun == .state)
      return(paste0("A_0(", .cmt, ")"))
    } else if (any(.fun == c("cmt", "dvid"))) {
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
        .ret <- paste0(
          "(1.0-(", .ret0[2], "))*(", .p,
          ")+(", .ret0[2], ")"
        )
      } else if (length(.ret0) == 3) {
        .ret0 <- unlist(.ret0)
        .p <- paste0("PHI(", .ret0[1], ")")
        .ret <- paste0(
          "((", .ret0[3], ")-(", .ret0[2], "))*(", .p,
          ")+(", .ret0[2], ")"
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
        .p <- paste0(
          "((", .ret0[1], ")-(", .ret0[2], "))/(1.0-",
          "(", .ret0[2], "))"
        )
        .ret <- paste0("-DLOG(1/(", .p, ")-1)")
      } else if (length(.ret0) == 3) {
        .ret0 <- unlist(.ret0)
        ## (x-low)/(high-low)
        .p <- paste0(
          "((", .ret0[1], ")-(", .ret0[2],
          "))/((", .ret0[3], ")-(", .ret0[2], "))"
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
        .p <- paste0("1/(1+exp(-(", .ret0[1], ")))")
        ## return (high-low)*p+low;
        .ret <- paste0(
          "(1.0-(", .ret0[2], "))*(", .p,
          ")+(", .ret0[2], ")"
        )
      } else if (length(.ret0) == 3) {
        .ret0 <- unlist(.ret0)
        .p <- paste0("1/(1+DEXP(-(", .ret0[1], ")))")
        .ret <- paste0(
          "((", .ret0[3], ")-(", .ret0[2], "))*(", .p,
          ")+(", .ret0[2], ")"
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
    return(.rxToNonmemHandleNamesOrAtomic(x))
  } else if (is.call(x)) {
    return(.rxToNonmemHandleCall(x, ui))
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
    }
    return(.rxToNonmem(x, ui=ui))
  }
  return(.rxToNonmem(eval(parse(text = paste0("quote({", x, "})"))),
                      ui=ui))
}
