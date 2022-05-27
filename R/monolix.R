##
## observationTypes (list): A list giving the type of each observation present in the data file. If there is only one y-type, the corresponding observation name can be omitted.
## The possible observation types are "continuous", "discrete", and "event".
##
## nbSSDoses [optional](int): Number of doses (if there is a SS column).

.monolixErrs <- c()


################################################################################
# https://cran.r-project.org/doc/manuals/r-release/R-exts.html
.rxMcnt <- c(
  # "band"
  # "bsmm"
  # "categorical"
  # "categories"
  ########
  # "amtDose"
  # "inftDose"
  tlast="tDose",
  time="t",
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

.rxMbad <- c("NA", "NaN", "Inf", "newind", "NEWIND")

.rxMbadF <- c("digamma", "trigamma", "tetragamma", "pentagamma", "psigamma", "choose", "lchoose")

.rxMsingle <- list(
  "gammafn" = c("exp(gammaln(", "))"),
  "lgammafn" = c("gammaln(", ")"),
  "lgamma" = c("gammaln(", ")"),
  "loggamma" = c("gammaln(", ")"),
  "cospi" = c("cos(3.141592653589793115998*(", "))"),
  "sinpi" = c("sin(3.141592653589793115998*(", "))"),
  "tanpi" = c("tan(3.141592653589793115998*(", "))"),
  "log1p" = c("log(1+", ")"),
  "expm1" = c("(exp(", ")-1)"),
  "lfactorial" = c("factln(", ")"),
  "lgamma1p" = c("gammaln(", "+1)"),
  "expm1" = c("(exp(", ")-1)"),
  "log10" = c("log10(", ")"),
  "log2" = c("(log(", ")*1.442695040888963387005)"),
  "log1pexp" = c("log(1+exp(", "))", "log1pexp"),
  "phi" = c("normcdf(", ")"),
  "pnorm" = c("normcdf(", ")"),
  "qnorm"=c("probit(", ")"),
  "fabs"=c("abs(", ")")
)

.rxMeq <- c("sqrt"=1,
            "exp"=1,
            "abs"=1,
            "log"=1,
            "log10"=1,
            "normcdf"=1,
            "sin"=1,
            "cos"=1,
            "tan"=1,
            "asin"=1,
            "acos"=1,
            "atan"=1,
            "sinh"=1,
            "cosh"=1,
            "atan2"=2,
            "floor"=1,
            "ceil"=1,
            "factorial"=1
            )

#' Get the monolix administration info
#'
#'
#' @param ui rxode2 user interface
#' @return internal adm dataset that comes from `bblDatToMonolix()`
#' @author Matthew L. Fidler
#' @noRd
.monolixGetAdm <- function(ui) {
  rxode2::rxGetControl(ui, ".adm",
                       data.frame(adm=1L,
                                  cmt=1L,
                                  type=factor("bolus", levels=c("empty", "modelRate", "modelDur", "infusion", "bolus")),
                                  f=NA_character_,
                                  dur=NA_character_,
                                  lag=NA_character_,
                                  rate=NA_character_))
}

#' Set all the administration types for each cmt for monolix
#'
#' @param ui rxode2 user interface
#' @param state which value this administration property is applied to
#' @param param The parameter in nlmixr that is applying this effect
#' @param type the type of administration
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.monolixSetAdm <- function(ui, state, param, type="f") {
  .adm <- .monolixGetAdm(ui)
  .state <- rxode2::rxState(ui)
  .cmt <- which(state == .state)
  .w <- which(.adm$cmt == .cmt)
  if (length(.w) == 0L) return(invisible())
  if (type == "f") {
    .adm[.w, "f"] <- param
  } else if (type == "dur") {
    .adm[.w, "dur"] <- param
  } else if (type == "lag") {
    .adm[.w, "lag"] <- param
  } else if (type == "rate") {
    .adm[.w, "rate"] <- param
  }
  rxode2::rxAssignControlValue(ui, ".adm", .adm)
}

.rxToMonolixHandleBinaryOperator <- function(x, ui) {
  if (identical(x[[1]], quote(`/`))) {
    .x2 <- x[[2]]
    .x3 <- x[[3]]
    ## df(%s)/dy(%s)
    if (identical(.x2, quote(`d`)) &&
          identical(.x3[[1]], quote(`dt`))) {
      if (length(.x3[[2]]) == 1) {
        .state <- as.character(.x3[[2]])
      } else {
        .state <- .rxToMonolix(.x3[[2]], ui=ui)
      }
      return(paste0("ddt_", .state))
    } else {
      if (length(.x2) == 2 && length(.x3) == 2) {
        if (identical(.x2[[1]], quote(`df`)) &&
              identical(.x3[[1]], quote(`dy`))) {
          stop('df()/dy() is not supported in monolix conversion', call.=FALSE)
        }
      }
      .ret <- paste0(
        .rxToMonolix(.x2, ui=ui),
        as.character(x[[1]]),
        .rxToMonolix(.x3, ui=ui)
      )
    }
  } else {
    .ret <- paste0(
      .rxToMonolix(x[[2]], ui=ui),
      as.character(x[[1]]),
      .rxToMonolix(x[[3]], ui=ui)
    )
  }
  return(.ret)
}

.rxToMonolixIndent <- function(ui) {
  rxode2::rxAssignControlValue(ui, ".mIndent",
                               rxode2::rxGetControl(ui, ".mIndent", 0) + 2)
}

.rxToMonolixUnIndent <- function(ui) {
  rxode2::rxAssignControlValue(ui, ".mIndent",
                               max(0, rxode2::rxGetControl(ui, ".mIndent", 0) - 2))
}

.rxToMonolixGetIndent <- function(ui, ind=NA) {
  if (is.na(ind)) {
  } else if (ind) {
    .rxToMonolixIndent(ui)
  } else {
    .rxToMonolixUnIndent(ui)
  }
  .nindent <- rxode2::rxGetControl(ui, ".mIndent", 0)
  paste(vapply(seq(1, .nindent), function(x) " ", character(1), USE.NAMES=FALSE), collapse="")
}


.rxToMonolixHandleIfExpressions <- function(x, ui) {
  .ret <- paste0(.rxToMonolixGetIndent(ui), "if ", .rxToMonolix(x[[2]], ui=ui), "\n")
  .rxToMonolixIndent(ui)
  .ret <- paste0(.ret, .rxToMonolix(x[[3]], ui=ui))
  x <- x[-c(1:3)]
  if (length(x) == 1) x <- x[[1]]
  while(identical(x[[1]], quote(`if`))) {
    .ret <- paste0(.ret, "\n",
                   .rxToMonolixGetIndent(ui, FALSE), "elseif ", .rxToMonolix(x[[2]], ui=ui), "\n")
    .rxToMonolixIndent(ui)
    .ret <- paste0(.ret, .rxToMonolix(x[[3]], ui=ui))
    x <- x[-c(1:3)]
    if (length(x) == 1) x <- x[[1]]
  }
  if (is.null(x)) {
    .ret <- paste0(.ret, "\n",
                   .rxToMonolixGetIndent(ui, FALSE), "end\n")
  }  else {
    .ret <- paste0(.ret, "\n",
                   .rxToMonolixGetIndent(ui, FALSE), "else\n")
    .rxToMonolixIndent(ui)
    .ret <- paste0(.ret, .rxToMonolix(x, ui=ui),
                   "\n",
                   .rxToMonolixGetIndent(ui, FALSE), "end\n")
  }
  return(.ret)
}

.rxToMonolix <- function(x, ui) {
  if (is.name(x) || is.atomic(x)) {
    if (is.character(x)) {
      stop("strings in nlmixr<->monolix are not supported", call.=FALSE)
    } else {
      .ret <- as.character(x)
      if (is.na(.ret) | (.ret %in% .rxMbad)) {
        stop("'", .ret, "' cannot be translated to monolix", call.=FALSE)
      }
      .v <- .rxMcnt[.ret]
      if (is.na(.v)) {
        if (is.numeric(.ret)) {
          return(.ret)
        } else if (regexpr("^(?:-)?(?:(?:0|(?:[1-9][0-9]*))|(?:(?:[0-9]+\\.[0-9]*)|(?:[0-9]*\\.[0-9]+))(?:(?:[Ee](?:[+\\-])?[0-9]+))?|[0-9]+[Ee](?:[\\-+])?[0-9]+)$",
                           .ret, perl=TRUE) != -1) {
          return(.ret)
        } else {
          return(gsub("[.]", "__", .ret))
        }
      } else {
        return(.v)
      }
    }
  } else if (is.call(x)) {
    if (identical(x[[1]], quote(`(`))) {
      return(paste0("(", .rxToMonolix(x[[2]], ui=ui), ")"))
    } else if (identical(x[[1]], quote(`{`))) {
      .x2 <- x[-1]
      .ret <- paste(lapply(.x2, function(x) {
        .rxToMonolix(x, ui=ui)
      }), collapse = "\n")
      return(.ret)
    } else if (.rxIsPossibleBinaryOperator(x[[1]])) {
      if (length(x) == 3) {
        return(.rxToMonolixHandleBinaryOperator(x, ui))
      } else {
        ## Unary Operators
        return(paste(
          as.character(x[[1]]),
          .rxToMonolix(x[[2]], ui=ui)
        ))
      }
    } else if (identical(x[[1]], quote(`if`))) {
      return(.rxToMonolixHandleIfExpressions(x, ui))
    } else if (.rxIsLogicalOperator(x[[1]])) {
        ## Use "preferred" monolix syntax
      return(paste0(.rxToMonolix(x[[2]], ui=ui), as.character(x[[1]]), .rxToMonolix(x[[3]], ui=ui)))
    } else if (identical(x[[1]], quote(`!`)) ) {
      return(paste0("~", .rxToMonolix(x[[2]], ui=ui)))
    } else if (identical(x[[1]], quote(`**`)) ) {
      return(paste(.rxToMonolix(x[[2]], ui=ui), "^", .rxToMonolix(x[[3]], ui=ui)))
    } else if (.rxIsAssignmentOperator(x[[1]])) {
      if (any(as.character(x[[2]])[1] == c("alag", "lag", "F", "f", "rate", "dur"))) {
        if (any(as.character(x[[2]])[1] == c("alag", "lag"))) {
          .state <- as.character(x[[2]][[2]])
          if (length(x[[3]]) == 1L) {
            .extra <- .rxToMonolix(x[[3]], ui=ui)
            .monolixSetAdm(ui, .state, .extra, type="lag")
          } else {
            stop("the complex lag time is not supported by babelmixr2",
                 call.=FALSE)
          }
        }
        if (any(as.character(x[[2]])[1] == c("F", "f"))) {
          .state <- as.character(x[[2]][[2]])
          if (length(x[[3]]) == 1L) {
            .extra <- .rxToMonolix(x[[3]], ui=ui)
            .monolixSetAdm(ui, .state, .extra, type="f")
          } else {
            stop("the complex F is not supported by babelmixr2",
                 call.=FALSE)
          }
        }
        if (as.character(x[[2]])[1] == "rate") {
          .state <- as.character(x[[2]][[2]])
          if (length(x[[3]]) == 1L) {
            .extra <- .rxToMonolix(x[[3]], ui=ui)
            .monolixSetAdm(ui, .state, .extra, type="rate")
          } else {
            stop("the complex rate is not supported by babelmixr2",
                 call.=FALSE)
          }
        }
        if (as.character(x[[2]])[1] == "dur") {
          .state <- as.character(x[[2]][[2]])
          if (length(x[[3]]) == 1L) {
            .extra <- .rxToMonolix(x[[3]], ui=ui)
            .monolixSetAdm(ui, .state, .extra, type="dur")
          } else {
            stop("the complex dur is not supported by babelmixr2",
                 call.=FALSE)
          }
        }
        return(paste0(.rxToMonolixGetIndent(ui),
                      ";", as.character(x[[2]])[1], " defined in PK section"))
      }
      .var <- .rxToMonolix(x[[2]], ui=ui)
      return(paste(.rxToMonolixGetIndent(ui),
                   .var, "=", .rxToMonolix(x[[3]], ui=ui)))
    } else if (identical(x[[1]], quote(`[`))) {
      .type <- toupper(as.character(x[[2]]))
      if (any(.type == c("THETA", "ETA"))) {
        stop("'THETA'/'ETA' not supported by monolix", call.=FALSE);
      }
    } else if (identical(x[[1]], quote(`log1pmx`))) {
      if (length(x == 2)) {
        .a <- .rxToMonolix(x[[2]], ui=ui)
        return(paste0("(log(1+", .a, ")-(", .a, "))"))
      } else {
        stop("'log1pmx' only takes 1 argument", call. = FALSE)
      }
    } else if ((identical(x[[1]], quote(`pnorm`))) |
                 (identical(x[[1]], quote(`normcdf`))) |
                 (identical(x[[1]], quote(`phi`)))) {
      if (length(x) == 4) {
        .q <- .rxToMonolix(x[[2]], ui=ui)
        .mean <- .rxToMonolix(x[[3]], ui=ui)
        .sd <- .rxToMonolix(x[[4]], ui=ui)
        return(paste0("normcdf(((", .q, ")-(", .mean, "))/(", .sd, "))"))
      } else if (length(x) == 3) {
        .q <- .rxToMonolix(x[[2]], ui=ui)
        .mean <- .rxToMonolix(x[[3]], ui=ui)
        return(paste0("normcdf(((", .q, ")-(", .mean, ")))"))
      } else if (length(x) == 2) {
        .q <- .rxToMonolix(x[[2]], ui=ui)
        return(paste0("normcdf(", .q, ")"))
      } else {
        stop("'pnorm' can only take 1-3 arguments", call. = FALSE)
      }
    } else {
      if (length(x[[1]]) == 1) {
        .x1 <- as.character(x[[1]])
        .xc <- .rxMsingle[[.x1]]
        if (!is.null(.xc)) {
          if (length(x) == 2) {
            .ret <- paste0(
              .xc[1], .rxToMonolix(x[[2]], ui=ui),
              .xc[2])
            return(.ret)
          } else {
            stop(sprintf("'%s' only acceps 1 argument", .x1), call. = FALSE)
          }
        }
      }
      .ret0 <- c(list(as.character(x[[1]])), lapply(x[-1], .rxToMonolix, ui=ui))
      .SEeq <- .rxMeq
      .curName <- paste(.ret0[[1]])
      .nargs <- .SEeq[.curName]
      if (!is.na(.nargs)) {
        if (.nargs == length(.ret0) - 1) {
          .ret <- paste0(.ret0[[1]], "(")
          .ret0 <- .ret0[-1]
          .ret <- paste0(.ret, paste(unlist(.ret0), collapse = ","), ")")
          if (.ret == "exp(1)") {
            return("2.718281828459045090796")
          }
          return(.ret)
        } else {
          stop(sprintf(
            gettext("'%s' takes %s arguments (has %s)"),
            paste(.ret0[[1]]),
            .nargs, length(.ret0) - 1
          ), call. = FALSE)
        }
      } else {
        .fun <- paste(.ret0[[1]])
        .ret0 <- .ret0[-1]
        .ret <- paste0("(", paste(unlist(.ret0), collapse = ","), ")")
        if (.ret == "(0)") {
          return(paste0(.fun, "_0"))
        } else if (any(.fun == c("cmt", "dvid"))) {
          return("")
        } else if (any(.fun == c("max", "min"))) {
          ## Not sure but I think that max/min only supports 2 arguments in monolix
          .ret0 <- unlist(.ret0)
          if (length(.ret0) != 2) {
            stop("'", .fun, "' in monolix can only have 2 arguments", call.=FALSE)
          }
          .ret <- paste0(.fun, "(", paste(.ret0, collapse = ","), ")")
        } else if (.fun == "sum") {
          .ret <- paste0("(", paste(paste0("(", unlist(.ret0), ")"), collapse = "+"), ")")
        } else if (.fun == "prod") {
          .ret <- paste0("(", paste(paste0("(", unlist(.ret0), ")"), collapse = "*"), ")")
        } else if (.fun == "probitInv") {
          ##erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1 (probitInv=pnorm)
          if (length(.ret0) == 1) {
            .ret <- paste0("normcdf(", unlist(.ret0)[1], ")")
          } else if (length(.ret0) == 2) {
            .ret0 <- unlist(.ret0)
            .p <- paste0("normcdf(", .ret0[1], ")")
            ## return (high-low)*p+low;
            .ret <- paste0(
              "(1.0-(", .ret0[2], "))*(", .p,
              ")+(", .ret0[2], ")"
            )
          } else if (length(.ret0) == 3) {
            .ret0 <- unlist(.ret0)
            .p <- paste0("normcdf(", .ret0[1], ")")
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
          if (length(.ret0) == 1) {
            .ret <- paste0("probit(", unlist(.ret0), ")")
          } else if (length(.ret0) == 2) {
            .ret0 <- unlist(.ret0)
            .p <- paste0(
              "((", .ret0[1], ")-(", .ret0[2], "))/(1.0-",
              "(", .ret0[2], "))"
            )
            .ret <- paste0("probit(", .p, ")")
          } else if (length(.ret0) == 3) {
            .ret0 <- unlist(.ret0)
            ## (x-low)/(high-low)
            .p <- paste0(
              "((", .ret0[1], ")-(", .ret0[2],
              "))/((", .ret0[3], ")-(", .ret0[2], "))"
            )
            .ret <- paste0("probit(", .p, ")")
          } else {
            stop("'probit' requires 1-3 arguments",
              call. = FALSE
            )
          }
        } else if (.fun == "logit") {
          if (length(.ret0) == 1) {
            .ret <- paste0("-log(1/(", unlist(.ret0), ")-1)")
          } else if (length(.ret0) == 2) {
            .ret0 <- unlist(.ret0)
            .p <- paste0(
              "((", .ret0[1], ")-(", .ret0[2], "))/(1.0-",
              "(", .ret0[2], "))"
            )
            .ret <- paste0("-log(1/(", .p, ")-1)")
          } else if (length(.ret0) == 3) {
            .ret0 <- unlist(.ret0)
            ## (x-low)/(high-low)
            .p <- paste0(
              "((", .ret0[1], ")-(", .ret0[2],
              "))/((", .ret0[3], ")-(", .ret0[2], "))"
            )
            .ret <- paste0("-log(1/(", .p, ")-1)")
          } else {
            stop("'logit' requires 1-3 arguments",
              call. = FALSE
            )
          }
        } else if (any(.fun == c("expit", "invLogit", "logitInv"))) {
          if (length(.ret0) == 1) {
            .ret <- paste0("1/(1+exp(-(", unlist(.ret0)[1], ")))")
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
            .p <- paste0("1/(1+exp(-(", .ret0[1], ")))")
            .ret <- paste0(
              "((", .ret0[3], ")-(", .ret0[2], "))*(", .p,
              ")+(", .ret0[2], ")"
            )
          } else {
            stop("'expit' requires 1-3 arguments",
              call. = FALSE
            )
          }
        } else {
          stop(sprintf(gettext("function '%s' is not supported in monolix<->nlmixr"), .fun),
            call. = FALSE
          )
        }
      }
    }
  }
}

#' Convert RxODE syntax to monolix syntax
#'
#' @param x Expression
#' @param ui rxode2 ui
#' @return Monolix syntax
#' @author Matthew Fidler
#' @export
rxToMonolix <- function(x, ui) {
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
            return(.rxToMonolix(.val2, ui=ui))
          } else if (inherits(.val2, "numeric") || inherits(.val2, "integer")) {
            return(sprintf("%s", .val2))
          }
        }
      }
    }
    return(.rxToMonolix(x, ui=ui))
  }
  return(.rxToMonolix(eval(parse(text = paste0("quote({", x, "})"))),
                      ui=ui))
}

