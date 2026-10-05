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
  # probitInv() is normalized to erf(); erf(x) = 2*normcdf(sqrt(2)*x) - 1
  "erf"=c("(2*normcdf(1.414213562373095145475*(", "))-1)"),
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
  .cmt <- .rxGetCmtNumber(state, ui, error=FALSE)
  if (is.na(.cmt)) return(invisible())
  # the administrations dosing this compartment; a property of a
  # compartment without doses has no effect
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

#' Defaults of generated compartment property variables
#'
#' Lines like `rx_f_depot = 1` that go at the top of Monolix's
#' `EQUATION:` so a property only set inside a conditional keeps
#' rxode2's default otherwise.
#'
#' @noRd
.monolixCmtPropDefaults <- NULL

#' Translate a compartment property (f, alag, rate, dur) to Monolix
#'
#' Monolix sets these in the `PK:` macros, which take a single
#' variable.  The property is assigned in `EQUATION:` to a generated
#' variable (like `rx_f_depot`) that the macro then uses, so
#' expressions (issue #115) and conditional assignments translate.
#' When the property is set inside a conditional, the variable starts
#' at rxode2's default (1 for `f()`, 0 for `alag()`).
#'
#' @param x assignment expression, like `f(depot) <- exp(lfdepot)`
#' @param prop the property name (`"f"`, `"F"`, `"alag"`, `"lag"`,
#'   `"rate"` or `"dur"`)
#' @param ui rxode2 ui
#' @return Monolix `EQUATION:` line(s)
#' @author Matthew L. Fidler
#' @noRd
.rxToMonolixCmtProp <- function(x, prop, ui) {
  .type <- switch(prop, alag = "lag", F = "f", prop)
  .state <- as.character(x[[2]][[2]])
  .var <- paste0("rx_", .type, "_", gsub("[.]", "__", .state))
  .monolixSetAdm(ui, .state, .var, type = .type)
  .default <- switch(.type, f = "1", lag = "0", NA_character_)
  if (
    !is.na(.default) &&
      rxode2::rxGetControl(ui, ".mIndent", 0) > 0
  ) {
    assignInMyNamespace(
      ".monolixCmtPropDefaults",
      unique(c(.monolixCmtPropDefaults, paste0("   ", .var, " = ", .default)))
    )
  }
  .val <- .rxToMonolix(x[[3]], ui = ui)
  paste0(
    .rxToMonolixGetIndent(ui),
    ";",
    prop,
    " defined in PK section\n",
    .rxToMonolixFlushPrefixLines(ui),
    paste(.rxToMonolixGetIndent(ui), .var, "=", .val)
  )
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
  # top-level statements (.mIndent 0) are indented 2 spaces, like nested ones
  strrep(" ", max(2, .nindent))
}


#' Get (and clear) the lines that need to be written before a line
#'
#' @param ui rxode2 ui
#' @return the prefix lines followed by a new line, or `""` when there
#'   are no prefix lines
#' @noRd
#' @author Matthew L. Fidler
.rxToMonolixFlushPrefixLines <- function(ui) {
  # the indicator variables are only reused within the same line
  rxode2::rxAssignControlValue(ui, ".mLogicalDf", NULL)
  .prefixLines <- rxode2::rxGetControl(ui, ".mPrefixLines", NULL)
  if (is.null(.prefixLines)) {
    return("")
  }
  rxode2::rxAssignControlValue(ui, ".mPrefixLines", NULL)
  paste0(paste(.prefixLines, collapse = "\n"), "\n")
}

#' Translate a condition (like in `if ()`) to Monolix
#'
#' @param x rxode2 condition expression
#' @param ui rxode2 ui
#' @return Monolix logical expression
#' @noRd
#' @author Matthew L. Fidler
.rxToMonolixCondition <- function(x, ui) {
  if (is.call(x)) {
    if (identical(x[[1]], quote(`(`))) {
      return(paste0("(", .rxToMonolixCondition(x[[2]], ui), ")"))
    }
    if (identical(x[[1]], quote(`!`))) {
      return(paste0("~", .rxToMonolixCondition(x[[2]], ui)))
    }
    .op <- as.character(x[[1]])
    if (length(x) == 3L && .op %in% c("&&", "||", "&", "|")) {
      return(paste0(
        .rxToMonolixCondition(x[[2]], ui),
        .op,
        .rxToMonolixCondition(x[[3]], ui)
      ))
    }
    if (length(x) == 3L && .rxIsLogicalOperator(x[[1]])) {
      ## Use "preferred" monolix syntax
      return(paste0(
        .rxToMonolix(x[[2]], ui = ui),
        .op,
        .rxToMonolix(x[[3]], ui = ui)
      ))
    }
  }
  .rxToMonolix(x, ui = ui)
}

#' Translate a logical expression used as a number to Monolix
#'
#' Pruned `if`/`else` branches use logical expressions (like
#' `(WT > 70)`) as numbers.  This writes a 0/1 indicator variable in
#' the lines before the current line and uses that variable instead.
#'
#' @param x rxode2 logical expression
#' @param ui rxode2 ui
#' @return Monolix indicator variable name
#' @noRd
#' @author Matthew L. Fidler
.rxToMonolixLogicalIndicator <- function(x, ui) {
  .cond <- .rxToMonolixCondition(x, ui)
  .df <- rxode2::rxGetControl(
    ui,
    ".mLogicalDf",
    data.frame(cond = character(0), nm = character(0))
  )
  .w <- which(.df$cond == .cond)
  if (length(.w) == 1L) {
    return(.df$nm[.w])
  }
  .num <- rxode2::rxGetControl(ui, ".mVarLNum", 1)
  .newVar <- sprintf("rx_l%03d", .num)
  rxode2::rxAssignControlValue(ui, ".mVarLNum", .num + 1)
  .indent <- .rxToMonolixGetIndent(ui)
  .prefixLines <- c(
    rxode2::rxGetControl(ui, ".mPrefixLines", NULL),
    paste0(.indent, .newVar, " = 0"),
    paste0(.indent, "if ", .cond),
    paste0(.indent, "  ", .newVar, " = 1"),
    paste0(.indent, "end")
  )
  rxode2::rxAssignControlValue(ui, ".mPrefixLines", .prefixLines)
  rxode2::rxAssignControlValue(
    ui,
    ".mLogicalDf",
    rbind(.df, data.frame(cond = .cond, nm = .newVar))
  )
  .newVar
}

.rxToMonolixHandleIfExpressions <- function(x, ui) {
  .cond <- .rxToMonolixCondition(x[[2]], ui)
  .ret <- paste0(.rxToMonolixFlushPrefixLines(ui),
                 .rxToMonolixGetIndent(ui), "if ", .cond, "\n")
  .rxToMonolixIndent(ui)
  .ret <- paste0(.ret, .rxToMonolix(x[[3]], ui=ui))
  x <- x[-c(1:3)]
  if (length(x) == 1) x <- x[[1]]
  while(identical(x[[1]], quote(`if`))) {
    .cond <- .rxToMonolixCondition(x[[2]], ui)
    if (!identical(rxode2::rxGetControl(ui, ".mPrefixLines", NULL), NULL)) {
      stop(
        "an `else if` condition cannot use a logical expression as a ",
        "number in Monolix; prune with `monolixControl(prune=TRUE)`",
        call. = FALSE
      )
    }
    .ret <- paste0(.ret, "\n",
                   .rxToMonolixGetIndent(ui, FALSE), "elseif ", .cond, "\n")
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

#' Mu-referenced theta to Monolix variable name map
#'
#' Monolix variable names cannot contain `.`; `.rxToMonolix()` spells
#' them with `__` in the equations, so every other Monolix section
#' (and the output readers) must use the same spelling.
#'
#' @param ui rxode2 ui
#' @return named character vector; names are the thetas, values are
#'   the Monolix variable names
#' @author Matthew L. Fidler
#' @noRd
.monolixMuRef <- function(ui) {
  .split <- ui$getSplitMuModel
  .muRef <- c(.split$pureMuRef, .split$taintMuRef)
  setNames(gsub("[.]", "__", .muRef), names(.muRef))
}

#' Transformations of the Monolix individual parameters
#'
#' Like `ui$muRefCurEval`, but a "tainted" mu-referenced parameter (like
#' `tcl` in `cl <- exp(tcl + eta.cl) * (CRCL/100)^cl.crcl`) is linear:
#' Monolix only gets `rx__tcl <- tcl` and the model keeps the `exp()`,
#' so a log-normal distribution would apply the `exp()` twice.
#'
#' @param ui rxode2 ui
#' @return `muRefCurEval` data frame
#' @noRd
.monolixMuRefCurEval <- function(ui) {
  .ret <- ui$muRefCurEval
  .taint <- names(ui$getSplitMuModel$taintMuRef)
  if (length(.taint) == 0L) {
    return(.ret)
  }
  .mrt <- ui$muRefTable
  .eta <- .mrt$eta[.mrt$theta %in% .taint]
  .w <- .ret$parameter %in% c(.taint, .eta)
  .ret$curEval[.w] <- ""
  .ret$low[.w] <- NA_real_
  .ret$hi[.w] <- NA_real_
  .ret
}

.rxToMonolix <- function(x, ui) {
  ui <- rxode2::rxUiDecompress(ui)
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
    if (length(x) == 3L &&
          (identical(x[[1]], quote(`Rx_pow_di`)) || identical(x[[1]], quote(`Rx_pow`)))) {
      # rxode2's normalized powers
      return(paste0("(", .rxToMonolix(x[[2]], ui=ui), ")^(", .rxToMonolix(x[[3]], ui=ui), ")"))
    }
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
    } else if (.rxIsLogicalOperator(x[[1]]) || identical(x[[1]], quote(`!`))) {
      # a logical expression used as a number
      return(.rxToMonolixLogicalIndicator(x, ui))
    } else if (identical(x[[1]], quote(`**`)) ) {
      return(paste(.rxToMonolix(x[[2]], ui=ui), "^", .rxToMonolix(x[[3]], ui=ui)))
    } else if (.rxIsAssignmentOperator(x[[1]])) {
      .prop <- as.character(x[[2]])[1]
      if (is.call(x[[2]]) &&
            any(.prop == c("alag", "lag", "F", "f", "rate", "dur"))) {
        return(.rxToMonolixCmtProp(x, .prop, ui))
      }
      .var <- .rxToMonolix(x[[2]], ui=ui)
      .val <- .rxToMonolix(x[[3]], ui = ui)
      return(paste0(.rxToMonolixFlushPrefixLines(ui),
                    paste(.rxToMonolixGetIndent(ui), .var, "=", .val)))
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
            stop(sprintf("'%s' only accepts 1 argument", .x1), call. = FALSE)
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
      .frames <- seq_len(sys.nframe())
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

#' Replace variables in a model expression
#'
#' @param x model expression
#' @param map named list of new names (symbols), by old name
#' @return the expression with the variables replaced
#' @noRd
.monolixSubst <- function(x, map) {
  if (length(map) == 0L) {
    return(x)
  }
  if (is.name(x)) {
    .n <- as.character(x)
    if (.n %in% names(map)) {
      return(map[[.n]])
    }
    return(x)
  }
  if (is.call(x)) {
    for (.i in seq_along(x)[-1]) {
      if (!is.null(x[[.i]])) x[[.i]] <- .monolixSubst(x[[.i]], map)
    }
  }
  x
}

#' Variables assigned (with `<-` or `=`) in a model expression
#'
#' @param x model expression
#' @return character vector of the assigned variables
#' @noRd
.monolixAssigned <- function(x) {
  if (!is.call(x)) {
    return(character(0))
  }
  .ret <- unlist(lapply(as.list(x)[-1], .monolixAssigned))
  if ((identical(x[[1]], quote(`<-`)) || identical(x[[1]], quote(`=`))) &&
        is.name(x[[2]])) {
    .ret <- c(as.character(x[[2]]), .ret)
  }
  unique(.ret)
}

#' Give every reassignment a new variable for Monolix
#'
#' Monolix variables are assigned once ("Conflicting variable
#' definition"), so `cl <- cl * 1.2` after `cl` is defined (or is an
#' individual parameter) becomes `cl_rx1 <- cl * 1.2`, and the rest of
#' the model uses `cl_rx1`.  `if`/`else` statements that change a
#' defined variable are pruned first (see `.bblIfReassigns()`).
#'
#' @param lstExpr model lines
#' @param defined variables already defined (individual parameters,
#'   regressors)
#' @param keep variables never renamed (the endpoints)
#' @return model lines
#' @noRd
.monolixSsa <- function(lstExpr, defined, keep = character(0)) {
  .map <- list()
  .count <- integer(0)
  for (.i in seq_along(lstExpr)) {
    .e <- lstExpr[[.i]]
    if (identical(.e, quote(`_drop`))) next
    if ((identical(.e[[1]], quote(`<-`)) || identical(.e[[1]], quote(`=`))) &&
          is.name(.e[[2]])) {
      .v <- as.character(.e[[2]])
      .rhs <- .monolixSubst(.e[[3]], .map)
      if (.v %in% defined && !(.v %in% keep)) {
        .count[.v] <- ifelse(is.na(.count[.v]), 1L, .count[.v] + 1L)
        .new <- as.name(paste0(.v, "_rx", .count[.v]))
        .map[[.v]] <- .new
        .e <- as.call(list(.e[[1]], .new, .rhs))
      } else {
        .e[[3]] <- .rhs
        defined <- c(defined, .v)
      }
    } else {
      .e <- .monolixSubst(.e, .map)
      defined <- c(defined, .monolixAssigned(.e))
    }
    lstExpr[[.i]] <- .e
  }
  lstExpr
}
