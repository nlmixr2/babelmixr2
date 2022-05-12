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
  "log10" = c("log10(", ""),
  "log2" = c("log(", ")*1.442695040888963387005"),
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

# "logit" in monolix isn't as flexible as nlmixr
# "probit" in monolix isn't as flexible as nlmixr

.monolixTlag <- c()
.monolixP <- c()
#' Get the monolix administration info
#'
#'
#' @param ui rxode2 user interface
#' @return internal adm dataset that comes from `nlmixr2extra::nlmixrDataToMonolix`
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
        } else if (regexpr("(?:-)?(?:(?:0|(?:[1-9][0-9]*))|(?:(?:[0-9]+\\.[0-9]*)|(?:[0-9]*\\.[0-9]+))(?:(?:[Ee](?:[+\\-])?[0-9]+))?|[0-9]+[Ee](?:[\\-+])?[0-9]+)",
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
      assignInMyNamespace(".monolixTlag", c())
      assignInMyNamespace(".monolixP", c())
      .x2 <- x[-1]
      .ret <- paste(lapply(.x2, function(x) {
        .rxToMonolix(x, ui=ui)
      }), collapse = "\n")
      return(.ret)
    } else if (identical(x[[1]], quote(`*`)) ||
      identical(x[[1]], quote(`^`)) ||
      identical(x[[1]], quote(`+`)) ||
      identical(x[[1]], quote(`-`)) ||
      identical(x[[1]], quote(`/`))) {
      if (length(x) == 3) {
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
      } else {
        ## Unary Operators
        return(paste(
          as.character(x[[1]]),
          .rxToMonolix(x[[2]], ui=ui)
        ))
      }
    } else if (identical(x[[1]], quote(`if`))) {
      .ret <- paste0("if ", .rxToMonolix(x[[2]], ui=ui), "\n",
                     "  ", .rxToMonolix(x[[3]]), ui=ui)
      x <- x[-c(1:3)]
      if (length(x) == 1) x <- x[[1]]
      while(identical(x[[1]], quote(`if`))) {
        .ret <- paste0(.ret, "\nelseif ", .rxToMonolix(x[[2]], ui=ui), "\n",
                       "  ", .rxToMonolix(x[[3]], ui=ui))
        x <- x[-c(1:3)]
        if (length(x) == 1) x <- x[[1]]
      }
      if (is.null(x)) {
        .ret <- paste0(.ret, "\nend\n")
      }  else {
        .ret <- paste0(.ret, "\nelse \n",
                       "  ", .rxToMonolix(x, ui=ui),
                       "\nend\n")
      }
      return(.ret)
    } else if (identical(x[[1]], quote(`==`)) ||
                 identical(x[[1]], quote(`>`)) ||
                 identical(x[[1]], quote(`<`)) ||
                 identical(x[[1]], quote(`<=`)) ||
                 identical(x[[1]], quote(`>=`)) ||
                 identical(x[[1]], quote(`!=`)) ||
                 identical(x[[1]], quote(`&&`)) ||
                 identical(x[[1]], quote(`||`)) ||
                 identical(x[[1]], quote(`|`)) ||
                 identical(x[[1]], quote(`&`))) {
        ## Use "preferred" monolix syntax
      return(paste0(.rxToMonolix(x[[2]], ui=ui), as.character(x[[1]]), .rxToMonolix(x[[3]], ui=ui)))
    } else if (identical(x[[1]], quote(`!`)) ) {
      ## Use "preferred" monolix syntax
      return(paste0("~", .rxToMonolix(x[[2]], ui=ui)))
    } else if (identical(x[[1]], quote(`**`)) ) {
      return(paste(.rxToMonolix(x[[2]], ui=ui), "^", .rxToMonolix(x[[3]], ui=ui)))
    } else if (identical(x[[1]], quote(`=`)) ||
      identical(x[[1]], quote(`<-`)) ||
        identical(x[[1]], quote(`~`))) {
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
        return(paste0(";", as.character(x[[2]])[1], " defined in PK section"))
      }
      .var <- .rxToMonolix(x[[2]], ui=ui)
      return(paste(.var, "=", .rxToMonolix(x[[3]], ui=ui)))
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
    }  else {
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

##' Get the last babelmix nlmixr project exported
##'
##' @export
nlmixrMonolixLastProject <- function(){
  return(.nlmixrMonolixLastProject)
}
##' Converts nlmixr to monolix
##'
##' @param uif User interface function or parsed function from nlmixr
##' @param data Data for nlmixr
##' @param control Monolix control
##' @return Nothing; Writes files for monolix to run
##' @author Matthew Fidler
##' @export
nlmixrToMonolix <- function(uif, data, control=monolixControl()){
  name <- as.character(substitute(uif))
  if (!inherits(uif, "rxUi")) {
    uif <- nlmixr2::nlmixr2(uif)
  }
  data <- as.data.frame(data)
  if(!any(tolower(names(data)) == "dv")) stop("need dv in data", call.=FALSE)
  .wid <- which(tolower(names(data)) == "id")
  names(data)[.wid] <- "id"
  .waddl <- which(tolower(names(data)) == "addl")
  if (length(.waddl) == 1) {
    data[, .waddl] <- ifelse(is.na(data[, .waddl]), 0, is.na(data[, .waddl]))
  }
  .lst <- monolixModelTxt(uif, data, control=control, name=name)
  .rds <- paste0(.lst$file, ".rds")
  if (file.exists(.rds)){
    .lst2 <- readRDS(.rds)
    if (.lst$digest == .lst2$lst$digest) {
      .nlmixr <- paste0(.lst$file, "-nlmixr.rds")
      if (file.exists(.nlmixr)) {
        return(readRDS(.nlmixr))
      }
      .populationParameters <- file.path(getwd(), .lst$file, "populationParameters.txt")
      if (file.exists(.populationParameters)) {
        .populationParameters <- data.table::fread(.populationParameters)
        .individualParameters <- data.table::fread(file.path(getwd(), .lst$file, "IndividualParameters", "estimatedRandomEffects.txt"))
        ## Can be SA or Lin
        .covarianceEstimation <- file.path(getwd(), .lst$file, "FisherInformation", "covarianceEstimatesLin.txt")
        if (file.exists(.covarianceEstimation)) {
          .covarianceEstimation <- data.table::fread(.covarianceEstimation)
        } else {
          .covarianceEstimation <- file.path(getwd(), .lst$file, "FisherInformation", "covarianceEstimatesSA.txt")
          if (file.exists(.covarianceEstimation)) {
            .covarianceEstimation <- data.table::fread(.covarianceEstimation)
          }
        }
        ## Monolix final parameters are not on the original nlmixr scale
        .f1 <- setNames(as.data.frame(.populationParameters)[,c("parameter", "value")], c("typical", "thetaEstF"))
        .f2 <- setNames(as.data.frame(.populationParameters)[,c("parameter", "value")], c("sd", "sdEstF"))
        .final <- merge(merge(.lst2$lst$df, .f1, all.x=TRUE), .f2, all.x=TRUE)
        .final$varF <- .final$sdEstF^2
        .final$thetaF <- sapply(seq_along(.final$theta), function(i){
          .trans <- .final[i, "trans"]
          if (.trans == "logNormal") {
            return(log(.final$thetaEstF[i]))
          } else if (.trans == "logitNormal") {
            return(RODE::logit(.final$thetaEstF[i], .final$low[i], final$hi[i]))
          } else if (.trans == "probitNormal") {
            return(RODE::probit(.final$thetaEstF[i]))
          }
        })
        ## Now that they are merged with the original translation
        ## table, transform them to be on the same scale
        .ini <- as.data.frame(uif$ini)
        .ini$est <- sapply(seq_along(.ini$ntheta), function(i) {
          .name <- .ini$name[i]
          .w <- which(.final$theta == .name)
          if (length(.w) == 1) {
            return(.final$thetaF[.w])
          }
          .w <- which(.final$eta == .name)
          if (length(.w) == 1) {
            return(.final$varF[.w])
          }
          .namem <- eval(parse(text=paste0("rxToMonolix(", .name, ", ui=uif)")))
          .w <- which(.namem == .populationParameters$parameter)
          if (length(.w) == 1){
            return(.populationParameters$value[.w])
          }
          message("cannot find ", .name)
          return(NA_real_)
        })
        class(.ini) <- c("nlmixrBounds", "data.frame")
        .uif <- uif
        .uif$ini <- .ini

        # ETAS in  eta_parameterName_SAEM
        .etas <- .final$sd
        .etas <- .etas[!is.na(.etas)]
        .etas <- c("id", gsub("omega_(.*)", "eta_\\1_SAEM", .etas))
        .etas <- as.data.frame(.individualParameters)[, .etas]
        names(.etas) <- sapply(names(.etas), function(.n) {
          if (.n == "id") return("id")
          .n <- sub("eta_(.*)_SAEM", "omega_\\1", .n)
          .w <- which(.final$sd == .n)
          if (length(.w) == 1){
            return(.final$eta[.w])
          }
          stop("cannot determine nlmixr ETAs from monolix output")
        })
        .wid <- which(tolower(names(data)) == "id")
        if (length(.wid) == 1) {
          .etas <- merge(data.frame(id=unique(data[, .wid])), .etas, by="id", all.x=TRUE)
          .etas[is.na(.etas)] <- 0
          .etas <- .etas[order(.etas$id), ]
        }
        .etaMat <- as.matrix(.etas[, paste(.ini$name[which(.ini$neta1 == .ini$neta2)])])
        .nid <- length(unique(data[, .wid]))
        if (.nid != dim(.etaMat)[1]) {
          message("possibly corrupted run, remove ", .rds, " and associated directory to possibly fix")
          stop("number of IDs (", .nid, ") not equal to number of rows of eta matrix (", dim(.etaMat)[1])
        }

        .ctl <- list()
        .ctl$covMethod <- ""
        .ctl$maxOuterIterations <- 0
        .ctl$maxInnerIterations <- 0
        .ctl$scaleTo <- 0
        .ctl$addProp <- control$addProp
        .ctl$method <- "liblsoda"
        .ctl <- do.call(nlmixr2::foceiControl, .ctl)
        if (any(names(.ctl) == "singleOde")){
          if (control$singleOde) {
            .mod <- uif$focei.rx1
            .pars <- NULL
          } else {
            .mod <- uif$rxode.pred
            .pars <- uif$theta.pars
          }
        } else {
          .mod <- uif$rxode.pred
          .pars <- uif$theta.pars
        }
        .env <- new.env(parent = emptyenv())
        .env$method <- "Monolix"
        .env$extra <- " (babelmixr2)"
        .data2 <- data
        .allCovs <- tolower(.uif$all.covs)
        names(.data2) <- sapply(names(.data2), function(x){
          .w <- which(tolower(x) == .allCovs)
          if (length(.w) == 1) return(.uif$all.covs[.w])
          if (tolower(x) == "ytype") return("ytypeOld")
          return(x)
        })
        .nlmixrFit <- foceiFit(
          data = .data2, inits = .uif$focei.inits, PKpars = .pars,
          model = .mod, pred = function() {
            return(nlmixr_pred)
          }, err = .uif$error,
          lower = .uif$focei.lower,
          upper = .uif$focei.upper,
          thetaNames = .uif$focei.names,
          etaNames = .uif$eta.names,
          etaMat = .etaMat,
          env = .env,
          fixed = uif$focei.fixed,
          control = .ctl
        )
        saveRDS(.nlmixrFit, .nlmixr)
        return(.nlmixrFit)
      }
      message("the monolix model is current with the nlmixr model; remove ", .rds," to regenerate")
      .status <- getOption("babelmixr2.monolix.status", "")
      if (.status != "") {
        system(.status)
        # nvs specific errors
        .files <- list.files(pattern=paste0(name, "[.]e"))
        if (length(.files) > 0) {
          for (.f in .files){
            message("file:", .f)
            message("================================================================================\n")
            message(paste(suppressWarnings(readLines(.f)), collapse="\n"))
          }
        }
      }
      return(invisible())
    } else {
      message("regenerating monolix files because data or model have changed")
    }
  }
  .mlx <- paste0(.lst$file, ".mlxtran")
  assignInMyNamespace(".nlmixrMonolixLastProject", .mlx)
  ## cli::alert_info("writing mlxtran file to {.mlx}")
  message("writing mlxtran file to ", .mlx)
  writeLines(.lst$mlxtran, .mlx)
  .txt <- paste0(.lst$file, ".txt")
  ## cli::alert_info("writing monolix model txt file to {.txt}")
  message(paste0("writing monolix model txt file to ", .txt))
  writeLines(.lst$txt, .txt)
  .data <- paste0(.lst$data.md5, ".csv")
  if (!file.exists(.data)) {
    message(paste0("writing data to ", .data))
    ## cli::alert_info("writing data to {.data}")
    data.table::fwrite(data, .data)
  }
  .rds <- paste0(.lst$file, ".rds")
  message(paste0("writing monolix translation information to ", .rds))
  saveRDS(list(uif=uif, data=data, control=control, lst=.lst), .rds)
  if (control$runCommand != "") {
    system(sprintf(control$runCommand, .mlx))
  }
  return(invisible())
}

## lixoftConnectors::newProject(data = list(dataFile = "/path/to/data/file.txt",
##                                          headerTypes = c("IGNORE","OBSERVATION"),
##                                          observationTypes = "continuous"),
##                              modelFile = "/path/to/model/file.txt")

## Example with warfarin_data.txt from demos and oral1_1cpt_kaVCl.txt from libraries in the current directory
## data = list(dataFile= "./warfarin_data.txt",
##             headerTypes =c("id", "time", "amount", "observation", "obsid", "contcov", "catcov", "ignore"),
##             observationTypes = list(y1 = "continuous", y2 = "continuous" ),
##             mapping = list("1" = "y1", "2" = "y2"))
## modelFile <- './oral1_1cpt_kaVCl.txt'
## newProject(modelFile = modelFile, data = data)

