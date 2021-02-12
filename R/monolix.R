
##
## observationTypes (list): A list giving the type of each observation present in the data file. If there is only one y-type, the corresponding observation name can be omitted.
## The possible observation types are "continuous", "discrete", and "event".
##
## nbSSDoses [optional](int): Number of doses (if there is a SS column).


##' Monolix Controller for nlmixr
##'
##' @param nbSSDoses Number of steady state doses (default 7)
##' @param stiff boolean for using the stiff ODE solver
##' @param exploratoryautostop logical to turn on or off exploratory
##'   phase auto-stop of SAEM (default 250)
##' @param exploratoryiterations Number of iterations for exploratory
##'   phase (default 250)
##' @param exploratoryinterval Minimum number of interation in the
##'   exploratory phase (default 200)
##' @param simulatedannealingiterations Number of burn in iterations
##' @param runCommand is a shell command to run monolix; You can
##'   specfy the default by
##'   \code{options("babelmixr.monolix"="runMonolix \%s")} where the \code{"\%s"}
##'   represents the monolix project file.
##' @inheritParams nlmixr::foceiControl
##' @return A monolix control object
##' @author Matthew Fidler
##' @export
##' @importFrom nlmixr nlmixr
##' @importFrom methods is
##' @importFrom stats na.omit setNames
##' @importFrom utils assignInMyNamespace
monolixControl <- function(nbSSDoses=7,
                           stiff=FALSE,
                           addProp = c("combined2", "combined1"),
                           exploratoryautostop=FALSE,
                           smoothingautostop=FALSE,
                           simulatedannealing=TRUE,
                           burniniterations=5,
                           smoothingiterations=200,
                           exploratoryiterations=250,
                           simulatedannealingiterations=250,
                           exploratoryinterval=200,
                           exploratoryalpha=0.0,
                           omegatau=0.95,
                           errormodeltau=0.95,
                           optimizationiterations=20,
                           optimizationtolerance=0.0001,
                           variability=c("none", "firstStage", "decreasing"),
                           runCommand=getOption("babelmixr.monolix", "")) {

  checkmate::assertLogical(stiff, max.len=1)
  checkmate::assertLogical(exploratoryautostop, max.len=1)
  checkmate::assertLogical(smoothingautostop, max.len=1)
  checkmate::assertLogical(simulatedannealing, max.len=1)

  checkmate::assertIntegerish(burniniterations, max.len=1, lower=1)
  checkmate::assertIntegerish(exploratoryiterations, max.len=1, lower=1)
  checkmate::assertIntegerish(simulatedannealingiterations, max.len=1, lower=1)
  checkmate::assertIntegerish(nbSSDoses, lower=1, max.len=1)
  checkmate::assertIntegerish(exploratoryiterations, max.len=1, lower=1)
  checkmate::assertIntegerish(exploratoryinterval, max.len=1, lower=1)
  checkmate::assertIntegerish(smoothingiterations, max.len=1, lower=1)
  checkmate::assertIntegerish(optimizationiterations, max.len=1, lower=1)

  checkmate::assertNumeric(exploratoryalpha, lower=0.0, upper=1.0)
  checkmate::assertNumeric(omegatau, lower=0.0, upper=1.0)
  checkmate::assertNumeric(errormodeltau, lower=0.0, upper=1.0)
  checkmate::assertNumeric(optimizationtolerance, lower=0.0)
  if (optimizationtolerance == 0) stop("'optimizationtolerance' has to be above zero")

  if (runCommand != "") checkmate::assertCharacter(runCommand, pattern="%s", min.len=1, max.len=1)

  .ret <- list(nbSSDoses=as.integer(nbSSDoses), stiff=stiff,
               exploratoryautostop=exploratoryautostop,
               smoothingautostop=smoothingautostop,
               addProp=match.arg(addProp),
               burniniterations=burniniterations,
               simulatedannealingiterations=simulatedannealingiterations,
               exploratoryinterval=exploratoryinterval,
               smoothingiterations=smoothingiterations,
               exploratoryalpha=exploratoryalpha,
               omegatau=omegatau,
               errormodeltau=errormodeltau,
               exploratoryiterations=exploratoryiterations,
               optimizationiterations=optimizationiterations,
               optimizationtolerance=optimizationtolerance,
               variability=match.arg(variability),
               runCommand=runCommand
               )
  class(.ret) <- "monolixControl"
  .ret
}

.monolixErrs <- c()

##' This constructs the data->monolix header mapping & regressors
##'
##' @param data Input dataset
##' @param uif  Parsed nlmixr user interface function
##' @return list with (header=monolix header specification, regressor=monolix regressor specification)
##' @author Matthew Fidler
monolixMapData <- function(data, uif) {
  ## Monolix header types "id", "time", "observation", "amount",
  ## "contcov", "catcov", "occ", "evid", "mdv", "obsid", "cens",
  ## "limit", "regressor","admid", "rate", "tinf", "ss", "ii", "addl",
  ## "date", "ignore".
  ##' Generate the monolix input line
  regressors <- RxODE::rxModelVars(uif$rxode)$params
  covariates <- uif$saem.all.covs
  .env <- environment()
  .env$.regressor <- c(paste0("input={", paste(regressors, collapse=", "), "}"))
  .headers <- sapply(names(data), function(x) {
    if (tolower(x) == "id") return("id")
    if (tolower(x) == "dv") return("observation")
    if (tolower(x) == "amt") return("amount")
    if (tolower(x) == "evid") return("evid")
    if (tolower(x) == "mdv") return("mdv")
    if (tolower(x) == "cmt") return("obsid")
    if (tolower(x) == "cens") return("cens")
    if (tolower(x) == "limit") return("limit")
    if (tolower(x) == "rate") return("rate")
    if (tolower(x) == "dur") return("tinf")
    if (tolower(x) == "ss") return("ss")
    if (tolower(x) == "ii") return("ii")
    if (tolower(x) == "addl") return("addl")
    if (tolower(x) == "occ") return("occ")
    if (tolower(x) == "time") return("time")
    .w <- which(tolower(x) == tolower(regressors))
    if (length(.w) == 1) {
      .env$.regressor <- c(.env$.regressor,
                           paste0(regressors[.w], " = {use=regressor}"))
      return("regressor")
    }
    if (any(x == covariates)) {
      if (inherits(x, "factor") || inherits(x, "character")) {
        return("catcov")
      } else {
        return("contcov")
      }
    }
    return("ignore")
  })
  return(list(headerType=.headers, regressors=paste(.env$.regressor, collapse="\n")))
}

##' Monolix estimation routine
##'
##' @param env nlmixr environment
##'
##' @param ... Other parameters
##'
nlmixrEst.monolix <- function(env, ...){
  with(env, {
    ## obj$env$.curTv needs to be reset
    if (length(uif$noMuEtas) > 0) {
      stop(sprintf("Cannot run Monolix since some of the parameters are not mu-referenced (%s)", paste(uif$noMuEtas, collapse = ", ")),
           call.=FALSE)
    }
  })
}

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

.rxToMonolix <- function(x) {
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
      return(paste0("(", .rxToMonolix(x[[2]]), ")"))
    } else if (identical(x[[1]], quote(`{`))) {
      assignInMyNamespace(".monolixTlag", c())
      assignInMyNamespace(".monolixP", c())
      .x2 <- x[-1]
      .ret <- paste(lapply(.x2, function(x) {
        .rxToMonolix(x)
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
              .state <- .rxToMonolix(.x3[[2]])
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
              .rxToMonolix(.x2),
              as.character(x[[1]]),
              .rxToMonolix(.x3)
            )
          }
        } else {
          .ret <- paste0(
            .rxToMonolix(x[[2]]),
            as.character(x[[1]]),
            .rxToMonolix(x[[3]])
          )
        }
        return(.ret)
      } else {
        ## Unary Operators
        return(paste(
          as.character(x[[1]]),
          .rxToMonolix(x[[2]])
        ))
      }
    } else if (identical(x[[1]], quote(`if`))) {
      .ret <- paste0("if ", .rxToMonolix(x[[2]]), "\n",
                     "  ", .rxToMonolix(x[[3]]))
      x <- x[-c(1:3)]
      if (length(x) == 1) x <- x[[1]]
      while(identical(x[[1]], quote(`if`))) {
        .ret <- paste0(.ret, "\nelseif ", .rxToMonolix(x[[2]]), "\n",
                       "  ", .rxToMonolix(x[[3]]));
        x <- x[-c(1:3)]
        if (length(x) == 1) x <- x[[1]]
      }
      if (is.null(x)) {
        .ret <- paste0(.ret, "\nend\n")
      }  else {
        .ret <- paste0(.ret, "\nelse \n",
                       "  ", .rxToMonolix(x),
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
      return(paste0(.rxToMonolix(x[[2]]), as.character(x[[1]]), .rxToMonolix(x[[3]])))
    } else if (identical(x[[1]], quote(`!`)) ) {
      ## Use "preferred" monolix syntax
      return(paste0("~", .rxToMonolix(x[[2]])))
    } else if (identical(x[[1]], quote(`**`)) ) {
      return(paste(.rxToMonolix(x[[2]]), "^", .rxToMonolix(x[[3]])))
    } else if (identical(x[[1]], quote(`=`)) ||
      identical(x[[1]], quote(`<-`)) ||
        identical(x[[1]], quote(`~`))) {
      if (any(as.character(x[[2]])[1] == c("alag", "lag", "F", "f", "rate", "dur"))) {
        if (any(as.character(x[[2]])[1] == c("alag", "lag"))) {
          .state <- as.character(x[[2]][[2]])
          .extra <- .rxToMonolix(x[[3]])
          assignInMyNamespace(".monolixTlag", c(.monolixTlag, setNames(.extra, .state)))
        }
        if (any(as.character(x[[2]])[1] == c("F", "f"))) {
          .state <- as.character(x[[2]][[2]])
          .extra <- .rxToMonolix(x[[3]])
          assignInMyNamespace(".monolixP", c(.monolixP, setNames(.extra, .state)))
        }
        return(paste0(";", as.character(x[[2]])[1], " defined in PK section"))
      }
      .var <- .rxToMonolix(x[[2]])
      return(paste(.var, "=", .rxToMonolix(x[[3]])))
    } else if (identical(x[[1]], quote(`[`))) {
      .type <- toupper(as.character(x[[2]]))
      if (any(.type == c("THETA", "ETA"))) {
        stop("'THETA'/'ETA' not supported by monolix", call.=FALSE);
      }
    } else if (identical(x[[1]], quote(`log1pmx`))) {
      if (length(x == 2)) {
        .a <- .rxToMonolix(x[[2]])
        return(paste0("(log(1+", .a, ")-(", .a, "))"))
      } else {
        stop("'log1pmx' only takes 1 argument", call. = FALSE)
      }
    } else if ((identical(x[[1]], quote(`pnorm`))) |
                 (identical(x[[1]], quote(`normcdf`))) |
                 (identical(x[[1]], quote(`phi`)))) {
      if (length(x) == 4) {
        .q <- .rxToMonolix(x[[2]])
        .mean <- .rxToMonolix(x[[3]])
        .sd <- .rxToMonolix(x[[4]])
        return(paste0("normcdf(((", .q, ")-(", .mean, "))/(", .sd, "))"))
      } else if (length(x) == 3) {
        .q <- .rxToMonolix(x[[2]])
        .mean <- .rxToMonolix(x[[3]])
        return(paste0("normcdf(((", .q, ")-(", .mean, ")))"))
      } else if (length(x) == 2) {
        .q <- .rxToMonolix(x[[2]])
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
              .xc[1], .rxToMonolix(x[[2]]),
              .xc[2])
            return(.ret)
          } else {
            stop(sprintf("'%s' only acceps 1 argument", .x1), call. = FALSE)
          }
        }
      }
      .ret0 <- c(list(as.character(x[[1]])), lapply(x[-1], .rxToMonolix))
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

##' Convert RxODE syntax to monolix syntax
##'
##' @param x Expression
##' @return Monolix syntax
##' @author Matthew Fidler
##' @export
rxToMonolix <- function(x) {
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
            return(.rxToMonolix(.val2))
          } else if (inherits(.val2, "numeric") || inherits(.val2, "integer")) {
            return(sprintf("%s", .val2))
          }
        }
      }
    }
    return(.rxToMonolix(x))
  }
  return(.rxToMonolix(eval(parse(text = paste0("quote({", x, "})")))))
}

## FIXME these can be saved and retreived if rxToMonolix(.v) is run first
monolixTlag <- function(states) {
  sapply(states, function(state) {
    .cur <- .monolixTlag[state]
    if (length(.cur) == 0) return("0.0")
    if (is.na(.cur)) return("0.0")
    .cur
  })
}

monolixP <- function(states){
  sapply(states, function(state) {
    .cur <- .monolixP[state]
    if (length(.cur) == 0) return("1.0")
    if (is.na(.cur)) return("1.0")
    .cur
  })
}

.monolixGetErr <- list()

monolixGetErr0 <- function(cond, type, uif, control) {
  .ini <- as.data.frame(uif$ini)
  .ini <- .ini[which(.ini$condition == cond), ]
  if (type == 1) { # Constant
    .tmp <- as.character(.ini$name)
    .tmp <- eval(parse(text=paste0("rxToMonolix(", .tmp, ")")))
    assignInMyNamespace(".monolixErrs", c(.monolixErrs, .tmp))
    .lst <- .monolixGetErr
    .lst[[length(.lst) + 1]] <- data.frame(name=.ini$name, errName=.tmp)
    assignInMyNamespace(".monolixGetErr", .lst)
    return(paste0("constant(", .tmp, ")"))
  } else if (type == 2) { # Proportional
    .tmp <- as.character(.ini$name)
    .tmp <- eval(parse(text=paste0("rxToMonolix(", .tmp, ")")))
    .lst <- .monolixGetErr
    .lst[[length(.lst) + 1]] <- data.frame(name=.ini$name, errName=.tmp)
    assignInMyNamespace(".monolixGetErr", .lst)
    assignInMyNamespace(".monolixErrs", c(.monolixErrs, .tmp))
    return(paste0("proportional(", .tmp, ")"))
  } else if (type == 3) { # additive + proportional
    .add <- .ini[.ini$err == "add", "name"]
    .add <- eval(parse(text=paste0("rxToMonolix(", .add, ")")))
    .lst <- .monolixGetErr
    .lst[[length(.lst) + 1]] <- data.frame(name=.ini[.ini$err == "add", "name"], errName=.add)
    .prop <- .ini[.ini$err == "prop", "name"]
    .prop <- eval(parse(text=paste0("rxToMonolix(", .prop, ")")))
    .lst[[length(.lst) + 1]] <- data.frame(name=.ini[.ini$err == "prop", "name"], errName=.prop)
    assignInMyNamespace(".monolixErrs", c(.monolixErrs, .add, .prop))
    assignInMyNamespace(".monolixGetErr", .lst)
    return(paste0(control$addProp, "(", .add, ",", .prop, ")"))
  } else if (type == 4) { # additive + power
    stop("distribution not supported in monolix<->nlmixr")
  } else if (type == 5) { # pow
    stop("distribution not supported in monolix<->nlmixr")
  } else if (type >= 6) { ## + lambda
    stop("distribution not supported in monolix")
  }
}

monolixGetErr <- function(resMod, uif, control) {
  assignInMyNamespace(".monolixErrs", c())
  assignInMyNamespace(".monolixGetErr", list())
  .idx <- seq_along(resMod)
  sapply(.idx, function(.i){
    monolixGetErr0(names(resMod)[.i], setNames(resMod[.i], NULL), uif, control)
  })
}

.toMonolixDef <- list()
.toMonolixDefinition <- function(x, mu.ref) {
  if (is.call(x)) {
    if (identical(x[[1]], quote(`{`))) {
      assignInMyNamespace(".toMonolixDef", list())
      .x2 <- x[-1]
      .ret <- list(paste0("DEFINITION:\n", paste(lapply(.x2, function(x) {
        .toMonolixDefinition(x, mu.ref)
      }), collapse = "\n")), do.call(rbind, .toMonolixDef))
      return(.ret)
    } else if (identical(x[[1]], quote(`=`))) {
      .var <- as.character(x[[2]])
      .start <- paste0(.var, " = ")
      .x3 <- x[[3]]
      .muRef <- unlist(mu.ref)
      .eta <- names(.muRef)
      .muRef <- setNames(.muRef, NULL)
      if (length(.x3) >= 2){
        if (identical(.x3[[1]], quote(`exp`))) {
          .theta <- as.character(.x3[[2]])
          .w <- which(.muRef == .theta)
          .start <- paste0(.start, "{distribution=logNormal, typical=", .var, "_pop, ")
          if (length(.w) == 1) {
            .md <- .toMonolixDef
            .md[[length(.md) + 1]] <- data.frame(var=.var, theta=.theta, typical=paste0(.var, "_pop"),
                                                 eta=.eta[.w], sd=paste0("omega_", .var),
                                                 trans="logNormal", low=0, hi=1, stringsAsFactors = FALSE)
            assignInMyNamespace(".toMonolixDef", .md)
            .start <- paste0(.start, "sd=omega_", .var, "}")
            return(.start)
          } else {
            .md <- .toMonolixDef
            .md[[length(.md) + 1]] <- data.frame(var=.var, theta=.theta, typical=paste0(.var, "_pop"), eta=NA_character_,
                                                 sd=NA_character_, trans="logNormal", low=0, hi=1, stringsAsFactors = FALSE)
            assignInMyNamespace(".toMonolixDef", .md)
            .start <- paste0(.start, "no-variability}")
            return(.start)
          }
        } else if (identical(.x3[[1]], quote(`expit`))) {
          .theta <- as.character(.x3[[2]])
          .w <- which(.muRef == .theta)
          .start <- paste0(.start, "{distribution=logitNormal, min=",
                           ifelse(length(.x3) >= 3, paste(eval(.x3[[3]])), "0.0"))
          .start <- paste0(.start, ", max=", ifelse(length(.x3) >= 4, paste(eval(.x3[[4]])), "1.0"))
          .start <- paste0(.start, ", typical=", .var, "_pop, ")
          if (length(.w) == 1) {
            .md <- .toMonolixDef
            .md[[length(.md) + 1]] <- data.frame(var=.var, theta=.theta, typical=paste0(.var, "_pop"), eta=.eta[.w], sd=paste0("omega_", .var),
                                                 trans="logitNormal", low=ifelse(length(.x3) >= 3, eval(.x3[[3]]), 0.0),
                                                 hi=ifelse(length(.x3) >= 4, eval(.x3[[4]]), 1.0), stringsAsFactors = FALSE)
            assignInMyNamespace(".toMonolixDef", .md)
            .start <- paste0(.start, "sd=omega_", .var, "}")
            return(.start)
          } else {
            .md <- .toMonolixDef
            .md[[length(.md) + 1]] <- data.frame(var=.var, theta=.theta, typical=paste0(.var, "_pop"), eta=NA_character_, sd=NA_character_,
                                                 trans="logitNormal", low=ifelse(length(.x3) >= 3, eval(.x3[[3]]), 0.0),
                                                 hi=ifelse(length(.x3) >= 4, eval(.x3[[4]]), 1.0),
                                                 stringsAsFactors = FALSE)
            assignInMyNamespace(".toMonolixDef", .md)
            .start <- paste0(.start, "no-variability}")
            return(.start)
          }
          return(.start)
        } else if (identical(.x3[[1]], quote(`probitInv`))) {
          if (length(.x3) != 2) stop("'probitInv' cannot have more than 1 argument for a monolix model")
          .theta <- as.character(.x3[[2]])
          .w <- which(.muRef == .theta)
          .start <- paste0(.start, "{distribution=probitNormal, typical=", .var, "_pop, ")
          if (length(.w) == 1) {
            .md <- .toMonolixDef
            .md[[length(.md) + 1]] <- data.frame(var=.var, theta=.theta, typical=paste0(.var, "_pop"), eta=.eta[.w], sd=paste0("omega_", .var),
                                                 trans="probitNormal", low=0.0, hi=1.0,
                                                 stringsAsFactors = FALSE)
            assignInMyNamespace(".toMonolixDef", .md)
            .start <- paste0(.start, "sd=omega_", .var, "}")
            return(.start)
          } else {
            .md <- .toMonolixDef
            .md[[length(.md) + 1]] <- data.frame(var=.var, theta=.theta, typical=paste0(.var, "_pop"), eta=NA_character_, sd=NA_character_,
                                                 trans="probitNormal", low=0.0, hi=1.0,
                                                 stringsAsFactors = FALSE)
            assignInMyNamespace(".toMonolixDef", .md)
            .start <- paste0(.start, "no-variability}")
            return(.start)
          }
          return(.start)
        }
      } else {
        .start <- paste0(.start, "{distribution=normal, typical=", .var, "_pop, ")
        .theta <- as.character(.x3)
        .w <- which(.muRef == .theta)
        if (length(.w) == 1) {
          .md <- .toMonolixDef
          .md[[length(.md) + 1]] <- data.frame(var=.var, theta=.theta, typical=paste0(.var, "_pop"), eta=.eta[.w], sd=paste0("omega_", .var),
                                               trans="normal", low=0.0, hi=1.0,
                                               stringsAsFactors = FALSE)
          assignInMyNamespace(".toMonolixDef", .md)
          .start <- paste0(.start, "sd=omega_", .var, "}")
          return(.start)
        } else {
          .md <- .toMonolixDef
          .md[[length(.md) + 1]] <- data.frame(var=.var, theta=.theta, typical=paste0(.var, "_pop"), eta=NA_character_, sd=NA_character_,
                                               trans="normal", low=0.0, hi=1.0,
                                               stringsAsFactors = FALSE)
          assignInMyNamespace(".toMonolixDef", .md)
          .start <- paste0(.start, "no-variability}")
          return(.start)
        }
      }
    }
  }
}

monolixModelParameter <- function(uif) {
  .mu <- uif$nmodel$mu.ref
  .ret <- do.call(rbind, lapply(names(.mu), function(x) {
    data.frame(theta = .mu[[x]], eta = x, stringsAsFactors = FALSE)
  }))
  .ret0 <- merge(data.frame(name=.ret$theta), as.data.frame(uif$ini)[, c("name", "est")])
  names(.ret0)[1] <- "theta"
  return(.ret0)
}

monolixDataContent <- function(lst, uif, data, control=monolixControl()) {
  .headerType <- lst$headerType
  .obs <- ""
  .env <- environment()
  .ret <- paste(sapply(seq_along(.headerType), function(.i) {
    .type <- setNames(.headerType[.i], NULL)
    .col <- names(.headerType)[.i]
    if (.type == "id") {
      return(paste0(.col, " = {use=identifier}"))
    } else if (.type == "amount") {
      return(paste0(.col, " = {use=amount}"))
    } else if (.type == "evid") {
      return(paste0(.col, " = {use=eventidentifier}"))
    } else if (.type == "mdv") {
      return(paste0(.col, " = {use=missingdependentvariable}"))
    } else if (.type == "occ") {
      return(paste0(.col, " = {use=occasion}"))
    } else if (.type == "cens") {
      return(paste0(.col, " = {use=censored}"))
    } else if (.type == "limit") {
      return(paste0(.col, " = {use=limit}"))
    } else if (.type == "catcov") {
      return(paste0(.col, " = {use=covariate, type=categorical}"))
    } else if (.type == "contcov") {
      return(paste0(.col, " = {use=covariate, type=continuous}"))
    } else if (.type == "regressor") {
      return(paste0(.col, " = {use=regressor}"))
    } else if (.type == "admid") {
      return(paste0(.col, " = {use=administration}"))
    } else if (.type == "ii") {
      return(paste0(.col, " = {use=interdoseinterval}"))
    } else if (.type == "addl") {
      return(paste0(.col, " = {use=additionaldose}"))
    } else if (.type == "ss") {
      return(paste0(.col, " = {use=steadystate, nbdoses=", control$nbSSDoses, "}"))
    } else if (.type == "rate") {
      return(paste0(.col, " = {use=rate}"))
    } else if (.type == "tinf") {
      return(paste0(.col, " = {use=infusiontime}"))
    } else if (.type == "obsid") {
      return(paste0(.col, " = {use=observationtype}"))
    } else if (.type == "time") {
      return(paste0(.col, " = {use=time}"))
    } else if (.type == "observation") {
      # With one observation
      .predDf <- uif$nmodel$predDf
      if (length(.predDf$cmt) > 1) {
        return(paste0(.col, " = {use=observation, name={", paste(paste0("y_", .predDf$cmt), collapse=", "), "},yname={",
                      paste(paste0("'", .predDf$cmt, "'"), collapse=", "), "},type={",
                      paste(rep("continuous", length(.predDf$cmt)), collapse=", "), "}}"))
      }
      # With more than one observation
      #return(paste0(.col, " = {use=observation, name={", .col, "}, yname={'2','6'},type={continuous}}"))
      assign(".obs", .col, .env)
      .w <- which(tolower(names(data)) == "evid")
      .wc <- which(tolower(names(data)) == "cmt")
      if (length(.wc) == 1) {
        if (length(.w) == 1) {
          .cmt <- unique(data[data[, .w] == 0, .wc])
          if (length(.cmt) == 1){
            return(paste0(.col, " = {use=observation, name=", .col, ", yname='", .cmt, "', type=continuous}"))
          }
        } else {
          .w <- which(tolower(names(data)) == "mdv")
          if (length(.w) == 1) {
            .cmt <- unique(data[data[, .w] == 0, .wc])
            if (length(.cmt) == 1){
              return(paste0(.col, " = {use=observation, name=", .col, ", yname='", .cmt, "', type=continuous}"))
            }
          }
        }
        stop("more than one compartment for observations, should be a multiple endpoint model", call.=FALSE)
      }
      return(paste0(.col, " = {use=observation, name=", .col, ", type=continuous}"))
    } else {
      return("")
    }
  }), collapse="\n")
  return(list(obs=.obs, ret=.ret))
}

monolixDataFile <- function(lst, uif, data, control=monolixControl()) {
  .lst <- lst
  .cnt <- monolixDataContent(lst, uif, data, control)
  .lst$obs <- .cnt$obs
  .lst$datafile <- paste0("<DATAFILE>\n\n[FILEINFO]\n",
                          "file='", lst$data.md5, ".csv'\n",
                          "delimiter = comma\n",
                          "header = {", paste(names(data), collapse=", "), "}\n\n",
                          "[CONTENT]\n",
                          gsub("\n+", "\n", .cnt$ret)
                          )
  return(.lst)
}


#

monolixModelParameter <- function(.df, .dfError) {
  ## FIXME Covariance estimates?
  #F FIXME Covariate estimates?
  paste0("\n\n<PARAMETER>\n",
         paste(setNames(sapply(seq_along(.df$theta), function(.i){
           .typical <- .df$typical[.i]
           .val <- .df$thetaEst[.i]
           .fixed <- .df$thetaFixed[.i]
           .ret <- paste0(.typical, " = {value=", .val)
           if (.fixed) {
             .ret <- paste0(.ret, ", method=FIXED}")
           } else {
             .ret <- paste0(.ret, ", method=MLE}")
           }
           .sd <- .df$sd[.i]
           if (!is.na(.sd)) {
             .sdEst <- .df$sdEst[.i]
             .fixed <- .df$sdFixed[.i]
             .ret <- paste0(.ret, "\n", .sd, " = {value=", .sdEst)
             if (.fixed) {
               .ret <- paste0(.ret, ", method=FIXED}")
             } else {
               .ret <- paste0(.ret, ", method=MLE}")
             }
           }
           return(.ret)
         }), NULL), collapse="\n"),"\n",
         paste(setNames(sapply(seq_along(.dfError$name), function(.i){
           .errName <- .dfError$errName[.i]
           .val <- .dfError$est[.i]
           .fixed <- .dfError$fix[.i]
           .ret <- paste0(.errName, " = {value=", .val)
           if (.fixed) {
             .ret <- paste0(.ret, ", method=FIXED}")
           } else {
             .ret <- paste0(.ret, ", method=MLE}")
           }
           return(.ret)
         }), NULL), collapse="\n"),
         "\n\n")
}

monolixModelFit <- function(uif, obs) {
  .predDf <- uif$nmodel$predDf
  if (length(.predDf$cmt) > 1) {
    return(paste0("\n\n<FIT>\n",
                  "data = {", paste(paste0("y_", .predDf$cmt), collapse=", "), "}\n",
                  "model = {", paste(paste0(.predDf$cond, "_pred"), collapse=", "), "}\n"))
  }
  paste0("\n\n<FIT>\n",
         "data = ", obs, "\n",
         "model = ", paste0(.predDf$cond, "_pred"), "\n")

}

monolixModelTxt <- function(uif, data, control=monolixControl(), name=NULL) {
  if (inherits(control, "monolixControl")) {
    control <-do.call(monolixControl, control)
  }
  .v <- uif$rxode
  .mv <- RxODE::rxModelVars(.v)
  ## Run this before monolixMapData to populate errors in input={}
  .resMod <- uif$saem.res.mod
  .def <- paste(paste0(names(.resMod), "_pred= {distribution = normal, prediction = ", names(.resMod), ", errorModel=",
                       monolixGetErr(.resMod, uif, control), "}"), collapse="\n")
  ## Now map data
  .map <- monolixMapData(data, uif)
  .mod <- rxToMonolix(.v)

  .lst <- .map
  .lst$data.md5 <- digest::digest(data)
  .lst$file <- paste0(ifelse(missing(name), uif$model.name, name))
  .lst$digest <- digest::digest(list(control, .lst$data.md5, control))
  .lst$txt <- paste0("DESCRIPTION:\n",
         paste0("model translated from babelmixr and nlmixr function ", uif$model.name, " to ", .lst$file, ".txt\n\n"),
         "[LONGITUDINAL]\n",
         .map$regressors,
         ifelse(control$stiff, "\n\nodeType = stiff", ""),
         "\n\nPK:\n; Define compartment(s)\n",
         paste(paste0("compartment(cmt=", seq_along(.mv$state), ", amount=", .mv$state, ")"), collapse="\n"),
         paste("\n\n;Define depot compartment information\n"),
         paste(paste0("depot(type=1, target=", .mv$state, ", Tlag=", monolixTlag(.mv$state), ", p=", monolixP(.mv$state), ")"), collapse="\n"),
         "\n\nEQUATION:\n",gsub("\n\n+", "\n",.mod),
         "\n\nOUTPUT:\n",
         "output={", paste(names(.resMod), collapse=", "), "}\n"
         )

  .lst <- monolixDataFile(.lst, uif, data, control=control)
  .definition <- .toMonolixDefinition(body(uif$saem.pars), uif$mu.ref)
  .df <- .definition[[2]]
  .dft <- as.data.frame(uif$ini)[, c("name", "est", "fix")]
  names(.dft) <- c("theta", "thetaEst", "thetaFixed")
  .df <- merge(.df, .dft, by="theta")
  .w <- which(.df$trans == "logNormal")
  if (length(.w) > 0) .df$thetaEst[.w] <- exp(.df$thetaEst[.w])
  .w <- which(.df$trans == "logitNormal")
  if (length(.w) > 0) .df$thetaEst[.w] <- sapply(.w, function(.i){RxODE::expit(.df$thetaEst[.i], .df$low[.i], .df$hi[.i])})
  .w <- which(.df$trans == "probitNormal")
  if (length(.w) > 0) .df$thetaEst[.w] <- sapply(.w, function(.i){RxODE::probitInv(.df$thetaEst[.i])})
  .dft <- as.data.frame(uif$ini)[, c("name", "est", "fix")]
  names(.dft) <- c("eta", "sdEst", "sdFixed")
  .df <- merge(.df, .dft, by="eta", all.x=TRUE)
  .df$sdEst <- sqrt(.df$sdEst)
  .lst$df <- .df
  .dfError <- do.call(rbind, .monolixGetErr)
  .dft <- as.data.frame(uif$ini)[, c("name", "est", "fix")]
  .dfError <- merge(.dfError, .dft, by="name")
  .lst$dfError <- .dfError
  .lst$parameter <- monolixModelParameter(.df, .dfError)
  .vals <- c(.df$typical, .df$sd)
  .vals <- .vals[!is.na(.vals)]
  .lst$model <- paste0("<MODEL>\n\n[INDIVIDUAL]\ninput = {", paste(.vals, collapse=", "), "}\n\n",
                       .definition[[1]], "\n\n[LONGITUDINAL]\ninput={", paste(.monolixErrs, collapse=", "), "}\n\nfile='", .lst$file,
                       ".txt'\n\nDEFINITION:\n",
                       .def)
  .lst$monolix <- paste0("<MONOLIX>\n\n[TASKS]\npopulationParameters()\nindividualParameters(method = {conditionalMean, conditionalMode})\nfim(method = Linearization)\nlogLikelihood(method = Linearization)\nplotResult(method = {outputplot, indfits, obspred, residualsscatter, residualsdistribution, parameterdistribution, covariatemodeldiagnosis, randomeffects, covariancemodeldiagnosis, saemresults })\n\n[SETTINGS]\nGLOBAL:\nexportpath = '", file.path(getwd(), .lst$file), "'\n\nPOPULATION:\nexploratoryautostop = ",
                         ifelse(control$exploratoryautostop, "yes", "no"),
                         "\nsmoothingautostop = ", ifelse(control$smoothingautostop, "yes", "no"),
                         "\nburniniterations = ", control$burniniterations,
                         "\nexploratoryiterations = ", control$exploratoryiterations,
                         "\nsimulatedannealingiterations = ", control$simulatedannealingiterations,
                         "\nsmoothingiterations = ", control$smoothingiterations,
                         "\nexploratoryalpha = ", control$exploratoryalpha,
                         "\nsimulatedannealingiterations = ", control$simulatedannealingiterations,
                         ifelse(control$variability == "none", "", paste0("\nvariability = ", control$variability)),
                         "\nexploratoryinterval = ", control$exploratoryinterval,
                         "\nomegatau = ", control$omegatau,
                         "\nerrormodeltau = ", control$errormodeltau,
                         "\n")

  .lst$fit <- monolixModelFit(uif, .lst$obs)

  .lst$mlxtran <- paste0(.lst$datafile, "\n", .lst$model, .lst$fit,
                         .lst$parameter, .lst$monolix)
  return(.lst)
}
.nlmixrMonolixLastProject <- ""

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
  if (!inherits(uif, "nlmixrUI")) {
    uif <- nlmixr::nlmixr(uif)
  }
  if(!any(tolower(names(data)) == "dv")) stop("need dv in data", call.=FALSE)
  .wid <- which(tolower(names(data)) == "id")
  names(data)[.wid] <- "id"
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
          .namem <- eval(parse(text=paste0("rxToMonolix(", .name, ")")))
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
        .ctl <- do.call(nlmixr::foceiControl, .ctl)
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
        .env$extra <- " (babelmixr)"
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
      .status <- getOption("babelmixr.monolix.status", "")
      if (.status != "") {
        system(.status)
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
