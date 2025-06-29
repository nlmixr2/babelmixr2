#' Monolix Controller for nlmixr2
#'
#' @param nbSSDoses Number of steady state doses (default 7)
#' @param stiff boolean for using the stiff ODE solver
#' @param exploratoryAutoStop logical to turn on or off exploratory
#'   phase auto-stop of SAEM (default 250)
#' @param exploratoryIterations Number of iterations for exploratory
#'   phase (default 250)
#' @param exploratoryInterval Minimum number of iterations in the
#'   exploratory phase (default 200)
#' @param exploratoryAlpha Convergence memory in the exploratory phase
#'   (only used when `exploratoryAutoStop` is `TRUE`)
#' @param simulatedAnnealingIterations Number of simulating annealing
#'   iterations
#' @param burnInIterations Number of burn in iterations
#' @param smoothingIterations Number of smoothing iterations
#' @param smoothingAutoStop Boolean indicating if the smoothing should
#'   automatically stop (default `FALSE`)
#' @param useLinearization Use linearization for log likelihood and
#'   fim.
#' @inheritParams nlmixr2est::foceiControl
#' @param omegaTau Proportional rate on variance for simulated
#'   annealing
#' @param errorModelTau Proportional rate on error model for simulated
#'   annealing
#' @param variability This describes the methodology for parameters
#'   without variability.  It could be: - Fixed throughout (none) -
#'   Variability in the first stage (firstStage) - Decreasing until it
#'   reaches the fixed value (decreasing)
#' @param runCommand is a shell command or function to run monolix; You can specify
#'   the default by
#'   \code{options("babelmixr2.monolix"="runMonolix")}. If it is empty
#'   and 'lixoftConnectors' is available, use lixoftConnectors to run
#'   monolix. See details for function usage.
#' @param absolutePath Boolean indicating if the absolute path should
#'   be used for the monolix runs
#'
#' @param run Should monolix be run and the results be imported to nlmixr2?  (Default is TRUE)
#' @inheritParams nonmemControl
#' @inheritParams nlmixr2est::saemControl
#' @return A monolix control object
#' @author Matthew Fidler
#' @details
#'
#' If \code{runCommand} is given as a string, it will be called with the
#' \code{system()} command like:
#'
#' \code{runCommand mlxtran}.
#'
#' For example, if \code{runCommand="'/path/to/monolix/mlxbsub2021' -p "} then the command line
#' used would look like the following:
#'
#' \code{'/path/to/monolix/mlxbsub2021' monolix.mlxtran}
#'
#' If \code{runCommand} is given as a function, it will be called as
#' \code{FUN(mlxtran, directory, ui)} to run Monolix.  This allows you to run Monolix
#' in any way that you may need, as long as you can write it in R.  babelmixr2
#' will wait for the function to return before proceeding.
#'
#' If \code{runCommand} is \code{NA}, \code{nlmixr()} will stop after writing
#' the model files and without starting Monolix.
#'
#' Note that you can get the translated monolix components from a
#' parsed/compiled rxode2 ui object with `ui$monolixModel` and `ui$mlxtran`
#'
#' @export
#' @importFrom nlmixr2est nlmixr2
#' @importFrom methods is
#' @importFrom stats na.omit setNames
#' @importFrom utils assignInMyNamespace read.csv write.csv
#' @importFrom rxode2 `model<-`
monolixControl <- function(nbSSDoses=7,
                           useLinearization=FALSE,
                           stiff=FALSE,
                           addProp = c("combined2", "combined1"),
                           exploratoryAutoStop=FALSE,
                           smoothingAutoStop=FALSE,
                           burnInIterations=5,
                           smoothingIterations=200,
                           exploratoryIterations=250,
                           simulatedAnnealingIterations=250,
                           exploratoryInterval=200,
                           exploratoryAlpha=0.0,
                           omegaTau=0.95,
                           errorModelTau=0.95,
                           variability=c("none", "firstStage", "decreasing"),
                           runCommand=getOption("babelmixr2.monolix", ""),
                           rxControl=NULL,
                           sumProd = FALSE,
                           optExpression = TRUE,
                           calcTables = TRUE,
                           compress = TRUE,
                           ci = 0.95,
                           sigdigTable=NULL,
                           absolutePath=FALSE,
                           modelName=NULL,
                           muRefCovAlg=TRUE,
                           run=TRUE,
                           ...) {
  checkmate::assertLogical(stiff, max.len=1, any.missing=FALSE)
  checkmate::assertLogical(exploratoryAutoStop, max.len=1, any.missing=FALSE)
  checkmate::assertLogical(smoothingAutoStop, max.len=1, any.missing=FALSE)
  checkmate::assertLogical(absolutePath, max.len=1, any.missing=FALSE)
  checkmate::assertLogical(useLinearization, max.len=1, any.missing=FALSE)

  checkmate::assertIntegerish(burnInIterations, max.len=1, lower=1)
  checkmate::assertIntegerish(exploratoryIterations, max.len=1, lower=1)
  checkmate::assertIntegerish(simulatedAnnealingIterations, max.len=1, lower=1)
  checkmate::assertIntegerish(nbSSDoses, lower=7, max.len=1)
  checkmate::assertIntegerish(exploratoryInterval, max.len=1, lower=1)
  checkmate::assertIntegerish(smoothingIterations, max.len=1, lower=1)

  checkmate::assertNumeric(exploratoryAlpha, lower=0.0, upper=1.0)
  checkmate::assertNumeric(omegaTau, lower=0.0, upper=1.0)
  checkmate::assertNumeric(errorModelTau, lower=0.0, upper=1.0)
  checkmate::assertLogical(run, len=1, any.missing=FALSE)

  if (!is.null(modelName)) {
    checkmate::assertCharacter(modelName, len=1, any.missing=FALSE)
  }

  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% c("genRxControl"))]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste
    (paste0("'", .bad, "'", sep=""), collapse=", "),
    call.=FALSE)
  }

  if (checkmate::testIntegerish(addProp, lower=1, upper=1, len=1)) {
    addProp <- c("combined1", "combined2")[addProp]
  } else {
    addProp <- match.arg(addProp)
  }
  checkmate::assertLogical(compress, any.missing=FALSE, len=1)

  if (!is.null(.xtra$genRxControl)) {
    genRxControl <- .xtra$genRxControl
  } else {
    genRxControl <- FALSE
    if (is.null(rxControl)) {
      rxControl <- rxode2::rxControl(
        maxSS=nbSSDoses + 1,
        minSS=nbSSDoses,
        ssAtol=100,
        ssRtol=100,
        atol=ifelse(stiff, 1e-9, 1e-6),
        rtol=ifelse(stiff, 1e-6, 1e-3),
        method=ifelse(stiff, "liblsoda", "dop853")
      )
      genRxControl <- TRUE
    } else if (is.list(rxControl)) {
        rxControl$maxSS <- nbSSDoses + 1
        rxControl$minSS <- nbSSDoses
        rxControl$ssAtol <- 100
        rxControl$ssRtol <- 100
        rxControl$atol <- ifelse(stiff, 1e-9, 1e-6)
        rxControl$rtol <- ifelse(stiff, 1e-6, 1e-3)
        rxControl$method <- ifelse(stiff, "liblsoda", "dop853")
      rxControl <- do.call(rxode2::rxControl, rxControl)
    }
    if (!inherits(rxControl, "rxControl")) {
      stop("rxControl needs to be ode solving options from rxode2::rxControl()",
           call.=FALSE)
    }
  }

  checkmate::assertLogical(sumProd, any.missing=FALSE, len=1)
  checkmate::assertLogical(optExpression, any.missing=FALSE, len=1)
  checkmate::assertNumeric(ci, any.missing=FALSE, len=1, lower=0, upper=1)
  checkmate::assertLogical(calcTables, len=1, any.missing=FALSE)

  if (!identical(runCommand, "")) {
    if (!(checkmate::testCharacter(runCommand, len=1) ||
            checkmate::testFunction(runCommand, args=c("ctl", "directory", "ui")))) {
      stop("runCommand must be a character string or a function with arguments 'ctl', 'directory', and 'ui'")
    }
  }

  .ret <- list(nbSSDoses=as.integer(nbSSDoses), stiff=stiff,
               exploratoryAutoStop=exploratoryAutoStop,
               smoothingAutoStop=smoothingAutoStop,
               addProp=addProp,
               burnInIterations=burnInIterations,
               simulatedAnnealingIterations=simulatedAnnealingIterations,
               exploratoryInterval=exploratoryInterval,
               smoothingIterations=smoothingIterations,
               exploratoryAlpha=exploratoryAlpha,
               omegaTau=omegaTau,
               errorModelTau=errorModelTau,
               exploratoryIterations=exploratoryIterations,
               variability=match.arg(variability),
               runCommand=runCommand,
               rxControl=rxControl,
               sumProd = sumProd,
               optExpression=optExpression,
               calcTables = calcTables,
               compress = compress,
               ci = ci,
               sigdigTable=sigdigTable,
               genRxControl=genRxControl,
               useLinearization=useLinearization,
               modelName=modelName,
               muRefCovAlg=muRefCovAlg,
               run=run)
  class(.ret) <- "monolixControl"
  .ret
}

rxUiDeparse.monolixControl <- function(object, var) {
  .default <- monolixControl()
  .w <- nlmixr2est::.deparseDifferent(.default, object, "genRxControl")
  nlmixr2est::.deparseFinal(.default, object, .w, var)
}

.monolixControlToFoceiControl <- function(env, assign = TRUE) {
  .monolixControl <- env$monolixControl
  .ui <- env$ui
  .foceiControl <- nlmixr2est::foceiControl(rxControl = env$monolixControl$rxControl,
                                            maxOuterIterations = 0L, maxInnerIterations = 0L, covMethod = 0L,
                                            etaMat = env$etaMat, sumProd = .monolixControl$sumProd,
                                            optExpression = .monolixControl$optExpression, scaleTo = 0,
                                            calcTables = .monolixControl$calcTables,
                                            addProp = .monolixControl$addProp,
                                            skipCov = .ui$foceiSkipCov, interaction = 1L,
                                            compress = .monolixControl$compress,
                                            ci = .monolixControl$ci,
                                            sigdigTable = .monolixControl$sigdigTable)
  if (assign)
    env$control <- .foceiControl
  .foceiControl
}

#' @export
nmObjGetFoceiControl.monolix <- function(x, ...) {
  .monolixControlToFoceiControl(x[[1]])
}

#' @export
getValidNlmixrCtl.monolix <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- monolixControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("monolixControl", .ctl)
  if (!inherits(.ctl, "monolixControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- monolixControl()
  } else {
    .ctl <- do.call(monolixControl, .ctl)
  }
  .ctl
}

#' @export
nmObjHandleControlObject.monolixControl <- function(control, env) {
  assign("monolixControl", control, envir=env)
}

#' @export
nmObjGetControl.monolix <- function(x, ...) {
  .env <- x[[1]]
  if (exists("monolixControl", .env)) {
    .control <- get("monolixControl", .env)
    if (inherits(.control, "monolixControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "monolixControl")) return(.control)
  }
  stop("cannot find monolix related control object", call.=FALSE)
}
