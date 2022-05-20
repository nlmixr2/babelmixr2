#' Monolix Controller for nlmixr
#'
#' @param nbSSDoses Number of steady state doses (default 7)
#' @param stiff boolean for using the stiff ODE solver
#' @param exploratoryautostop logical to turn on or off exploratory
#'   phase auto-stop of SAEM (default 250)
#' @param exploratoryiterations Number of iterations for exploratory
#'   phase (default 250)
#' @param exploratoryinterval Minimum number of interation in the
#'   exploratory phase (default 200)
#' @param simulatedannealingiterations Number of burn in iterations
#' @param runCommand is a shell command to run monolix; You can
#'   specfy the default by
#'   \code{options("babelmixr2.monolix"="runMonolix \%s")} where the \code{"\%s"}
#'   represents the monolix project file.
#' @inheritParams nlmixr2est::foceiControl
#' @return A monolix control object
#' @author Matthew Fidler
#' @export
#' @importFrom nlmixr2 nlmixr2
#' @importFrom methods is
#' @importFrom stats na.omit setNames
#' @importFrom utils assignInMyNamespace
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
                           runCommand=getOption("babelmixr2.monolix", ""),
                           adjObf=TRUE,
                           rxControl=NULL,
                           sumProd = FALSE,
                           optExpression = TRUE,
                           calcTables = TRUE,
                           compress = TRUE,
                           ci = 0.95,
                           sigdigTable=NULL, ...) {

  checkmate::assertLogical(stiff, max.len=1)
  checkmate::assertLogical(exploratoryautostop, max.len=1)
  checkmate::assertLogical(smoothingautostop, max.len=1)
  checkmate::assertLogical(simulatedannealing, max.len=1)

  checkmate::assertIntegerish(burniniterations, max.len=1, lower=1)
  checkmate::assertIntegerish(exploratoryiterations, max.len=1, lower=1)
  checkmate::assertIntegerish(simulatedannealingiterations, max.len=1, lower=1)
  checkmate::assertIntegerish(nbSSDoses, lower=7, max.len=1)
  checkmate::assertIntegerish(exploratoryiterations, max.len=1, lower=1)
  checkmate::assertIntegerish(exploratoryinterval, max.len=1, lower=1)
  checkmate::assertIntegerish(smoothingiterations, max.len=1, lower=1)
  checkmate::assertIntegerish(optimizationiterations, max.len=1, lower=1)

  checkmate::assertNumeric(exploratoryalpha, lower=0.0, upper=1.0)
  checkmate::assertNumeric(omegatau, lower=0.0, upper=1.0)
  checkmate::assertNumeric(errormodeltau, lower=0.0, upper=1.0)
  checkmate::assertNumeric(optimizationtolerance, lower=0.0)
  if (optimizationtolerance == 0) stop("'optimizationtolerance' has to be above zero")

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
      rxControl <- rxode2::rxControl()
      genRxControl <- TRUE
    } else if (is.list(rxControl)) {
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
               runCommand=runCommand,
               adjObf=adjObf,
               rxControl=rxControl,
               sumProd = sumProd,
               optExpression=optExpression,
               calcTables = calcTables,
               compress = compress,
               ci = ci,
               sigdigTable=sigdigTable,
               genRxControl=genRxControl)
  class(.ret) <- "monolixControl"
  .ret
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
                                            sigdigTable = .monolixControl$sigdigTable,
                                            maxSS=.monolixControl$nbSSDoses + 1,
                                            minSS=.monolixControl$nbSSDoses,
                                            ssAtol=100,
                                            ssRtol=100,
                                            atol=ifelse(.monolixControl$stiff, 1e-9, 1e-6),
                                            rtol=ifelse(.monolixControl$stiff, 1e-6, 1e-3),
                                            method=ifelse(.monolixControl$stiff, "liblsoda",
                                                          "dop853"))
  if (assign)
    env$control <- .foceiControl
  .foceiControl
}

#' @importFrom nlmixr2est nmObjGetFoceiControl
#' @export
nlmixr2est::nmObjGetFoceiControl

#' @export
#' @rdname nmObjGetFoceiControl
nmObjGetFoceiControl.monolix <- function(x, ...) {
  .monolixControlToFoceiControl(x[[1]])
}


#' @importFrom nlmixr2est getValidNlmixrCtl
#' @export
nlmixr2est::getValidNlmixrCtl

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

#' @importFrom nlmixr2est nmObjHandleControlObject
#' @export
nlmixr2est::nmObjHandleControlObject

#' @export
nmObjHandleControlObject.monolixControl <- function(control, env) {
  assign("monolixControl", control, envir=env)
}

#' @importFrom nlmixr2est nmObjGetControl
#' @export
nlmixr2est::nmObjGetControl

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
