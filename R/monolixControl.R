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
#' @inheritParams nlmixr2::foceiControl
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
                           runCommand=getOption("babelmixr2.monolix", "")) {

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
