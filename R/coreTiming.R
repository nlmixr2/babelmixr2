#' Warn that coreTiming() cannot set an external engine's threads
#'
#' NONMEM and Monolix take their parallel setup from `runCommand`, not from
#' rxode2, so only the R side of a timed run changes with `cores`.
#'
#' @param engine engine name for the message
#' @return a `function(control, cores)` for `coreTimingRegister(cores=)`
#' @noRd
#' @author Matthew L. Fidler
.bblCoreTimingExternalCores <- function(engine) {
  force(engine)
  function(control, cores) {
    if (is.list(control$rxControl)) {
      control$rxControl$cores <- as.integer(cores)
    }
    warning(engine, " sets its own threads through 'runCommand'; ",
            "coreTiming() varies only the R-side threads",
            call. = FALSE)
    control
  }
}

#' Register babelmixr2's estimation methods with nlmixr2extra::coreTiming()
#'
#' Each registration names the control options that cap the method's
#' iterations, so `coreTiming()` can time short runs.  Skipped when the
#' installed nlmixr2extra predates `coreTimingRegister()`.
#'
#' @return `TRUE` when registered, `FALSE` otherwise, invisibly
#' @noRd
#' @author Matthew L. Fidler
.bblCoreTimingRegister <- function() {
  .ns <- asNamespace("nlmixr2extra")
  if (!exists("coreTimingRegister", envir = .ns, mode = "function", inherits = FALSE)) {
    return(invisible(FALSE))
  }
  .reg <- get("coreTimingRegister", envir = .ns, mode = "function")
  # NONMEM: evaluations of the classical methods; iterations of imp/its
  .reg("nonmemControl", args = c("maxeval", "niter"),
       cores = .bblCoreTimingExternalCores("NONMEM"))
  # Monolix: every phase of SAEM, with the auto-stops off so the cap holds
  .reg("monolixControl",
       fun = function(control, nIter) {
         .n <- as.integer(nIter)
         control$burnInIterations <- min(control$burnInIterations, .n)
         control$simulatedAnnealingIterations <- min(control$simulatedAnnealingIterations, .n)
         control$exploratoryIterations <- .n
         control$smoothingIterations <- .n
         control$exploratoryInterval <- min(control$exploratoryInterval, .n)
         control$exploratoryAutoStop <- FALSE
         control$smoothingAutoStop <- FALSE
         control
       },
       cores = .bblCoreTimingExternalCores("Monolix"))
  # saemix: c(K1, K2) SAEM iterations
  .reg("saemixControl", fun = function(control, nIter) {
    .n <- as.integer(nIter)
    control$nbiter.saemix <- c(.n, .n)
    control$nbiter.burn <- min(control$nbiter.burn, .n)
    control
  })
  # nlmer: optimizer evaluations; bobyqa needs more than 2 * npar + 1
  .reg("nlmerControl", fun = function(control, nIter) {
    .n <- max(as.integer(nIter), 50L)
    if (identical(control$optimizer, "nloptwrap")) {
      control$optCtrl$maxeval <- .n
    } else {
      control$optCtrl$maxfun <- .n
    }
    control
  })
  # FME MCMC: chain length, with burn-in and stored output no longer than it
  .reg("fmeMcmcControl", fun = function(control, nIter) {
    .n <- as.integer(nIter)
    control$niter <- .n
    control$burninlength <- min(control$burninlength, .n - 1L)
    control$outputlength <- min(control$outputlength, .n)
    control$updatecov <- min(control$updatecov, .n)
    control
  })
  .reg("pseudoOptimControl", args = "numiter")
  invisible(TRUE)
}
