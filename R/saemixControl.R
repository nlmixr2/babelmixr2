#' Control for saemix estimation method in nlmixr2
#'
#' @inheritParams nlmixr2est::iterPrintParams
#' @inheritParams nlmixr2est::foceiControl
#' @inheritParams nlmixr2est::saemControl
#'
#' @param map logical, whether to compute MAP estimates (default `TRUE`)
#' @param fim logical, whether to compute FIM (default `TRUE`)
#' @param ll.is logical, whether to compute log-likelihood by Importance Sampling (default `TRUE`)
#' @param ll.gq logical, whether to compute log-likelihood by Gaussian Quadrature (default `FALSE`)
#' @param nbiter.saemix integer vector of length 2, number of iterations for phase 1 (exploratory) and phase 2 (smoothing) of SAEM (default `c(300, 100)`)
#' @param nbiter.sa number of iterations for simulated annealing (default `NA`)
#' @param nbiter.burn number of iterations for burn-in (default `5`)
#' @param nbiter.map number of iterations for MAP estimation (default `5`)
#' @param nb.chains number of chains (default `1`)
#' @param fix.seed logical, whether to fix the random seed (default `TRUE`)
#' @param seed random seed (default `23456`)
#' @param nmc.is number of Monte Carlo samples for Importance Sampling (default `5000`)
#' @param nu.is number of degrees of freedom for the student distribution in Importance Sampling (default `4`)
#' @param print.is logical, whether to print progress during Importance Sampling (default `FALSE`)
#' @param nbdisplay number of iterations between progress displays (default `100`)
#' @param displayProgress logical, whether to display graphical progress (default `FALSE`)
#' @param print whether `saemix` prints algorithm progress (default
#'   `FALSE`).  Accepts the legacy logical/integer scalar or a
#'   pre-built [nlmixr2est::iterPrintControl()] object; internally the
#'   control stores a single `iterPrintControl` sub-list (matching the
#'   `nlmixr2est` iteration-print unification), and any nonzero
#'   `every` enables the `saemix` progress output
#' @param save logical, whether to save results to files (default `TRUE`)
#' @param save.graphs logical, whether to save graphs (default `TRUE`)
#' @param directory directory where results and graphs are saved (default `"newdir"`)
#' @param warnings logical, whether to show warnings (default `FALSE`)
#' @param nbiter.mcmc integer vector of length 4, number of iterations for MCMC kernel updates (default `c(2, 2, 2, 0)`)
#' @param proba.mcmc probability for MCMC kernel selection (default `0.4`)
#' @param stepsize.rw stepsize for random walk kernel (default `0.4`)
#' @param rw.init initial standard deviation for random walk kernel (default `0.5`)
#' @param alpha.sa parameter for simulated annealing (default `0.97`)
#' @param nnodes.gq number of nodes for Gaussian Quadrature (default `12`)
#' @param nsd.gq number of standard deviations for Gaussian Quadrature (default `4`)
#' @param maxim.maxiter maximum number of iterations for maximization step (default `100`)
#' @param nb.sim number of simulations for visual predictive check (default `1000`)
#' @param nb.simpred number of simulations for predictions (default `100`)
#' @param ipar.lmcmc parameter for L-MCMC (default `50`)
#' @param ipar.rmcmc parameter for R-MCMC (default `0.05`)
#'
#' @return saemix control structure
#' @export
#' @author Matthew L. Fidler & Antigravity
#' @examples
#' \donttest{
#' # Example showing control creation:
#' ctrl <- saemixControl(seed = 123456)
#' }
saemixControl <- function(map = TRUE,
                          fim = TRUE,
                          ll.is = TRUE,
                          ll.gq = FALSE,
                          nbiter.saemix = c(300, 100),
                          nbiter.sa = NA,
                          nbiter.burn = 5,
                          nbiter.map = 5,
                          nb.chains = 1,
                          fix.seed = TRUE,
                          seed = 23456,
                          nmc.is = 5000,
                          nu.is = 4,
                          print.is = FALSE,
                          nbdisplay = 100,
                          displayProgress = FALSE,
                          print = FALSE,
                          save = FALSE,
                          save.graphs = TRUE,
                          directory = "newdir",
                          warnings = FALSE,
                          nbiter.mcmc = c(2, 2, 2, 0),
                          proba.mcmc = 0.4,
                          stepsize.rw = 0.4,
                          rw.init = 0.5,
                          alpha.sa = 0.97,
                          nnodes.gq = 12,
                          nsd.gq = 4,
                          maxim.maxiter = 100,
                          nb.sim = 1000,
                          nb.simpred = 100,
                          ipar.lmcmc = 50,
                          ipar.rmcmc = 0.05,

                          # rxode2/nlmixr2est options
                          rxControl = NULL,
                          stickyRecalcN = 4,
                          maxOdeRecalc = 5,
                          odeRecalcFactor = 10^(0.5),

                          useColor = NULL,
                          printNcol = NULL,

                          normType = c("rescale2", "mean", "rescale", "std", "len", "constant"),
                          scaleType = c("none", "nlmixr2", "norm", "mult", "multAdd"),
                          scaleCmax = 1e5,
                          scaleCmin = 1e-5,
                          scaleC = NULL,
                          scaleTo = 1.0,

                          addProp = c("combined2", "combined1"),
                          calcTables = TRUE,
                          compress = TRUE,
                          ci = 0.95,
                          sigdigTable = NULL,
                          sigdig = 4,
                          ...) {

  checkmate::assertLogical(map, len = 1, any.missing = FALSE)
  checkmate::assertLogical(fim, len = 1, any.missing = FALSE)
  checkmate::assertLogical(ll.is, len = 1, any.missing = FALSE)
  checkmate::assertLogical(ll.gq, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(nbiter.saemix, len = 2, any.missing = FALSE)
  if (!is.na(nbiter.sa)) {
    checkmate::assertIntegerish(nbiter.sa, len = 1)
  }
  checkmate::assertIntegerish(nbiter.burn, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(nbiter.map, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(nb.chains, len = 1, any.missing = FALSE)
  checkmate::assertLogical(fix.seed, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(seed, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(nmc.is, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(nu.is, len = 1, any.missing = FALSE)
  checkmate::assertLogical(print.is, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(nbdisplay, len = 1, any.missing = FALSE)
  checkmate::assertLogical(displayProgress, len = 1, any.missing = FALSE)
  checkmate::assertLogical(save, len = 1, any.missing = FALSE)
  checkmate::assertLogical(save.graphs, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(directory, len = 1, any.missing = FALSE)
  checkmate::assertLogical(warnings, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(nbiter.mcmc, len = 4, any.missing = FALSE)
  checkmate::assertNumeric(proba.mcmc, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(stepsize.rw, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(rw.init, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(alpha.sa, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(nnodes.gq, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(nsd.gq, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(maxim.maxiter, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(nb.sim, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(nb.simpred, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(ipar.lmcmc, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(ipar.rmcmc, len = 1, any.missing = FALSE)

  checkmate::assertLogical(compress, any.missing = FALSE, len = 1)
  checkmate::assertLogical(calcTables, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(ci, any.missing = FALSE, len = 1, lower = 0, upper = 1)

  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% c("genRxControl", "iterPrintControl"))]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste(paste0("'", .bad, "'"), collapse = ", "), call. = FALSE)
  }

  # `print` was historically a logical for this method; coerce so the
  # shared absorb helper (which requires an integerish `every`) accepts it
  if (checkmate::testLogical(print, len = 1, any.missing = FALSE)) {
    print <- as.integer(print)
  }
  .iterPrintControl <- nlmixr2est::.absorbIterPrintControl(print = print,
                                                           printNcol = printNcol,
                                                           useColor = useColor,
                                                           iterPrintControl = .xtra$iterPrintControl)

  .genRxControl <- FALSE
  if (!is.null(.xtra$genRxControl)) {
    .genRxControl <- .xtra$genRxControl
  }
  if (is.null(rxControl)) {
    if (!is.null(sigdig)) {
      rxControl <- rxode2::rxControl(sigdig = sigdig)
    } else {
      rxControl <- rxode2::rxControl(atol = 1e-4, rtol = 1e-4)
    }
    .genRxControl <- TRUE
  } else if (inherits(rxControl, "rxControl")) {
  } else if (is.list(rxControl)) {
    rxControl <- do.call(rxode2::rxControl, rxControl)
  } else {
    stop("solving options 'rxControl' needs to be generated from 'rxode2::rxControl'", call. = FALSE)
  }

  if (checkmate::testIntegerish(addProp, lower = 1, upper = 2, len = 1)) {
    addProp <- c("combined2", "combined1")[addProp]
  } else {
    addProp <- match.arg(addProp)
  }

  if (is.null(sigdigTable)) {
    sigdigTable <- 3
  }
  checkmate::assertIntegerish(sigdigTable, lower = 1, len = 1, any.missing = FALSE)

  if (checkmate::testIntegerish(scaleType, len = 1, lower = 1, upper = 5, any.missing = FALSE)) {
    scaleType <- as.integer(scaleType)
  } else {
    .scaleTypeIdx <- c("norm" = 1L, "nlmixr2" = 2L, "mult" = 3L, "multAdd" = 4L, "none" = 5L)
    scaleType <- setNames(.scaleTypeIdx[match.arg(scaleType)], NULL)
  }

  .normTypeIdx <- c("rescale2" = 1L, "rescale" = 2L, "mean" = 3L, "std" = 4L, "len" = 5L, "constant" = 6L)
  if (checkmate::testIntegerish(normType, len = 1, lower = 1, upper = 6, any.missing = FALSE)) {
    normType <- as.integer(normType)
  } else {
    normType <- setNames(.normTypeIdx[match.arg(normType)], NULL)
  }
  checkmate::assertNumeric(scaleCmax, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(scaleCmin, lower = 0, any.missing = FALSE, len = 1)
  if (!is.null(scaleC)) {
    checkmate::assertNumeric(scaleC, lower = 0, any.missing = FALSE)
  }
  checkmate::assertNumeric(scaleTo, len = 1, lower = 0, any.missing = FALSE)

  .ret <- list(map = map,
               fim = fim,
               ll.is = ll.is,
               ll.gq = ll.gq,
               nbiter.saemix = nbiter.saemix,
               nbiter.sa = nbiter.sa,
               nbiter.burn = nbiter.burn,
               nbiter.map = nbiter.map,
               nb.chains = nb.chains,
               fix.seed = fix.seed,
               seed = seed,
               nmc.is = nmc.is,
               nu.is = nu.is,
               print.is = print.is,
               nbdisplay = nbdisplay,
               displayProgress = displayProgress,
               save = save,
               save.graphs = save.graphs,
               directory = directory,
               warnings = warnings,
               nbiter.mcmc = nbiter.mcmc,
               proba.mcmc = proba.mcmc,
               stepsize.rw = stepsize.rw,
               rw.init = rw.init,
               alpha.sa = alpha.sa,
               nnodes.gq = nnodes.gq,
               nsd.gq = nsd.gq,
               maxim.maxiter = maxim.maxiter,
               nb.sim = nb.sim,
               nb.simpred = nb.simpred,
               ipar.lmcmc = ipar.lmcmc,
               ipar.rmcmc = ipar.rmcmc,

               # rxode2/nlmixr2est options
               rxControl = rxControl,
               stickyRecalcN = as.integer(stickyRecalcN),
               maxOdeRecalc = as.integer(maxOdeRecalc),
               odeRecalcFactor = odeRecalcFactor,
               iterPrintControl = .iterPrintControl,
               normType = normType,
               scaleType = scaleType,
               scaleCmax = scaleCmax,
               scaleCmin = scaleCmin,
               scaleC = scaleC,
               scaleTo = scaleTo,
               addProp = addProp,
               calcTables = calcTables,
               compress = compress,
               ci = ci,
               sigdigTable = sigdigTable,
               sigdig = sigdig,
               genRxControl = .genRxControl)

  class(.ret) <- "saemixControl"
  .ret
}

rxUiDeparse.saemixControl <- function(object, var) {
  .default <- saemixControl()
  .w <- nlmixr2est::.deparseDifferent(.default, object, "genRxControl")
  nlmixr2est::.deparseFinal(.default, object, .w, var)
}

#' @export
nmObjHandleControlObject.saemixControl <- function(control, env) {
  assign("saemixControl", control, envir = env)
}

#' @export
nmObjGetControl.saemix <- function(x, ...) {
  .env <- x[[1]]
  if (exists("saemixControl", .env)) {
    .control <- get("saemixControl", .env)
    if (inherits(.control, "saemixControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "saemixControl")) return(.control)
  }
  stop("cannot find saemix related control object", call. = FALSE)
}

#' @export
getValidNlmixrCtl.saemix <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- saemixControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("saemixControl", .ctl)
  if (!inherits(.ctl, "saemixControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- saemixControl()
  } else {
    .ctl <- do.call(saemixControl, .ctl)
  }
  .ctl
}
