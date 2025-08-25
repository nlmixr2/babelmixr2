#' Control for fmeMcmc estimation method in nlmixr2
#'
#' @inheritParams nlmixr2est::foceiControl
#' @inheritParams nlmixr2est::saemControl
#'
#' @param npop Number of elements in the population. Defaults to
#'   max(5*length(p),50) which is calculated from the number of parameters in the model
#'
#' @param numiter Number of iterations to run the optimization.
#'   Defaults to 10000. The algorithm either stops when `numiter`
#'   iterations has been performed or when the remaining variation is
#'   less than `varleft`.
#'
#' @param centroid Number of elements from which to estimate a new
#'   parameter vector.  The default is 3.
#'
#' @param varleft relative variation remaining; if below this value,
#'   the algorithm stops.  Defaults to 1e-8.
#'
#' @param verbose If TRUE, print information about the optimization
#'   from `FME::pseudoOptim`.  Default is FALSE.
#'
#' @param returnPseudoOptim return the pseudoOptim output instead of
#'   the nlmixr2 fit
#'
#' @return pseudoOptim control structure
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' \donttest{
#' # A logit regression example with emax model
#'
#' dsn <- data.frame(i=1:1000)
#' dsn$time <- exp(rnorm(1000))
#' dsn$DV=rbinom(1000,1,exp(-1+dsn$time)/(1+exp(-1+dsn$time)))
#'
#' mod <- function() {
#'  ini({
#'    # This estimation method requires all parameters
#'    # to be bounded:
#'    E0 <- c(-100, 0.5, 100)
#'    Em <- c(0, 0.5, 10)
#'    E50 <- c(0, 2, 20)
#'    g <- fix(c(0.1, 2, 10))
#'  })
#'  model({
#'    v <- E0+Em*time^g/(E50^g+time^g)
#'    ll(bin) ~ DV * v - log(1 + exp(v))
#'  })
#' }
#'
#' fit2 <- nlmixr(mod, dsn, est="pseudoOptim")
#'
#' print(fit2)
#'
#'
#' }
pseudoOptimControl <- function(npop=NULL, # Number of elements in the pouplation
                               numiter=10000, # Number of iterations
                               centroid=3, # Centroid for the population
                               varleft=1e-8,
                               verbose=FALSE,
                           returnPseudoOptim=FALSE,

                           # rxode2/nlmixr2est options

                           stickyRecalcN=4,
                           maxOdeRecalc=5,
                           odeRecalcFactor=10^(0.5),

                           useColor = crayon::has_color(),
                           printNcol = floor((getOption("width") - 23) / 12), #
                           print = 1L, #

                           normType = c("rescale2", "mean", "rescale", "std", "len", "constant"), #
                           scaleType = c("none", "nlmixr2", "norm", "mult", "multAdd"), #
                           scaleCmax = 1e5, #
                           scaleCmin = 1e-5, #
                           scaleC=NULL,
                           scaleTo=1.0,

                           rxControl=NULL,
                           optExpression=TRUE, sumProd=FALSE,
                           literalFix=TRUE,
                           literalFixRes=TRUE,
                           addProp = c("combined2", "combined1"),
                           calcTables=TRUE, compress=TRUE,
                           covMethod=c("r", ""),
                           adjObf=TRUE, ci=0.95, sigdig=4, sigdigTable=NULL, ...) {


  checkmate::assertLogical(optExpression, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFix, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFixRes, len=1, any.missing=FALSE)
  checkmate::assertLogical(sumProd, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnPseudoOptim, len=1, any.missing=FALSE)
  checkmate::assertLogical(calcTables, len=1, any.missing=FALSE)
  checkmate::assertLogical(compress, len=1, any.missing=TRUE)
  checkmate::assertLogical(adjObf, len=1, any.missing=TRUE)

  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad == "genRxControl")]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste
    (paste0("'", .bad, "'", sep=""), collapse=", "),
    call.=FALSE)
  }

  checkmate::assertIntegerish(stickyRecalcN, any.missing=FALSE, lower=0, len=1)
  checkmate::assertIntegerish(maxOdeRecalc, any.missing=FALSE, len=1)
  checkmate::assertNumeric(odeRecalcFactor, len=1, lower=1, any.missing=FALSE)

  .genRxControl <- FALSE
  if (!is.null(.xtra$genRxControl)) {
    .genRxControl <- .xtra$genRxControl
  }
  if (is.null(rxControl)) {
    if (!is.null(sigdig)) {
      rxControl <- rxode2::rxControl(sigdig=sigdig)
    } else {
      rxControl <- rxode2::rxControl(atol=1e-4, rtol=1e-4)
    }
    .genRxControl <- TRUE
  } else if (inherits(rxControl, "rxControl")) {
  } else if (is.list(rxControl)) {
    rxControl <- do.call(rxode2::rxControl, rxControl)
  } else {
    stop("solving options 'rxControl' needs to be generated from 'rxode2::rxControl'", call=FALSE)
  }
  if (!is.null(sigdig)) {
    checkmate::assertNumeric(sigdig, lower=1, finite=TRUE, any.missing=TRUE, len=1)
    if (is.null(sigdigTable)) {
      sigdigTable <- round(sigdig)
    }
  }
  if (is.null(sigdigTable)) {
    sigdigTable <- 3
  }
  checkmate::assertIntegerish(sigdigTable, lower=1, len=1, any.missing=FALSE)

  checkmate::assertLogical(useColor, any.missing=FALSE, len=1)
  checkmate::assertIntegerish(print, len=1, lower=0, any.missing=FALSE)
  checkmate::assertIntegerish(printNcol, len=1, lower=0, any.missing=FALSE)
  if (checkmate::testIntegerish(scaleType, len=1, lower=1, upper=5, any.missing=FALSE)) {
    scaleType <- as.integer(scaleType)
  } else {
    .scaleTypeIdx <- c("norm" = 1L, "nlmixr2" = 2L, "mult" = 3L, "multAdd" = 4L,
                       "none"= 5L)
    scaleType <- setNames(.scaleTypeIdx[match.arg(scaleType)], NULL)
  }

  .normTypeIdx <- c("rescale2" = 1L, "rescale" = 2L, "mean" = 3L, "std" = 4L, "len" = 5L, "constant" = 6L)
  if (checkmate::testIntegerish(normType, len=1, lower=1, upper=6, any.missing=FALSE)) {
    normType <- as.integer(normType)
  } else {
    normType <- setNames(.normTypeIdx[match.arg(normType)], NULL)
  }
  checkmate::assertNumeric(scaleCmax, lower=0, any.missing=FALSE, len=1)
  checkmate::assertNumeric(scaleCmin, lower=0, any.missing=FALSE, len=1)
  if (!is.null(scaleC)) {
    checkmate::assertNumeric(scaleC, lower=0, any.missing=FALSE)
  }
  checkmate::assertNumeric(scaleTo, len=1, lower=0, any.missing=FALSE)

  .ret <- list(
    npop = npop,
    numiter= as.integer(numiter),
    centroid= as.integer(centroid),
    varleft = as.double(varleft),
    verbose = as.logical(verbose),
    covMethod=match.arg(covMethod),
    optExpression=optExpression,
    literalFix=literalFix,
    literalFixRes=literalFixRes,
    sumProd=sumProd,
    rxControl=rxControl,
    returnPseudoOptim=returnPseudoOptim,

    stickyRecalcN=as.integer(stickyRecalcN),
    maxOdeRecalc=as.integer(maxOdeRecalc),
    odeRecalcFactor=odeRecalcFactor,

    useColor=useColor,
    print=print,
    printNcol=printNcol,
    scaleType=scaleType,
    normType=normType,

    scaleCmax=scaleCmax,
    scaleCmin=scaleCmin,
    scaleC=scaleC,
    scaleTo=scaleTo,

    addProp=match.arg(addProp),
    calcTables=calcTables,
    compress=compress,
    ci=ci, sigdig=sigdig, sigdigTable=sigdigTable,
    genRxControl=.genRxControl)
  class(.ret) <- "pseudoOptimControl"
  .ret
}

#' Get the pseudoOptim family control
#'
#' @param env pseudoOptim optimization environment
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.pseudoOptimFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- babelmixr2::pseudoOptimControl()
  }
  if (!inherits(.control, "pseudoOptimControl")){
    .control <- do.call(babelmixr2::pseudoOptimControl, .control)
  }
  assign("control", .control, envir=.ui)
}

nmObjHandleControlObject.pseudoOptimControl <- function(control, env) {
  assign("pseudoOptimControl", control, envir=env)
}

nmObjGetControl.pseudoOptim <- function(x, ...) {
  .env <- x[[1]]
  if (exists("pseudoOptimControl", .env)) {
    .control <- get("pseudoOptimControl", .env)
    if (inherits(.control, "pseudoOptimControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "pseudoOptimControl")) return(.control)
  }
  stop("cannot find pseudoOptim related control object", call.=FALSE)
}

getValidNlmixrCtl.pseudoOptim <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- pseudoOptimControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("pseudoOptimControl", .ctl)
  if (!inherits(.ctl, "pseudoOptimControl")) {
    .minfo("invalid control for `est=\"pseudoOptim\"`, using default")
    .ctl <- pseudoOptimControl()
  } else {
    .ctl <- do.call(pseudoOptimControl, .ctl)
  }
  .ctl
}

.pseudoOptimControlToFoceiControl <- function(env, assign=TRUE) {
  .pseudoOptimControl <- env$pseudoOptimControl
  .ui <- env$ui
  .foceiControl <- foceiControl(rxControl=env$pseudoOptimControl$rxControl,
                                maxOuterIterations=0L,
                                maxInnerIterations=0L,
                                covMethod=0L,
                                sumProd=.pseudoOptimControl$sumProd,
                                optExpression=.pseudoOptimControl$optExpression,
                                literalFix=.pseudoOptimControl$literalFix,
                                literalFixRes=.pseudoOptimControl$literalFixRes,
                                scaleTo=0,
                                calcTables=.pseudoOptimControl$calcTables,
                                addProp=.pseudoOptimControl$addProp,
                                #skipCov=.ui$foceiSkipCov,
                                interaction=0L,
                                compress=.pseudoOptimControl$compress,
                                ci=.pseudoOptimControl$ci,
                                sigdigTable=.pseudoOptimControl$sigdigTable)
  if (assign) env$control <- .foceiControl
  .foceiControl
}
#' When the model is loaded, this function is called to evaluate the
#' -2* log-likelihood value for FME::modMCMC()
#'
#' This should not be called directly
#'
#' @param p parameter to evaluate
#' @return -2 * log-likelihood value
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.pseudoOptimF <- function(p) {
  # This returns -2* LL, as required from FME::modMCMC
  2*nlmixr2est::.nlmixrOptimFunC(p)
}

.pseudoOptimFitModel <- function(ui, dataSav) {
  # Use nlmEnv and function for DRY principle
  rxode2::rxReq("FME")
  .ctl <- ui$control
  .p <- setNames(ui$nlmParIni, ui$nlmParName)
  .mi <-  ui$nlmRxModel
  .env <- nlmixr2est::.nlmSetupEnv(.p, ui, dataSav, .mi, .ctl,
                                   lower=ui$optimParLower,
                                   upper=ui$optimParUpper)
  on.exit({nlmixr2est::.nlmFreeEnv()})
  if (is.null(.ctl$npop)) {
    .ctl$npop <- max(5*length(.env$par.ini),50)
  }
  .control <- list(npop=.ctl$npop,
                   numiter=.ctl$numiter,
                   centroid=.ctl$centroid,
                   varleft=.ctl$varleft,
                   verbose=.ctl$verbose)
  .ret <- bquote(FME::pseudoOptim(
    f=.(babelmixr2::.pseudoOptimF),
    p=.(.env$par.ini),
    lower=.(.env$lower),
    upper=.(.env$upper),
    control=.(.control)))
  .ret <- eval(.ret)
  nlmixr2est::.nlmFinalizeList(.env, .ret, par="par", printLine=TRUE,
                               hessianCov=TRUE)
}

#' Get the full theta for nlm methods
#'
#' @param po enhanced pseudoOpt return
#' @param ui ui object
#' @return named theta matrix
#' @author Matthew L. Fidler
#' @noRd
.pseudoOptimGetTheta <- function(po, ui) {
  .iniDf <- ui$iniDf
  stats::setNames(vapply(seq_along(.iniDf$name),
                         function(i) {
                           if (.iniDf$fix[i]) {
                             .iniDf$est[i]
                           } else {
                             po$par[.iniDf$name[i]]
                           }
                         }, double(1), USE.NAMES=FALSE),
                  .iniDf$name)
}

.pseudoOptimFamilyFit <- function(env, ...) {
  .ui <- env$ui
  .control <- .ui$control
  .data <- env$data
  .ret <- new.env(parent=emptyenv())
  # The environment needs:
  # - table for table options
  # - $origData -- Original Data
  # - $dataSav -- Processed data from .foceiPreProcessData
  # - $idLvl -- Level information for ID factor added
  # - $covLvl -- Level information for items to convert to factor
  # - $ui for ui fullTheta Full theta information
  # - $etaObf data frame with ID, etas and OBJI
  # - $cov For covariance
  # - $covMethod for the method of calculating the covariance
  # - $adjObf Should the objective function value be adjusted
  # - $objective objective function value
  # - $extra Extra print information
  # - $method Estimation method (for printing)
  # - $omega Omega matrix
  # - $theta Is a theta data frame
  # - $model a list of model information for table generation.  Needs a `predOnly` model
  # - $message Message for display
  # - $est estimation method
  # - $ofvType (optional) tells the type of ofv is currently being used
  # When running the focei problem to create the nlmixr object, you also need a
  #  foceiControl object
  .ret$table <- env$table
  nlmixr2est::.foceiPreProcessData(.data, .ret, .ui, .control$rxControl)
  .pseudoOptim <- nlmixr2est::.collectWarn(.pseudoOptimFitModel(.ui, .ret$dataSav), lst = TRUE)
  .ret$pseudoOptim <- .pseudoOptim[[1]]
  .ret$parHistData <- .ret$pseudoOptim$parHistData
  .ret$pseudoOptim$parHistData <- NULL
  .ret$message <- .ret$pseudoOptim$message
  if (rxode2::rxGetControl(.ui, "returnpseudoOptim", FALSE)) {
    return(.ret$pseudoOptim)
  }
  .ret$ui <- .ui
  .ret$adjObf <- rxode2::rxGetControl(.ui, "adjObf", TRUE)
  .ret$fullTheta <- .pseudoOptimGetTheta(.ret$pseudoOptim, .ui)
  .ret$cov <- .ret$pseudoOptim$cov
  .ret$covMethod <- .ret$pseudoOptim$covMethod
  .ret$control <- .control
  .ret$extra <- ""
  nlmixr2est::.nlmixr2FitUpdateParams(.ret)
  nlmixr2est::nmObjHandleControlObject(.ret$control, .ret)
  if (exists("control", .ui)) {
    rm(list="control", envir=.ui)
  }
  .ret$est <- "pseudoOptim"
  # There is no parameter history for nlme
  .ret$objective <- as.numeric(.ret$pseudoOptim$cost)
  .ret$model <- .ui$ebe
  .ret$ofvType <- "pseudoOptim"
  .pseudoOptimControlToFoceiControl(.ret)
  .ret$theta <- .ret$ui$saemThetaDataFrame
  .ret <- nlmixr2est::nlmixr2CreateOutputFromUi(.ret$ui, data=.ret$origData, control=.ret$control, table=.ret$table, env=.ret, est="pseudoOptim")
  .env <- .ret$env
  .env$method <- "pseudoOptim"
  .ret
}

nlmixr2Est.pseudoOptim <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiPopulationOnly(.ui, " for the estimation routine 'pseudoOptim', try 'focei'", .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'pseudoOptim'", .var.name=.ui$modelName)
  .ctl <- env$control
  .iniDf <- .ui$iniDf
  if (.ctl$literalFixRes) {
    # Drop fixed parameters; they do not need boundaries
    .iniDf <- .iniDf[!is.na(.iniDf$err) & .iniDf$fix, , drop=FALSE]
  }
  if (.ctl$literalFix) {
    .iniDf <- .iniDf[is.na(.iniDf$err) & .iniDf$fix, , drop=FALSE]
  }
  if (!all(is.finite(.iniDf$lower)) ||
      !all(is.finite(.iniDf$upper))) {
    stop("pseudoOptim requires all parameters to have finite lower and upper bounds",
         call.=FALSE)
  }
  if (.ctl$literalFix)
  .pseudoOptimFamilyControl(env, ...)
  on.exit({if (exists("control", envir=.ui)) rm("control", envir=.ui)}, add=TRUE)
  .pseudoOptimFamilyFit(env,  ...)
}
attr(nlmixr2Est.pseudoOptim, "covPresent") <- TRUE
