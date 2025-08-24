#' Control for fmeMcmc estimation method in nlmixr2
#'
#' @inheritParams nlmixr2est::foceiControl
#' @inheritParams nlmixr2est::saemControl
#'
#' @inheritParams FME::modMCMC
#'
#' @param returnFmeMcmc return the fmeMcmc output instead of the nlmixr2
#'   fit
#'
#' @return fmeMcmc control structure
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
#'    E0 <- 0.5
#'    Em <- 0.5
#'    E50 <- 2
#'    g <- fix(2)
#'  })
#'  model({
#'    v <- E0+Em*time^g/(E50^g+time^g)
#'    ll(bin) ~ DV * v - log(1 + exp(v))
#'  })
#' }
#'
#' fit2 <- nlmixr(mod, dsn, est="fmeMcmc")
#'
#' print(fit2)
#'
#' # you can also get the FME modMCMC output with
#'
#' # fit2$fmeMcmc
#'
#' # And use it in the summaries from FME, i.e.
#'
#' summary(fit2$fmeMcmc)
#'
#' pairs(fit2$fmeMcmc)
#'
#' # and you can also use the coda package with `as.mcmc()`
#' coda::raftery.diag(coda::as.mcmc(fit2))
#'
#' }
fmeMcmcControl <- function(jump=NULL,
                           prior = NULL,
                           niter = 1000L,
                           outputlength= niter,
                           burninlength = 0,
                           updatecov = niter,
                           covscale = NULL, # 2.4^2/length(p),
                           ntrydr = 1,
                           drscale = NULL,
                           verbose = FALSE,
                           returnFmeMcmc=FALSE,

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
                           covMethod=c("mcmc", "r", ""),
                           adjObf=TRUE, ci=0.95, sigdig=4, sigdigTable=NULL, ...) {


  checkmate::assertLogical(optExpression, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFix, len=1, any.missing=FALSE)
  checkmate::assertLogical(literalFixRes, len=1, any.missing=FALSE)
  checkmate::assertLogical(sumProd, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnFmeMcmc, len=1, any.missing=FALSE)
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
    jump=jump,
    prior = prior,
    niter = as.integer(niter),
    outputlength=as.integer(outputlength),
    burninlength = as.integer(burninlength),
    updatecov = updatecov,
    covscale = covscale,
    ntrydr = as.integer(ntrydr),
    drscale = drscale,
    verbose = as.logical(verbose),
    covMethod=match.arg(covMethod),
    optExpression=optExpression,
    literalFix=literalFix,
    literalFixRes=literalFixRes,
    sumProd=sumProd,
    rxControl=rxControl,
    returnFmeMcmc=returnFmeMcmc,

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
  class(.ret) <- "fmeMcmcControl"
  .ret
}

#' Get the fmeMcmc family control
#'
#' @param env fmeMcmc optimization environment
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.fmeMcmcFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- babelmixr2::fmeMcmcControl()
  }
  if (!inherits(.control, "fmeMcmcControl")){
    .control <- do.call(babelmixr2::fmeMcmcControl, .control)
  }
  assign("control", .control, envir=.ui)
}

nmObjHandleControlObject.fmeMcmcControl <- function(control, env) {
  assign("fmeMcmcControl", control, envir=env)
}

nmObjGetControl.fmeMcmc <- function(x, ...) {
  .env <- x[[1]]
  if (exists("fmeMcmcControl", .env)) {
    .control <- get("fmeMcmcControl", .env)
    if (inherits(.control, "fmeMcmcControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "fmeMcmcControl")) return(.control)
  }
  stop("cannot find fmeMcmc related control object", call.=FALSE)
}

getValidNlmixrCtl.fmeMcmc <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- fmeMcmcControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("fmeMcmcControl", .ctl)
  if (!inherits(.ctl, "fmeMcmcControl")) {
    .minfo("invalid control for `est=\"fmeMcmc\"`, using default")
    .ctl <- fmeMcmcControl()
  } else {
    .ctl <- do.call(fmeMcmcControl, .ctl)
  }
  .ctl
}

.fmeMcmcControlToFoceiControl <- function(env, assign=TRUE) {
  .fmeMcmcControl <- env$fmeMcmcControl
  .ui <- env$ui
  .foceiControl <- foceiControl(rxControl=env$fmeMcmcControl$rxControl,
                                maxOuterIterations=0L,
                                maxInnerIterations=0L,
                                covMethod=0L,
                                sumProd=.fmeMcmcControl$sumProd,
                                optExpression=.fmeMcmcControl$optExpression,
                                literalFix=.fmeMcmcControl$literalFix,
                                literalFixRes=.fmeMcmcControl$literalFixRes,
                                scaleTo=0,
                                calcTables=.fmeMcmcControl$calcTables,
                                addProp=.fmeMcmcControl$addProp,
                                #skipCov=.ui$foceiSkipCov,
                                interaction=0L,
                                compress=.fmeMcmcControl$compress,
                                ci=.fmeMcmcControl$ci,
                                sigdigTable=.fmeMcmcControl$sigdigTable)
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
.fmeMcmcF <- function(p) {
  # This returns -2* LL, as required from FME::modMCMC
  2*nlmixr2est::.nlmixrOptimFunC(p)
}

.fmeMcmcFitModel <- function(ui, dataSav) {
  # Use nlmEnv and function for DRY principle
  rxode2::rxReq("FME")
  .ctl <- ui$control

  ## .keep <- c("npt", "rhobeg", "rhoend", "iprint", "maxfun")
  ## .keep <- .keep[vapply(.keep, function(opt) {
  ##   !is.null(.ctl[[opt]])
  ## }, logical(1), USE.NAMES = FALSE)]

  ## .oCtl <- setNames(lapply(.keep, function(x) {.ctl[[x]]}), .keep)
  ## class(.ctl) <- NULL
  .p <- setNames(ui$nlmParIni, ui$nlmParName)
  .mi <-  ui$nlmRxModel
  .env <- nlmixr2est::.nlmSetupEnv(.p, ui, dataSav, .mi, .ctl,
                                   lower=ui$optimParLower, upper=ui$optimParUpper)
  on.exit({nlmixr2est::.nlmFreeEnv()})
  if (is.null(.ctl$covscale)) {
    .ctl$covscale <- 2.4^2 / length(.env$par.ini)
  }
  .ret <- bquote(FME::modMCMC(
    f=.(babelmixr2::.fmeMcmcF),
    p=.(.env$par.ini),
    lower=.(.env$lower),
    upper=.(.env$upper),
    prior=.(.ctl$prior),
    niter=.(.ctl$niter),
    outputlength=.(.ctl$outputlength),
    burninlength=.(.ctl$burninlength),
    updatecov=.(.ctl$updatecov),
    covscale=.(.ctl$covscale),
    ntrydr=.(.ctl$ntrydr),
    drscale=.(.ctl$drscale),
    verbose=.(.ctl$verbose)))
  .ret <- eval(.ret)
  if (.ctl$covMethod == "mcmc") {
    .ret$cov.unscaled <- cov(.ret$pars)
    .ret$cov <- nlmixr2est::.nlmAdjustCov(.ret$cov.unscaled, .ret$bestpar)
  }
  .ret <- nlmixr2est::.nlmFinalizeList(.env, .ret, par="bestpar", printLine=TRUE,
                                       hessianCov=(.ctl$covMethod == "r"))
  if (.ctl$covMethod == "mcmc") {
    .n <- names(.ret$bestpar)
    dimnames(.ret$cov) <- list(.n, .n)
    dimnames(.ret$cov.unscaled) <- list(.n, .n)
    .ret$covMethod <- "mcmc"
  }
  .ret
}

#' Get the full theta for nlm methods
#'
#' @param fme enhanced fme return
#' @param ui ui object
#' @return named theta matrix
#' @author Matthew L. Fidler
#' @noRd
.fmeMcmcGetTheta <- function(fme, ui) {
  .iniDf <- ui$iniDf
  stats::setNames(vapply(seq_along(.iniDf$name),
                  function(i) {
                    if (.iniDf$fix[i]) {
                      .iniDf$est[i]
                    } else {
                      fme$bestpar[.iniDf$name[i]]
                    }
                  }, double(1), USE.NAMES=FALSE),
           .iniDf$name)
}

.fmeMcmcFamilyFit <- function(env, ...) {
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
  .fmeMcmc <- nlmixr2est::.collectWarn(.fmeMcmcFitModel(.ui, .ret$dataSav), lst = TRUE)
  .ret$fmeMcmc <- .fmeMcmc[[1]]
  .ret$parHistData <- .ret$fmeMcmc$parHistData
  .ret$fmeMcmc$parHistData <- NULL
  .ret$message <- .ret$fmeMcmc$message
  if (rxode2::rxGetControl(.ui, "returnFmeMcmc", FALSE)) {
    return(.ret$fmeMcmc)
  }
  .ret$ui <- .ui
  .ret$adjObf <- rxode2::rxGetControl(.ui, "adjObf", TRUE)
  .ret$fullTheta <- .fmeMcmcGetTheta(.ret$fmeMcmc, .ui)
  .ret$cov <- .ret$fmeMcmc$cov
  .ret$covMethod <- .ret$fmeMcmc$covMethod
  #.ret$etaMat <- NULL
  #.ret$etaObf <- NULL
  #.ret$omega <- NULL
  .ret$control <- .control
  .ret$extra <- ""
  nlmixr2est::.nlmixr2FitUpdateParams(.ret)
  nlmixr2est::nmObjHandleControlObject(.ret$control, .ret)
  if (exists("control", .ui)) {
    rm(list="control", envir=.ui)
  }
  .ret$est <- "fmeMcmc"
  # There is no parameter history for nlme
  .ret$objective <- as.numeric(.ret$fmeMcmc$bestfunp)
  .ret$model <- .ui$ebe
  .ret$ofvType <- "fmeMcmc"
  .fmeMcmcControlToFoceiControl(.ret)
  .ret$theta <- .ret$ui$saemThetaDataFrame
  .ret <- nlmixr2est::nlmixr2CreateOutputFromUi(.ret$ui, data=.ret$origData, control=.ret$control, table=.ret$table, env=.ret, est="fmeMcmc")
  .env <- .ret$env
  .env$method <- "fmeMcmc"
  .ret
}

nlmixr2Est.fmeMcmc <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiPopulationOnly(.ui, " for the estimation routine 'fmeMcmc', try 'focei'", .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'fmeMcmc'", .var.name=.ui$modelName)
  .fmeMcmcFamilyControl(env, ...)
  on.exit({if (exists("control", envir=.ui)) rm("control", envir=.ui)}, add=TRUE)
  .fmeMcmcFamilyFit(env,  ...)
}
attr(nlmixr2Est.fmeMcmc, "covPresent") <- TRUE

as.mcmc.nlmixr2.fmeMcmc <- function(x, ...) {
  if (!inherits(x, "nlmixr2.fmeMcmc")) {
    stop("x must be a nlmixr2.fmeMcmc object", call.=FALSE)
  }
  coda::as.mcmc(x$fmeMcmc, ...)
}
