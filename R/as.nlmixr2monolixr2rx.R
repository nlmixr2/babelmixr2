#' @export
nmObjGetControl.monolix2rx <- function(x, ...) {
  .env <- x[[1]]
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "foceiControl")) return(.control)
  }
  if (exists("foceiControl0", .env)) {
    .control <- get("foceiControl0", .env)
    if (inherits(.control, "foceiControl")) return(.control)
  }
  stop("cannot find monolix2rx related control object", call.=FALSE)
}

.monolix2rxToFoceiControl <- function(env, model, assign=FALSE) {
  ## maxSS=nbSSDoses + 1,
  ## minSS=nbSSDoses,
  ## ssAtol=100,
  ## ssRtol=100,
  ## atol=ifelse(stiff, 1e-9, 1e-6),
  ## rtol=ifelse(stiff, 1e-6, 1e-3),
  ## method=ifelse(stiff, "liblsoda", "dop853")
  .nbSsDoses <- monolix2rx::.getNbdoses(model)
  .stiff <- monolix2rx::.getStiff(model)
  .rxControl <- rxode2::rxControl(covsInterpolation="locf",
                                  atol=ifelse(.stiff, 1e-9, 1e-6),
                                  rtol=ifelse(.stiff, 1e-6, 1e-3),
                                  ssRtol=100,
                                  ssAtol=100,
                                  maxSS=.nbSsDoses + 1,
                                  minSS=.nbSsDoses,
                                  method=ifelse(.stiff, "liblsoda", "dop853"),
                                  safeZero=FALSE)
  .foceiControl <- nlmixr2est::foceiControl(rxControl=.rxControl,
                                            maxOuterIterations = 0L, maxInnerIterations = 0L,
                                            etaMat = env$etaMat,
                                            covMethod=0L,
                                            interaction = 1L)
  if (assign)
    env$control <- .foceiControl
  .foceiControl
}

#' @export
as.nlmixr2.monolix2rx <- function(x, ..., table=nlmixr2est::tableControl(), rxControl=rxode2::rxControl()) {
  #need x$nonmemData
  # need x to have at least one endpoint
  # The environment needs:
  env <- new.env(parent=emptyenv())
  x <- rxode2::rxUiDecompress(x)
  nlmixr2est::nlmixrWithTiming("as.nlmixr2", {
    .ui <- new.env(parent=emptyenv())
    .oldUi <- x
    for (n in ls(envir=.oldUi, all.names=TRUE)) {
      assign(n, get(n, envir=.oldUi), envir=.ui)
    }
    class(.ui) <- class(.oldUi)
    # - $table for table options -- already present
    env$table <- table
    env$origData <- x$monolixData
    nlmixr2est::.foceiPreProcessData(env$origData, env, .ui, rxControl)
    # - $origData -- Original Data -- already present
    # - $dataSav -- Processed data from .foceiPreProcessData --already present
    # - $idLvl -- Level information for ID factor added -- already present
    env$ui <- .ui
    # - $ui for ui fullTheta Full theta information
    env$fullTheta <- .ui$monolixFullTheta
    # - $etaObf data frame with ID, etas and OBJI
    env$etaObf <- .ui$monolixEtaObf
    if (is.null(env$etaObf)) {
      .df <- data.frame(ID=unique(env$dataSav$ID))
      for (.n in .getEtaNames(.ui)) {
        .df[[.n]] <- 0
      }
      .df[["OBJI"]] <- NA_real_
      env$etaObf <- .df
      warning("since NONMEM did not output between subject variability, assuming all ETA(#) are zero",
              call.=FALSE)
    }
    # - $cov For covariance
    .cov <- .ui$monolixCovariance
    if (!is.null(.cov)) {
      env$cov <- .cov
      # - $covMethod for the method of calculating the covariance
      env$covMethod <- "monolix2rx"
    }
    # - $objective objective function value
    env$objective <- .ui$monolixObjf
    # - $extra Extra print information
    env$extra <- paste0(" reading Monolix ver ", env$ui$monolixOutputVersion)
    # - $method Estimation method (for printing)
    env$method <- "monolix2rx"
    # - $omega Omega matrix
    env$omega <- .ui$monolixOmega
    # - $theta Is a theta data frame
    env$theta <- .ui$monolixTheta
    # - $model a list of model information for table generation.  Needs a `predOnly` model
    env$model <- .ui$ebe
    # - $message Message for display
    env$message <- ""
    # - $est estimation method
    env$est <- "monolix2rx"
    # - $ofvType (optional) tells the type of ofv is currently being used
    #env$ofvType
    env$ofvType <- .ui$monolixObjfType
    # Add parameter history
    env$nobs <- x$dfObs
    env$nobs2<- x$dfObs
    # Run before converting to nonmemControl
    .objf <- .ui$monolixObjf
    # When running the focei problem to create the nlmixr object, you also need a
    #  foceiControl object
    .monolix2rxToFoceiControl(env, x, TRUE)
    .ret <- nlmixr2est::nlmixr2CreateOutputFromUi(env$ui, data=env$origData,
                                                  control=env$control, table=env$table,
                                                  env=env, est="monolix2rx")
    if (inherits(.ret, "nlmixr2FitData")) {
      .msg <- .monolixMergePredsAndCalcRelativeErr(.ret)
      .msg$message <- c(.msg$message)
      .tmp <- .ret$ui$monolixParHistory
      assign("message", paste(.msg$message, collapse="\n    "), envir=.ret$env)
      if (is.null(.tmp)) {
        .minfo("monolix parameter history needs exported charts, please export charts")
      } else {
        .tmp$type <- "Unscaled"
        assign("parHistData", .tmp, .ret$env)
        .minfo("monolix parameter history integrated into fit object")
      }
    }
    ## .time <- get("time", .ret$env)
    ## .time <- .time[,!(names(.time) %in% c("optimize", "covariance"))]
    ## assign("time",
    ##        cbind(.time, data.frame(NONMEM=.ui$nonmemRunTime)),
    ##        .ret$env)
    .ret
  }, env=env)
}
