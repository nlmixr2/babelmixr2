#' @export
nmObjGetControl.nonmem2rx <- function(x, ...) {
  .env <- x[[1]]
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "foceiControl")) return(.control)
  }
  if (exists("foceiControl0", .env)) {
    .control <- get("foceiControl0", .env)
    if (inherits(.control, "foceiControl")) return(.control)
  }
  stop("cannot find nonmem2rx related control object", call.=FALSE)
}


.nonmem2rxToFoceiControl <- function(env, model, assign=FALSE) {
  .rxControl <- rxode2::rxControl(covsInterpolation="nocb",
                                  atol=model$atol,
                                  rtol=model$rtol,
                                  ssRtol=model$ssRtol,
                                  ssAtol=model$ssAtol,
                                  method="lsoda")
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
as.nlmixr2.nonmem2rx <- function(x, ..., table=nlmixr2est::tableControl()) {
  #need x$nonmemData
  # need x to have at least one endpoint
  # The environment needs:
  env <- new.env(parent=emptyenv())
  .ui <- new.env(parent=emptyenv())
  .oldUi <- rxode2::rxUiDecompress(x)
  for (n in ls(envir=.oldUi, all.names=TRUE)) {
    assign(n, get(n, envir=.oldUi), envir=.ui)
  }
  class(.ui) <- class(.oldUi)
  # - $table for table options -- already present
  env$table <- table
  env$origData <- x$nonmemData
  nlmixr2est::.foceiPreProcessData(env$origData, env, .ui)
  # - $origData -- Original Data -- already present
  # - $dataSav -- Processed data from .foceiPreProcessData --already present
  # - $idLvl -- Level information for ID factor added -- already present
  env$ui <- .ui
  # - $ui for ui fullTheta Full theta information
  env$fullTheta <- .ui$nonmemFullTheta
  # - $etaObf data frame with ID, etas and OBJI
  env$etaObf <- .ui$nonmemEtaObf
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
  .cov <- .ui$nonmemCovariance
  if (!is.null(.cov)) {
    env$cov <- .cov
    # - $covMethod for the method of calculating the covariance
    env$covMethod <- "nonmem2rx"
  }
  # - $objective objective function value
  env$objective <- .ui$nonmemObjf
  # - $extra Extra print information
  env$extra <- paste0(" reading NONMEM ver ", env$ui$nonmemOutputVersion)
  # - $method Estimation method (for printing)
  env$method <- "nonmem2rx"
  # - $omega Omega matrix
  env$omega <- .ui$nonmemOutputOmega
  # - $theta Is a theta data frame
  env$theta <- .ui$nonmemThetaDf
  # - $model a list of model information for table generation.  Needs a `predOnly` model
  env$model <- .ui$ebe
  # - $message Message for display
  env$message <- ""
  # - $est estimation method
  env$est <- "nonmem2rx"
  # - $ofvType (optional) tells the type of ofv is currently being used
  #env$ofvType
  env$ofvType <- .ui$nonmemObjfType
  # Add parameter history
  env$parHist <- .ui$nonmemParHistory
  env$nobs <- x$dfObs
  env$nobs2<- x$dfObs
  # Run before converting to nonmemControl
  .objf <- .ui$nonmemObjf
  # When running the focei problem to create the nlmixr object, you also need a
  #  foceiControl object
  .nonmem2rxToFoceiControl(env, x, TRUE)
  .ret <- nlmixr2est::nlmixr2CreateOutputFromUi(env$ui, data=env$origData,
                                               control=env$control, table=env$table,
                                               env=env, est="nonmem2rx")
  if (inherits(.ret, "nlmixr2FitData")) {
    .msg <- .nonmemMergePredsAndCalcRelativeErr(.ret)
    .prderrPath <- file.path(x$nonmemExportPath, "PRDERR")
    .msg$message <- c(.ui$nonmemTransMessage,
                      .ui$nonmemTermMessage,
                      .msg$message)
    if (file.exists(.prderrPath)) {
      .prderr <- paste(readLines(.prderrPath), collapse="\n")
      .msg$message <- c(.msg$message,
                        "there are solving errors during optimization (see '$prderr')")
      assign("prderr", .prderr, envir=.ret$env)
    }
    .msg$message <- c(.msg$message, paste0("nonmem2rx model file: '", x$file, "'"))
    assign("message", paste(.msg$message, collapse="\n    "), envir=.ret$env)
  }
  .ret
}
