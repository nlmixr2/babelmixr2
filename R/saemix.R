#' @export
nlmixr2Est.saemix <- function(env, ...) {
  .model <- nlmixr2est::.uiApplyMu2(env)
  .ui <- env$ui

  # Decompress UI object if compressed
  .ui <- rxode2::rxUiDecompress(env$ui)
  nlmixr2est::nmObjUiSetCompressed(FALSE)
  on.exit({nlmixr2est::nmObjUiSetCompressed(TRUE)})
  assign("ui", .ui, envir = env)
  on.exit({
    assign("ui", rxode2::rxUiCompress(env$ui), envir = env)
  }, add = TRUE)

  # Assertions
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'saemix'", .var.name = .ui$modelName)

  # Setup control
  .saemixFamilyControl(env, ...)
  on.exit({
    if (exists("control", envir = .ui)) {
      rm("control", envir = .ui)
    }
  }, add = TRUE)

  nlmixr2est::.uiFinalizeMu2(.saemixFamilyFit(env, ...), .model)
}

.saemixFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- saemixControl()
  }
  if (!inherits(.control, "saemixControl")) {
    .control <- do.call(saemixControl, .control)
  }
  assign("control", .control, envir = .ui)
}

.saemixControlToFoceiControl <- function(env, assign = TRUE) {
  .saemixControl <- env$saemixControl
  .ui <- env$ui
  .foceiControl <- nlmixr2est::foceiControl(rxControl = env$saemixControl$rxControl,
                                            maxOuterIterations = 0L,
                                            maxInnerIterations = 0L,
                                            covMethod = 0L,
                                            calcTables = .saemixControl$calcTables,
                                            compress = .saemixControl$compress,
                                            ci = .saemixControl$ci,
                                            sigdigTable = .saemixControl$sigdigTable)
  if (assign) env$control <- .foceiControl
  .foceiControl
}

.saemixFamilyFit <- function(env, ...) {
  .ui <- env$ui
  .control <- .ui$control
  .data <- env$data

  .ret <- new.env(parent = emptyenv())
  .ret$table <- env$table
  .ret$saemixControl <- .control

  # Pre-process data
  nlmixr2est::.foceiPreProcessData(.data, .ret, .ui, .control$rxControl)
  # dataSav has ID, TIME, DV, AMT, EVID, CMT, II
  # Add origRow tracking column
  .ret$dataSav$origRow <- seq_len(nrow(.ret$dataSav))
  # Add mdv column (dosing rows have EVID != 0 or missing DV)
  .ret$dataSav$mdv <- ifelse(is.na(.ret$dataSav$DV), 1, 0)

  # Convert logical columns to numeric to satisfy saemixData validations
  for (.col in colnames(.ret$dataSav)) {
    if (is.logical(.ret$dataSav[[.col]])) {
      .ret$dataSav[[.col]] <- as.numeric(.ret$dataSav[[.col]])
    }
  }

  # Precompute row indices and counts per ID for fast rebuilding of evExpanded inside .saemixModelFunction
  .ret$rowsById <- split(seq_len(nrow(.ret$dataSav)), .ret$dataSav$ID)
  .ret$subjectRowCount <- sapply(.ret$rowsById, length)

  # Check if likelihood model
  .isLikelihood <- any(.ui$predDf$distribution == "LL")
  .modelType <- if (.isLikelihood) "likelihood" else "structural"

  # Extract structural parameters (thetas)
  # These are thetas not associated with residual error
  .structuralThetas <- .ui$iniDf[is.na(.ui$iniDf$neta1) & is.na(.ui$iniDf$err), ]
  .thetaNames <- .structuralThetas$name

  # Initial estimates (psi0)
  .psi0Vec <- .structuralThetas$est
  names(.psi0Vec) <- .thetaNames

  # Fixed parameter indicator (1=estimated, 0=fixed)
  .fixedEstim <- ifelse(.structuralThetas$fix, 0, 1)

  # Check which parameters have etas (mixed-effect parameters)
  .hasEta <- .thetaNames %in% .ui$muRefTable$theta

  # Covariance model (matrix of 1/0 for estimation/correlation of random effects)
  .covMatrix <- matrix(0, nrow = length(.thetaNames), ncol = length(.thetaNames))
  dimnames(.covMatrix) <- list(.thetaNames, .thetaNames)
  .omegaInit <- matrix(0, nrow = length(.thetaNames), ncol = length(.thetaNames))
  dimnames(.omegaInit) <- list(.thetaNames, .thetaNames)

  for (i in seq_along(.thetaNames)) {
    pI <- .thetaNames[i]
    if (.hasEta[i]) {
      etaI <- .ui$muRefTable$eta[.ui$muRefTable$theta == pI]
      for (j in seq_along(.thetaNames)) {
        pJ <- .thetaNames[j]
        if (.hasEta[j]) {
          etaJ <- .ui$muRefTable$eta[.ui$muRefTable$theta == pJ]
          # Look up covariance in ui$omega
          val <- .ui$omega[etaI, etaJ]
          .omegaInit[i, j] <- val
          if (val != 0) {
            .covMatrix[i, j] <- 1
          }
        }
      }
    }
  }

  # Error model setup (only for structural models)
  .errorModel <- NULL
  .errorInit <- NULL
  if (.modelType == "structural") {
    .errType <- .ui$predDf$errType[1]
    if (.errType == "add") {
      .errorModel <- "constant"
      .errRows <- .ui$iniDf[!is.na(.ui$iniDf$err) & .ui$iniDf$err == "add", ]
      .errorInit <- c(.errRows$est[1], 0)
    } else if (.errType == "prop") {
      .errorModel <- "proportional"
      .errRows <- .ui$iniDf[!is.na(.ui$iniDf$err) & .ui$iniDf$err == "prop", ]
      .errorInit <- c(0, .errRows$est[1])
    } else if (.errType %in% c("add + prop", "combined", "combined2", "combined1")) {
      .errorModel <- "combined"
      .addRows <- .ui$iniDf[!is.na(.ui$iniDf$err) & .ui$iniDf$err == "add", ]
      .propRows <- .ui$iniDf[!is.na(.ui$iniDf$err) & .ui$iniDf$err == "prop", ]
      .errorInit <- c(.addRows$est[1], .propRows$est[1])
    } else {
      .errorModel <- "constant"
      .errorInit <- c(1, 0)
    }
  }

  .saemixPredictors <- c("origRow", "TIME")
  .etaNamesUi <- .ui$muRefTable$eta
  if (is.null(.etaNamesUi)) .etaNamesUi <- character(0)

  .saemixModelFunction <- function(psi, id, xidep) {
    .oldSeed <- rxode2::.rxGetSeed()
    on.exit(rxode2::.rxSetSeed(.oldSeed))
    tryCatch({
      if (is.null(dim(psi))) {
        psi <- matrix(psi, nrow = 1)
        colnames(psi) <- .thetaNames
      }
      M <- nrow(psi)

      # Reconstruct original row indices and IDs
      origRowCol <- xidep[, 1]
      simOrigIds <- .ret$dataSav$ID[origRowCol[match(1:M, id)]]

      # Check if cached version of evExpanded exists and is identical
      hit <- FALSE
      if (!is.null(.ret$cachedSimOrigIds) &&
          length(.ret$cachedSimOrigIds) == length(simOrigIds) &&
          all(.ret$cachedSimOrigIds == simOrigIds)) {
        evExpanded <- .ret$cachedEvExpanded
        hit <- TRUE
      } else {
        # Rebuild evExpanded using precomputed indices
        indices <- unlist(.ret$rowsById[as.character(simOrigIds)], use.names = FALSE)
        evExpanded <- .ret$dataSav[indices, ]

        repCounts <- .ret$subjectRowCount[as.character(simOrigIds)]
        evExpanded$ID <- rep(1:M, times = repCounts)
        evExpanded$sim_id <- evExpanded$ID

        # Save to cache
        .ret$cachedSimOrigIds <- simOrigIds
        .ret$cachedEvExpanded <- evExpanded
      }

      # Construct params template data frame
      if (!is.null(.ret$cachedM) && .ret$cachedM == M) {
        params <- .ret$cachedParamsTemplate
      } else {
        params <- as.data.frame(matrix(0, nrow = M, ncol = length(.thetaNames) + length(.etaNamesUi)))
        colnames(params) <- c(.thetaNames, .etaNamesUi)
        params$ID <- 1:M
        .ret$cachedM <- M
        .ret$cachedParamsTemplate <- params
      }
      params[, .thetaNames] <- psi

      # Run the ODE solver
      res <- rxode2::rxSolve(.ui, events = evExpanded, params = params, keep = "origRow", addDosing = TRUE)

      # Get or compute matchedIdx
      if (hit && !is.null(.ret$cachedMatchedIdx)) {
        matchedIdx <- .ret$cachedMatchedIdx
      } else {
        # Match the rows
        resId <- if ("id" %in% names(res)) as.integer(res$id) else rep(1L, length(res$time))
        keyRes <- paste(resId, res$origRow, sep = "_")
        keyXidep <- paste(id, origRowCol, sep = "_")
        matchedIdx <- match(keyXidep, keyRes)

        # Save to cache
        .ret$cachedMatchedIdx <- matchedIdx
      }

      predCols <- .ui$predDf$cond
      if (.isLikelihood) {
        predCols <- rep("pred", length(predCols))
      }

      if (length(predCols) == 1) {
        predictions <- res[[predCols]][matchedIdx]
      } else {
        # Locate ytype column in .ret$dataSav$origRow mapping
        ytypeCol <- if ("YTYPE" %in% colnames(.ret$dataSav)) "YTYPE" else "ytype"
        ytypeVec <- .ret$dataSav[[ytypeCol]][origRowCol]

        predictions <- numeric(nrow(xidep))
        for (yVal in unique(ytypeVec)) {
          rows <- which(ytypeVec == yVal)
          predCol <- predCols[yVal]
          predictions[rows] <- res[[predCol]][matchedIdx[rows]]
        }
      }

      return(predictions)
    }, error = function(e) {
      cat("ERROR IN .saemixModelFunction:\n")
      print(e)
      stop(e)
    })
  }


  # Build SaemixData
  .dataForSaemix <- .ret$dataSav[.ret$dataSav$mdv == 0, ]
  .saemixCovariates <- setdiff(.ui$allCovs, "DV")
  if (length(.saemixCovariates) == 0) {
    .saemixCovariates <- NULL
  }
  if (is.null(.saemixCovariates)) {
    saemixData <- saemix::saemixData(name.data = .dataForSaemix,
                                     name.group = "ID",
                                     name.predictors = .saemixPredictors,
                                     name.response = "DV",
                                     name.X = "TIME",
                                     verbose = .control$warnings)
  } else {
    saemixData <- saemix::saemixData(name.data = .dataForSaemix,
                                     name.group = "ID",
                                     name.predictors = .saemixPredictors,
                                     name.response = "DV",
                                     name.covariates = .saemixCovariates,
                                     name.X = "TIME",
                                     verbose = .control$warnings)
  }

  # Build SaemixModel
  saemixModel <- saemix::saemixModel(model = .saemixModelFunction,
                                     modeltype = .modelType,
                                     psi0 = .psi0Vec,
                                     transform.par = rep(0, length(.thetaNames)),
                                     fixed.estim = .fixedEstim,
                                     covariance.model = .covMatrix,
                                     omega.init = .omegaInit,
                                     error.model = .errorModel,
                                     error.init = .errorInit,
                                     verbose = .control$warnings)

  # Build SaemixOptions
  saemixOptions <- list(map = .control$map,
                         fim = .control$fim,
                         ll.is = .control$ll.is,
                         ll.gq = .control$ll.gq,
                         nbiter.saemix = .control$nbiter.saemix,
                         nbiter.burn = .control$nbiter.burn,
                         nbiter.map = .control$nbiter.map,
                         nb.chains = .control$nb.chains,
                         fix.seed = .control$fix.seed,
                         seed = .control$seed,
                         nmc.is = .control$nmc.is,
                         nu.is = .control$nu.is,
                         print.is = .control$print.is,
                         nbdisplay = .control$nbdisplay,
                         displayProgress = .control$displayProgress,
                         print = .control$iterPrintControl$every > 0L,
                         save = FALSE, # handled by nlmixr2
                         save.graphs = FALSE, # handled by nlmixr2
                         warnings = .control$warnings,
                         nbiter.mcmc = .control$nbiter.mcmc,
                         proba.mcmc = .control$proba.mcmc,
                         stepsize.rw = .control$stepsize.rw,
                         rw.init = .control$rw.init,
                         alpha.sa = .control$alpha.sa,
                         nnodes.gq = .control$nnodes.gq,
                         nsd.gq = .control$nsd.gq,
                         maxim.maxiter = .control$maxim.maxiter,
                         nb.sim = .control$nb.sim,
                         nb.simpred = .control$nb.simpred,
                         ipar.lmcmc = .control$ipar.lmcmc,
                         ipar.rmcmc = .control$ipar.rmcmc)

  if (!is.na(.control$nbiter.sa)) {
    saemixOptions$nbiter.sa <- .control$nbiter.sa
  }

  # Run saemix
  fit <- withCallingHandlers({
    .fit <- NULL
    utils::capture.output({
      .fit <- saemix::saemix(saemixModel, saemixData, saemixOptions)
    })
    .fit
  }, error = function(e) {
    cat("\n=== INNER ERROR IN SAEMIX ===\n")
    print(e)
    calls <- sys.calls()
    for (i in seq_along(calls)) {
      cat(sprintf("[%2d] %s\n", i, paste(deparse(calls[[i]]), collapse = "\n")))
    }
    stop(e)
  })

  # Convert results back
  .ret$ui <- .ui
  .ret$saemix <- fit

  # 1. fullTheta
  # Initialize with all thetas from UI
  .fullTheta <- stats::setNames(.ui$iniDf$est[is.na(.ui$iniDf$neta1)], .ui$iniDf$name[is.na(.ui$iniDf$neta1)])
  # Update estimated/fixed structural parameters
  for (i in seq_along(.thetaNames)) {
    pName <- .thetaNames[i]
    .fullTheta[pName] <- fit@results@fixed.effects[i]
  }
  # Update residual parameters
  if (.modelType == "structural") {
    .errType <- .ui$predDf$errType[1]
    if (.errType == "add") {
      errName <- .ui$iniDf$name[!is.na(.ui$iniDf$err) & .ui$iniDf$err == "add"]
      .fullTheta[errName] <- fit@results@respar[1]
    } else if (.errType == "prop") {
      errName <- .ui$iniDf$name[!is.na(.ui$iniDf$err) & .ui$iniDf$err == "prop"]
      .fullTheta[errName] <- fit@results@respar[2]
    } else if (.errType %in% c("add + prop", "combined", "combined2", "combined1")) {
      addName <- .ui$iniDf$name[!is.na(.ui$iniDf$err) & .ui$iniDf$err == "add"]
      propName <- .ui$iniDf$name[!is.na(.ui$iniDf$err) & .ui$iniDf$err == "prop"]
      .fullTheta[addName] <- fit@results@respar[1]
      .fullTheta[propName] <- fit@results@respar[2]
    }
  }
  .ret$fullTheta <- .fullTheta

  # 2. etaObf
  etaNamesFit <- colnames(fit@results@map.eta)
  etaNamesUi <- sapply(etaNamesFit, function(name) {
    thetaName <- sub("^eta\\.", "", name)
    .ui$muRefTable$eta[.ui$muRefTable$theta == thetaName]
  })
  etaObf <- as.data.frame(fit@results@map.eta)
  colnames(etaObf) <- etaNamesUi
  etaObf$ID <- unique(.ret$dataSav$ID)
  etaObf <- etaObf[, c("ID", etaNamesUi), drop = FALSE]
  etaObf$OBJI <- NA_real_
  .ret$etaObf <- etaObf

  # 3. omega
  # Extract omega matrix corresponding to UI eta names
  .etaNamesUiAll <- .ui$iniDf$name[!is.na(.ui$iniDf$neta1)]
  .omega <- matrix(0, nrow = length(.etaNamesUiAll), ncol = length(.etaNamesUiAll))
  dimnames(.omega) <- list(.etaNamesUiAll, .etaNamesUiAll)
  for (i in seq_along(.etaNamesUiAll)) {
    etaI <- .etaNamesUiAll[i]
    thetaI <- .ui$muRefTable$theta[.ui$muRefTable$eta == etaI]
    if (length(thetaI) == 1) {
      for (j in seq_along(.etaNamesUiAll)) {
        etaJ <- .etaNamesUiAll[j]
        thetaJ <- .ui$muRefTable$theta[.ui$muRefTable$eta == etaJ]
        if (length(thetaJ) == 1) {
          .omega[etaI, etaJ] <- fit@results@omega[thetaI, thetaJ]
        }
      }
    }
  }
  .ret$omega <- .omega

  # 4. covariance matrix
  estThetaNames <- .thetaNames[!.structuralThetas$fix]
  covMatrix <- fit@results@MCOV
  if (!is.null(covMatrix) && nrow(covMatrix) == length(estThetaNames)) {
    dimnames(covMatrix) <- list(estThetaNames, estThetaNames)
    .ret$cov <- covMatrix
  }

  # 5. objective
  ll <- fit@results@ll.is
  if (is.null(ll) || is.na(ll)) {
    ll <- fit@results@LL
  }
  .ret$objective <- if (!is.null(ll) && !is.na(ll)) -2 * ll else NA_real_

  # 6. Metadata
  .ret$method <- "saemix"
  .ret$est <- "saemix"
  .ret$extra <- ""
  .ret$message <- ""
  .ret$model <- .ui$ebe
  .ret$ofvType <- if (.isLikelihood) "saemixLik" else "saemix"

  # Update fit parameters
  nlmixr2est::.nlmixr2FitUpdateParams(.ret)
  nmObjHandleControlObject(.ret$saemixControl, .ret)

  if (exists("control", .ui)) {
    rm(list = "control", envir = .ui)
  }

  # Generate output
  .saemixControlToFoceiControl(.ret)
  .ret$theta <- .ret$ui$saemThetaDataFrame
  .ret <- nlmixr2est::nlmixr2CreateOutputFromUi(.ret$ui, data = .ret$origData, control = .ret$control, table = .ret$table, env = .ret, est = "saemix")

  .env <- .ret$env
  .env$method <- "saemix"
  # Embed saemix object into the final fit environment
  assign("saemix", fit, envir = .env)

  .ret
}
