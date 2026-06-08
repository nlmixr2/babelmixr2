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
  # Add orig_row tracking column
  .ret$dataSav$orig_row <- seq_len(nrow(.ret$dataSav))
  # Add mdv column (dosing rows have EVID != 0 or missing DV)
  .ret$dataSav$mdv <- ifelse(is.na(.ret$dataSav$DV), 1, 0)

  # Convert logical columns to numeric to satisfy saemixData validations
  for (.col in colnames(.ret$dataSav)) {
    if (is.logical(.ret$dataSav[[.col]])) {
      .ret$dataSav[[.col]] <- as.numeric(.ret$dataSav[[.col]])
    }
  }

  # Check if likelihood model
  .isLikelihood <- any(.ui$predDf$distribution == "LL")
  .modeltype <- if (.isLikelihood) "likelihood" else "structural"

  # Extract structural parameters (thetas)
  # These are thetas not associated with residual error
  .structural_thetas <- .ui$iniDf[is.na(.ui$iniDf$neta1) & is.na(.ui$iniDf$err), ]
  .theta_names <- .structural_thetas$name

  # Initial estimates (psi0)
  .psi0_vec <- .structural_thetas$est
  names(.psi0_vec) <- .theta_names

  # Fixed parameter indicator (1=estimated, 0=fixed)
  .fixed.estim <- ifelse(.structural_thetas$fix, 0, 1)

  # Check which parameters have etas (mixed-effect parameters)
  .has_eta <- .theta_names %in% .ui$muRefTable$theta

  # Covariance model (matrix of 1/0 for estimation/correlation of random effects)
  .cov_matrix <- matrix(0, nrow = length(.theta_names), ncol = length(.theta_names))
  dimnames(.cov_matrix) <- list(.theta_names, .theta_names)
  .omega_init <- matrix(0, nrow = length(.theta_names), ncol = length(.theta_names))
  dimnames(.omega_init) <- list(.theta_names, .theta_names)

  for (i in seq_along(.theta_names)) {
    p_i <- .theta_names[i]
    if (.has_eta[i]) {
      eta_i <- .ui$muRefTable$eta[.ui$muRefTable$theta == p_i]
      for (j in seq_along(.theta_names)) {
        p_j <- .theta_names[j]
        if (.has_eta[j]) {
          eta_j <- .ui$muRefTable$eta[.ui$muRefTable$theta == p_j]
          # Look up covariance in ui$omega
          val <- .ui$omega[eta_i, eta_j]
          .omega_init[i, j] <- val
          if (val != 0) {
            .cov_matrix[i, j] <- 1
          }
        }
      }
    }
  }

  # Error model setup (only for structural models)
  .error.model <- NULL
  .error.init <- NULL
  if (.modeltype == "structural") {
    .errType <- .ui$predDf$errType[1]
    if (.errType == "add") {
      .error.model <- "constant"
      .err_rows <- .ui$iniDf[!is.na(.ui$iniDf$err) & .ui$iniDf$err == "add", ]
      .error.init <- c(.err_rows$est[1], 0)
    } else if (.errType == "prop") {
      .error.model <- "proportional"
      .err_rows <- .ui$iniDf[!is.na(.ui$iniDf$err) & .ui$iniDf$err == "prop", ]
      .error.init <- c(0, .err_rows$est[1])
    } else if (.errType %in% c("add + prop", "combined", "combined2", "combined1")) {
      .error.model <- "combined"
      .add_rows <- .ui$iniDf[!is.na(.ui$iniDf$err) & .ui$iniDf$err == "add", ]
      .prop_rows <- .ui$iniDf[!is.na(.ui$iniDf$err) & .ui$iniDf$err == "prop", ]
      .error.init <- c(.add_rows$est[1], .prop_rows$est[1])
    } else {
      .error.model <- "constant"
      .error.init <- c(1, 0)
    }
  }

  # Build predictors list
  .saemix_predictors <- c("orig_row", "TIME")

  .saemixModelFunction <- function(psi, id, xidep) {

    if (is.null(dim(psi))) {
      psi <- matrix(psi, nrow = 1)
      colnames(psi) <- .theta_names
    }
    M <- nrow(psi)

    # Reconstruct original row indices and IDs
    orig_row_col <- xidep[, 1]

    sim_orig_ids <- numeric(M)
    for (i in 1:M) {
      obs_idx <- which(id == i)[1]
      sim_orig_ids[i] <- .ret$dataSav$ID[orig_row_col[obs_idx]]
    }

    # Pre-split dataSav by ID
    data_by_id <- split(.ret$dataSav, .ret$dataSav$ID)

    # Rebuild ev_expanded
    ev_list <- lapply(1:M, function(i) {
      df_subj <- data_by_id[[as.character(sim_orig_ids[i])]]
      df_subj$ID <- i
      df_subj$sim_id <- i
      df_subj
    })
    ev_expanded <- do.call(rbind, ev_list)

    # Construct params for rxSolve
    params <- as.data.frame(psi)
    colnames(params) <- .theta_names
    params$ID <- 1:M

    # Set all model etas to 0
    .eta_names_ui <- .ui$muRefTable$eta
    for (eta_name in .eta_names_ui) {
      params[[eta_name]] <- 0
    }

    # Run the ODE solver
    res <- rxode2::rxSolve(.ui, events = ev_expanded, params = params, keep = "orig_row")

    # Match the rows
    res_id <- if ("id" %in% names(res)) as.integer(res$id) else rep(1L, length(res$time))
    key_res <- paste(res_id, res$orig_row, sep = "_")
    key_xidep <- paste(id, orig_row_col, sep = "_")
    matched_idx <- match(key_xidep, key_res)

    pred_cols <- .ui$predDf$cond
    if (.isLikelihood) {
      pred_cols <- rep("pred", length(pred_cols))
    }

    if (length(pred_cols) == 1) {
      predictions <- res[[pred_cols]][matched_idx]
    } else {
      # Locate ytype column in .ret$dataSav using orig_row_col mapping
      ytype_col <- if ("YTYPE" %in% colnames(.ret$dataSav)) "YTYPE" else "ytype"
      ytype_vec <- .ret$dataSav[[ytype_col]][orig_row_col]

      predictions <- numeric(nrow(xidep))
      for (y_val in unique(ytype_vec)) {
        rows <- which(ytype_vec == y_val)
        pred_col <- pred_cols[y_val]
        predictions[rows] <- res[[pred_col]][matched_idx[rows]]
      }
    }
    return(predictions)
  }

  # Build SaemixData
  .saemix_covariates <- .ui$allCovs
  if (length(.saemix_covariates) == 0) {
    .saemix_covariates <- NULL
  }
  if (is.null(.saemix_covariates)) {
    saemix.data <- saemix::saemixData(name.data = .ret$dataSav,
                                      name.group = "ID",
                                      name.predictors = .saemix_predictors,
                                      name.response = "DV",
                                      name.mdv = "mdv",
                                      name.X = "TIME",
                                      verbose = .control$warnings)
  } else {
    saemix.data <- saemix::saemixData(name.data = .ret$dataSav,
                                      name.group = "ID",
                                      name.predictors = .saemix_predictors,
                                      name.response = "DV",
                                      name.mdv = "mdv",
                                      name.covariates = .saemix_covariates,
                                      name.X = "TIME",
                                      verbose = .control$warnings)
  }

  # Build SaemixModel
  saemix.model <- saemix::saemixModel(model = .saemixModelFunction,
                                      modeltype = .modeltype,
                                      psi0 = .psi0_vec,
                                      transform.par = rep(0, length(.theta_names)),
                                      fixed.estim = .fixed.estim,
                                      covariance.model = .cov_matrix,
                                      omega.init = .omega_init,
                                      error.model = .error.model,
                                      error.init = .error.init,
                                      verbose = .control$warnings)

  # Build SaemixOptions
  saemix.options <- list(map = .control$map,
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
                         print = .control$print,
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
    saemix.options$nbiter.sa <- .control$nbiter.sa
  }

  # Run saemix
  fit <- saemix::saemix(saemix.model, saemix.data, saemix.options)

  # Convert results back
  .ret$ui <- .ui
  .ret$saemix <- fit

  # 1. fullTheta
  # Initialize with all thetas from UI
  .fullTheta <- stats::setNames(.ui$iniDf$est[is.na(.ui$iniDf$neta1)], .ui$iniDf$name[is.na(.ui$iniDf$neta1)])
  # Update estimated/fixed structural parameters
  for (i in seq_along(.theta_names)) {
    p_name <- .theta_names[i]
    .fullTheta[p_name] <- fit@results@fixed.effects[i]
  }
  # Update residual parameters
  if (.modeltype == "structural") {
    .errType <- .ui$predDf$errType[1]
    if (.errType == "add") {
      err_name <- .ui$iniDf$name[!is.na(.ui$iniDf$err) & .ui$iniDf$err == "add"]
      .fullTheta[err_name] <- fit@results@respar[1]
    } else if (.errType == "prop") {
      err_name <- .ui$iniDf$name[!is.na(.ui$iniDf$err) & .ui$iniDf$err == "prop"]
      .fullTheta[err_name] <- fit@results@respar[2]
    } else if (.errType %in% c("add + prop", "combined", "combined2", "combined1")) {
      add_name <- .ui$iniDf$name[!is.na(.ui$iniDf$err) & .ui$iniDf$err == "add"]
      prop_name <- .ui$iniDf$name[!is.na(.ui$iniDf$err) & .ui$iniDf$err == "prop"]
      .fullTheta[add_name] <- fit@results@respar[1]
      .fullTheta[prop_name] <- fit@results@respar[2]
    }
  }
  .ret$fullTheta <- .fullTheta

  # 2. etaObf
  eta_names_fit <- colnames(fit@results@map.eta)
  eta_names_ui <- sapply(eta_names_fit, function(name) {
    theta_name <- sub("^eta\\.", "", name)
    .ui$muRefTable$eta[.ui$muRefTable$theta == theta_name]
  })
  etaObf <- as.data.frame(fit@results@map.eta)
  colnames(etaObf) <- eta_names_ui
  etaObf$ID <- unique(.ret$dataSav$ID)
  etaObf <- etaObf[, c("ID", eta_names_ui), drop = FALSE]
  etaObf$OBJI <- NA_real_
  .ret$etaObf <- etaObf

  # 3. omega
  # Extract omega matrix corresponding to UI eta names
  .eta_names_ui_all <- .ui$iniDf$name[!is.na(.ui$iniDf$neta1)]
  .omega <- matrix(0, nrow = length(.eta_names_ui_all), ncol = length(.eta_names_ui_all))
  dimnames(.omega) <- list(.eta_names_ui_all, .eta_names_ui_all)
  for (i in seq_along(.eta_names_ui_all)) {
    eta_i <- .eta_names_ui_all[i]
    theta_i <- .ui$muRefTable$theta[.ui$muRefTable$eta == eta_i]
    if (length(theta_i) == 1) {
      for (j in seq_along(.eta_names_ui_all)) {
        eta_j <- .eta_names_ui_all[j]
        theta_j <- .ui$muRefTable$theta[.ui$muRefTable$eta == eta_j]
        if (length(theta_j) == 1) {
          .omega[eta_i, eta_j] <- fit@results@omega[theta_i, theta_j]
        }
      }
    }
  }
  .ret$omega <- .omega

  # 4. covariance matrix
  est_theta_names <- .theta_names[!.structural_thetas$fix]
  cov_matrix <- fit@results@MCOV
  if (!is.null(cov_matrix) && nrow(cov_matrix) == length(est_theta_names)) {
    dimnames(cov_matrix) <- list(est_theta_names, est_theta_names)
    .ret$cov <- cov_matrix
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
  .ret <- nlmixr2est::nlmixr2CreateOutputFromUi(.ret$ui, data = .ret$origData, control = .ret$control, table = .ret$table, env = .ret, est = "saemix")

  .env <- .ret$env
  .env$method <- "saemix"
  # Embed saemix object into the final fit environment
  assign("saemix", fit, envir = .env)

  .ret
}
