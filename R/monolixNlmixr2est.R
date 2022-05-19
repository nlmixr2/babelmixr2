#' Get the monolix control statement and install it into the ui
#'
#' @param env Environment with ui in it
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.monolixFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- monolixControl()
  }
  if (!inherits(.control, "monolixControl")){
    .control <- do.call(babelmixr2::monolixControl, .control)
  }
  assign("control", .control, envir=.ui)
}

.monolixFormatData <- function(data, ui) {
  .ret <- data
  .ret$SS <- ifelse(.ret$SS == 0, NA_real_, .ret$SS)
  .ret$YTYPE <- ifelse(.ret$YTYPE == 0, NA_real_, .ret$YTYPE)
  .ret$ADM <- ifelse(.ret$ADM == 0, NA_real_, .ret$ADM)
  if (all(is.na(.ret$SS))) {
    .ret <- .ret[, !(names(.ret) %in% c("SS", "II"))]
  }
  if (length(ui$predDf$cond) == 1L) {
    .ret <- .ret[, !(names(.ret) %in% c("YTYPE"))]
  }
  .n <- names(.ret)
  rxode2::rxAssignControlValue(ui, ".hasRate",
                               ifelse(any(.n == "RATE"), TRUE, ifelse(any(.n == "TINF"), FALSE, NA)))
  rxode2::rxAssignControlValue(ui, ".hasCens", any(.n == "CENS"))
  rxode2::rxAssignControlValue(ui, ".hasLimit", any(.n == "LIMIT"))
  rxode2::rxAssignControlValue(ui, ".hasIi", any(.n == "II"))
  rxode2::rxAssignControlValue(ui, ".hasSs", any(.n == "SS"))
  .ret
}

.monolixFinalizeEnv <- function(env, oldUi) {
  # The environment needs:
  .iniDf <- oldUi$monolixIniDf
  .ui <- new.env(parent=emptyenv())
  for (n in ls(envir=oldUi, all.names=TRUE)) {
    assign(n, get(n, envir=oldUi), envir=.ui)
  }
  assign("iniDf", .iniDf, envir=.ui)
  class(.ui) <- class(oldUi)
  # - $table for table options -- already present
  # - $origData -- Original Data -- already present
  # - $dataSav -- Processed data from .foceiPreProcessData --already present
  # - $idLvl -- Level information for ID factor added -- already present
  env$ui <- .ui
  # - $ui for ui fullTheta Full theta information
  env$fullTheta <- .ui$monolixFullTheta
  # - $etaObf data frame with ID, etas and OBJI
  env$etaObf <- .ui$monolixEtaObf
  # - $cov For covariance
  .cov <- .ui$monolixCovariance
  if (!is.null(.cov)) {
    env$cov <- .cov
    # - $covMethod for the method of calculating the covariance
    env$covMethod <- rxode2::rxGetControl(.ui, ".covMethod", "Monolix")
  }
  # - $adjObf Should the objective function value be adjusted
  env$adjObf <- rxode2::rxGetControl(.ui, "adjObf", TRUE)
  # - $objective objective function value
  env$objective <- .ui$monolixObjf
  # - $extra Extra print information
  env$extra <- paste0(" ver ", env$ui$monolixOutputVersion)
  # - $method Estimation method (for printing)
  env$method <- "Monolix"
  # - $omega Omega matrix
  env$omega <- .ui$monolixOmega
  # - $theta Is a theta data frame
  env$theta <- .ui$monolixTheta
  # - $model a list of model information for table generation.  Needs a `predOnly` model
  env$model <- .ui$ebe
  # - $message Message for display
  env$message <- ""
  # - $est estimation method
  env$est <- "monolix"
  # - $ofvType (optional) tells the type of ofv is currently being used
  #env$ofvType
  env$ofvType <- .ui$monolixObjfType
  # When running the focei problem to create the nlmixr object, you also need a
  #  foceiControl object
  .monolixControlToFoceiControl(env)
  env <- nlmixr2est::nlmixr2CreateOutputFromUi(env$ui, data=env$origData, control=env$control, table=env$table, env=env, est="monolix")
  .env <- env$env
  .env$method <- "monolix"
  env
}

.lixoftStarted <- NA
.lixoftNs <- NULL

.hasLixoftConnectors <- function() {
  if (is.na(.lixoftStarted)) {
    if (!requireNamespace("lixoftConnectors", quietly = TRUE)) {
      assignInMyNamespace(".lixoftStarted", FALSE)
      return(invisible(FALSE))
    }
    .l <- loadNamespace("lixoftConnectors")
    .x <- try(.l$initializeLixoftConnectors(software = "monolix"), silent=TRUE)
    if (inherits(.x, "try-error")) {
      assignInMyNamespace(".lixoftStarted", FALSE)
    } else {
      assignInMyNamespace(".lixoftStarted", TRUE)
      assignInMyNamespace(".lixoftNs", .l)
    }
  }
  invisible(.lixoftStarted)
}

.monolixFamilyFit <- function(env, ...) {
  .ui <- env$ui
  .control <- .ui$control
  .data <- env$data
  .ret <- new.env(parent=emptyenv())
  .ret$table <- env$table
  .ret$monolixControl <- .control
  .tmp  <- bblDatToMonolix(.ui, .data, table=env$table, env=.ret)
  .ret$monolixData <- .monolixFormatData(.tmp$monolix, .ui)
  .tmp <- .tmp$adm
  .tmp$f <- NA_real_
  .tmp$dur <- NA_real_
  .tmp$lag <- NA_real_
  .tmp$rate <- NA_real_
  rxode2::rxAssignControlValue(.ui, ".adm", .tmp)

  # Now make sure time varying covariates are not considered
  # mu-referenced items
  .et <- rxode2::etTrans(.ret$dataSav, .ui$mv0, addCmt=TRUE)
  .nTv <- attr(class(.et), ".rxode2.lst")$nTv
  if (is.null(.nTv)) {
    .tv <- names(.et)[-seq(1, 6)]
    .nTv <- length(.tv)
  } else {
    .tv <- character(0)
    if (.nTv != 0) {
      .tv <- names(.et)[-seq(1, 6)]
    }
  }
  .muRefCovariateDataFrame <- .ui$muRefCovariateDataFrame
  if (length(.tv) > 0) {
    # Drop time-varying covariates
    .muRefCovariateDataFrame <- .muRefCovariateDataFrame[!(.muRefCovariateDataFrame$covariate %in% .tv), ]
  }
  assign("muRefFinal", .muRefCovariateDataFrame, .ui)
  assign("timeVaryingCovariates", .tv, .ui)
  on.exit({
    if (exists("muRefFinal", envir=.ui)) {
      rm(list="muRefFinal", envir=.ui)
    }
    if (exists("timeVaryingCovariates", envir=.ui)) {
      rm(list="timeVaryingCovariates", envir=.ui)
    }
  })
  .qs <- .ui$monolixQs
  .exportPath <- .ui$monolixExportPath
  .csv <- .ui$monolixDataFile
  .model <- .ui$monolixModelFileName
  .mlxtran <- .ui$monolixMlxtranFile
  .runLock <- .ui$monolixRunLock

  if (file.exists(.qs)) {
    .minfo("load saved nlmixr2 object")
    return(qs::qread(.qs))
  } else if (!file.exists(.model)) {
    .minfo("writing monolix files")
    writeLines(text=.ui$monolixModel, con=.model)
    writeLines(text=.ui$mlxtran, con=.mlxtran)
    write.csv(.ret$monolixData, file=.csv, na = ".", row.names = FALSE)
    .minfo("done")
    .runLS <- FALSE
    .cmd <- rxode2::rxGetControl(.ui, "runCommand", "")
    if (.cmd != "") {
      .minfo(paste0("run monolix: ", sprintf(.cmd, .mlxtran)))
      system(sprintf(.cmd, .mlxtran))
    } else {
      if (.hasLixoftConnectors()) {
        .l <- .lixoftNs
        .x <- try(.l$loadProject(.mlxtran), silent=TRUE)
        if (inherits(.x, "try-error")) {
          stop("lixoftConnectors cannot load mlxtran",
               call.=FALSE)
        }
        .l$runScenario()
        .runLS <- TRUE
      } else {
        .minfo("run monolix manually or stop and setup monolix's run command")
      }
    }
  }
  if (!dir.exists(.exportPath)) {
    .minfo("waiting for monolix output")
    .i <- 0
    while (!dir.exists(.exportPath)) {
      .i <- .i + 1
      message(".", appendLF=FALSE)
      if (.i %% 50 == 0) {
        message(paste0(.i, "\n"), appendLF=TRUE)
      } else if (.i %% 10 == 0) {
        message("|", appendLF=FALSE)
        if (.runLS) {
          .status <- .lixoftNs$getLastRunStatus()
          if (.status$report != ""){
            message("")
            message(.status$report)
            break;
          }
          if (.status$status) {
            break;
          }
        }
      }
      Sys.sleep(1)
    }
    message("")
  }
  .ret <- .monolixFinalizeEnv(.ret, .ui)
  if (inherits(.ret, "nlmixr2FitData")) {
    qs::qsave(.ret, .qs)
  }
  return(.ret)
}

nlmixr2Est.monolix <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiTransformNormal(.ui, " for the estimation routine 'monolix'", .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'monolix'", .var.name=.ui$modelName)
  rxode2::assertRxUiEstimatedResiduals(.ui, " for the estimation routine 'monlix'", .var.name=.ui$modelName)
  .monolixFamilyControl(env, ...)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  }, add=TRUE)
  .monolixFamilyFit(env, ...)
}
