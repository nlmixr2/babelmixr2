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
    .ret <- .ret[, !(names(.ret) %in% "YTYPE")]
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
  # - $adjObf Should the objective function value be adjusted
  env$adjObf <- rxode2::rxGetControl(.ui, "adjObf", TRUE)
  # - $objective objective function value
  env$objective <- NA_real_
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
  # Last to try to ensure all files have been exported
  # - $cov For covariance
  .cov <- .ui$monolixCovariance
  if (!is.null(.cov)) {
    env$cov <- .cov
    # - $covMethod for the method of calculating the covariance
    env$covMethod <- rxode2::rxGetControl(.ui, ".covMethod", "Monolix")
  }

  # When running the focei problem to create the nlmixr object, you also need a
  #  foceiControl object
  .monolixControlToFoceiControl(env)
  env <- nlmixr2est::nlmixr2CreateOutputFromUi(env$ui, data=env$origData, control=env$control, table=env$table, env=env, est="monolix")
  .env <- env$env
  .env$method <- "monolix"

  .env$adj <- .env$nobs*log(2 * pi)
  .objf2 <- .ui$monolixObjf
  .objf <- .objf2 - .env$adj
  .llik <- -(.objf2) / 2
  attr(.llik, "df") <- attr(.env$logLik, "df")
  attr(.llik, "nobs") <- .env$nobs
  class(.llik) <- "logLik"
  .env$logLik <- .llik
  .tmp <- data.frame(
          OBJF = .objf, AIC = .objf2 + 2 * attr(get("logLik", .env), "df"),
          BIC = .objf2 + log(.env$nobs) * attr(get("logLik", .env), "df"),
          "Log-likelihood" = as.numeric(.llik), check.names = FALSE)
  nlmixr2est::nlmixrAddObjectiveFunctionDataFrame(env, .tmp, .env$ofvType)
  env
}

.lixoftStarted <- NA
.hasLixoftConnectors <- function() {
  if (is.na(.lixoftStarted)) {
    if (!requireNamespace("lixoftConnectors", quietly = TRUE)) {
      assignInMyNamespace(".lixoftStarted", FALSE)
      return(invisible(FALSE))
    }

    .x <- try(lixoftConnectors::initializeLixoftConnectors(software = "monolix", force=TRUE), silent=TRUE)
    if (inherits(.x, "try-error")) {
      assignInMyNamespace(".lixoftStarted", FALSE)
    } else {
      assignInMyNamespace(".lixoftStarted", TRUE)
    }
  }
  invisible(.lixoftStarted)
}

#' Run NONMEM using either the user-specified command or function
#'
#' @param ui The nlmixr2 UI object for running
#' @param monolix are we actually running monolix
#' @return NULL
#' @noRd
.monolixRunner <- function(ui) {
  cmd <- rxode2::rxGetControl(ui, "runCommand", "")
  if (is.character(cmd)) {
    cmd <- .monolixRunCommand
  } else if (is.na(cmd)) {
    .minfo("run Monolix manually and rerun nlmixr()")
    return(NULL)
  } else if (!is.function(cmd)) {
    stop("invalid value for monolixControl(runCommand=)",
         call.=FALSE)
  }
  cmd(mlxtran=ui$monolixMlxtranFile, directory=ui$monolixExportPath, ui=ui)
  NULL
}

.monolixFamilyFit <- function(env, ...) {
  .ui <- env$ui
  .control <- .ui$control
  .data <- env$data
  .ret <- new.env(parent=emptyenv())
  .ret$table <- env$table
  .ret$monolixControl <- .control
  .tmp  <- bblDatToMonolix(.ui, .data, table=env$table, rxControl=.control$rxControl, env=.ret)
  .ret$monolixData <- .monolixFormatData(.tmp$monolix, .ui)
  .tmp <- .tmp$adm
  if (length(.tmp$adm) == 0) {
    .tmp <- structure(list(adm = integer(0),
                           cmt = integer(0),
                           type = structure(integer(0), .Label = c("empty", "modelRate", "modelDur", "infusion", "bolus"), class = "factor"),
                           f = double(0),
                           dur=double(0),
                           lag=double(0),
                           rate=double(0)),
                      class = "data.frame", row.names = integer(0))
  } else {
    .tmp$f <- NA_real_
    .tmp$dur <- NA_real_
    .tmp$lag <- NA_real_
    .tmp$rate <- NA_real_
  }
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
  .modelText <- .ui$monolixModel
  .mlxtranText <- .ui$mlxtran
  .dataDf <- .ret$monolixData
  .hashMd5 <- digest::digest(list(.modelText, .mlxtranText, .dataDf))
  .foundModelName <- FALSE
  .hashFile <- .ui$monolixModelHashFileName
  while (!.foundModelName) {
    if (!file.exists(.hashFile)) {
      .foundModelName <- TRUE
    } else {
      if (readLines(.hashFile) == .hashMd5) {
        .foundModelName <- TRUE
      } else {
        .num <- rxode2::rxGetControl(.ui, ".modelNumber", 0) + 1
        rxode2::rxAssignControlValue(.ui, ".modelNumber", .num)
        .hashFile <- .ui$monolixModelHashFileName # regenerate hash file name
      }
    }
  }
  .csv <- .ui$monolixDataFile
  # Update if model name has changed
  .modelText <- .ui$monolixModel
  .mlxtranText <- .ui$mlxtran


  .qs <- .ui$monolixQs
  .exportPath <- .ui$monolixExportPath
  .model <- .ui$monolixModelFileName
  .mlxtran <- .ui$monolixMlxtranFile
  .runLock <- .ui$monolixRunLock

  .cmd <- rxode2::rxGetControl(.ui, "runCommand", "")
if (checkmate::testFileExists(.qs)) {
    .minfo("load saved nlmixr2 object")
    .ret <- qs::qread(.qs)
    if (!exists("parHistData", .ret$env)) {
      .tmp <- .ret$ui$monolixParHistory
      if (is.null(.tmp)) {
        .minfo("monolix parameter history needs exported charts, please export charts")
      } else {
        .tmp$type <- "Unscaled"
        assign("parHistData", .tmp, .ret$env)
        .minfo("monolix parameter history integrated into fit object")
        qs::qsave(.ret, .qs)
      }
    }
    return(.ret)
  } else if (!checkmate::testFileExists(.model)) {
    .minfo("writing monolix files")
    writeLines(text=.modelText, con=.model)
    writeLines(text=.mlxtranText, con=.mlxtran)
    writeLines(text=.hashMd5, con=.hashFile)
    write.csv(.dataDf, file=.csv, na = ".", row.names = FALSE)
    .minfo("done")
    if (!rxode2::rxGetControl(.ui, "run", TRUE)) {
      .minfo("only exported Monolix mlxtran, txt model and data")
      return(invisible())
    }
    .runLS <- FALSE
    if (!identical(.cmd, "")) {
      .monolixRunner(ui=.ui)
      if (is.na(.cmd)) {
        return(invisible())
      }
    } else {
      if (.hasLixoftConnectors()) {
        .x <- try(lixoftConnectors::loadProject(.mlxtran), silent=TRUE)
        if (inherits(.x, "try-error")) {
          stop("lixoftConnectors cannot load mlxtran",
               call.=FALSE)
        }
        .minfo("lixoftConnectors::runScenario()")
        lixoftConnectors::runScenario()
        .minfo("done")
        .runLS <- TRUE
      } else if (dir.exists(.exportPath)) { # needs to skip for tests
      } else if (!interactive()) {
        # Don't wait when running in a script or test
        print(.exportPath)
        stop("setup monolix's run command")
      } else {
        .minfo("run monolix manually or stop and setup monolix's run command")
      }
    }
  } else {
    if (is.na(.cmd)) {
      .minfo(paste0("leaving alone monolix files because '", .model, "' is present"))
      return(invisible())
    }
    .minfo(paste0("assuming monolix is running because '", .model, "' is present"))
  }
  if (!dir.exists(.exportPath)) {
    .minfo(paste0("waiting for monolix output (", .exportPath, ")"))
    .i <- 0
    while (!dir.exists(.exportPath)) {
      .i <- .i + 1
      message(".", appendLF=FALSE)
      if (.i %% 50 == 0) {
        message(paste0(.i, "\n"), appendLF=TRUE)
      } else if (.i %% 10 == 0) {
        message("|", appendLF=FALSE)
      }
      Sys.sleep(1)
    }
    message("")
  }
  .ret <- .monolixFinalizeEnv(.ret, .ui)
  if (inherits(.ret, "nlmixr2FitData")) {
    .msg <- .monolixMergePredsAndCalcRelativeErr(.ret)
    .msg$message <- c(.msg$message,
                      paste0("monolix model: '", .mlxtran, "'"))
    .tmp <- .ret$ui$monolixParHistory
    assign("message", paste(.msg$message, collapse="\n    "), envir=.ret$env)
    if (is.null(.tmp)) {
      .minfo("monolix parameter history needs exported charts, please export charts")
    } else {
      .tmp$type <- "Unscaled"
      assign("parHistData", .tmp, .ret$env)
      .minfo("monolix parameter history integrated into fit object")
      qs::qsave(.ret, .qs)
    }
    qs::qsave(.ret, .qs)
  }
  return(.ret)
}

.monolixRunCommand <- function(mlxtran, directory, ui) {
  cmd <- rxode2::rxGetControl(ui, "runCommand", "")
  if (cmd != "") {
    fullCmd <- paste(cmd, mlxtran)
    .minfo(paste0("run Monolix: ", fullCmd))
    system(fullCmd)
  } else {
    stop("run Monolix manually and rerun nlmixr() or setup Monolix's run command")
  }
}

#' @export
nlmixr2Est.monolix <- function(env, ...) {
  .model <- nlmixr2est::.uiApplyMu2(env)
  .ui <- env$ui
  rxode2::assertRxUiMuRefOnly(.ui, " for the estimation routine 'monolix'", .var.name=.ui$modelName)
  .ui <- rxode2::rxUiDecompress(env$ui)
  nlmixr2est::nmObjUiSetCompressed(FALSE)
  on.exit({nlmixr2est::nmObjUiSetCompressed(TRUE)})
  assign("ui", .ui, envir=env)
  on.exit({
    assign("ui", rxode2::rxUiCompress(env$ui), envir=env)
  })
  rxode2::assertRxUiTransformNormal(.ui, " for the estimation routine 'monolix'", .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'monolix'", .var.name=.ui$modelName)
  rxode2::assertRxUiEstimatedResiduals(.ui, " for the estimation routine 'monolix'", .var.name=.ui$modelName)
  .monolixFamilyControl(env, ...)
    nlmixr2est::nmObjUiSetCompressed(FALSE)

  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  }, add=TRUE)
  nlmixr2est::.uiFinalizeMu2(.monolixFamilyFit(env, ...), .model)
}
attr(nlmixr2Est.monolix, "covPresent") <- TRUE
