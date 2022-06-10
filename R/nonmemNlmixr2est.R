#' Get the nonmem control statement and install it into the ui
#'
#' @param env Environment with ui in it
#' @param ... Other arguments
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nonmemFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- nonmemControl()
  }
  if (!inherits(.control, "nonmemControl")){
    .control <- do.call(babelmixr2::nonmemControl, .control)
  }
  assign("control", .control, envir=.ui)
}

.nonmemFormatData <- function(data, ui) {
  .ret <- data
  .ret$SS <- ifelse(.ret$SS == 0, NA_real_, .ret$SS)
  if (all(is.na(.ret$SS))) {
    .ret <- .ret[, !(names(.ret) %in% c("SS", "II"))]
  }
  .n <- names(.ret)
  rxode2::rxAssignControlValue(ui, ".hasRate", any(.n == "RATE"))
  rxode2::rxAssignControlValue(ui, ".hasCens", any(.n == "CENS"))
  rxode2::rxAssignControlValue(ui, ".hasLimit", any(.n == "LIMIT"))
  rxode2::rxAssignControlValue(ui, ".hasIi", any(.n == "II"))
  rxode2::rxAssignControlValue(ui, ".hasSs", any(.n == "SS"))
  .ret
}

.nonmemFinalizeEnv <- function(env, oldUi) {
  # The environment needs:
  .iniDf <- oldUi$nonmemIniDf
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
  env$fullTheta <- .ui$nonmemFullTheta
  # - $etaObf data frame with ID, etas and OBJI
  env$etaObf <- .ui$nonmemEtaObf
  # - $cov For covariance
  .cov <- .ui$nonmemCovariance
  if (!is.null(.cov)) {
    env$cov <- .cov
    # - $covMethod for the method of calculating the covariance
    env$covMethod <- rxode2::rxGetControl(.ui, ".covMethod", "Nonmem")
  }
  # - $adjObf Should the objective function value be adjusted
  env$adjObf <- rxode2::rxGetControl(.ui, "adjObf", TRUE)
  # - $objective objective function value
  env$objective <- .ui$nonmemObjf
  # - $extra Extra print information
  env$extra <- paste0(" ver ", env$ui$nonmemOutputVersion)
  # - $method Estimation method (for printing)
  env$method <- "Nonmem"
  # - $omega Omega matrix
  env$omega <- .ui$nonmemOmega
  # - $theta Is a theta data frame
  env$theta <- .ui$nonmemTheta
  # - $model a list of model information for table generation.  Needs a `predOnly` model
  env$model <- .ui$ebe
  # - $message Message for display
  env$message <- ""
  # - $est estimation method
  env$est <- "nonmem"
  # - $ofvType (optional) tells the type of ofv is currently being used
  #env$ofvType
  env$ofvType <- .ui$nonmemObjfType
  # When running the focei problem to create the nlmixr object, you also need a
  #  foceiControl object
  .nonmemControlToFoceiControl(env)
  env <- nlmixr2est::nlmixr2CreateOutputFromUi(env$ui, data=env$origData, control=env$control, table=env$table, env=env, est="nonmem")
  .env <- env$env
  .env$method <- "nonmem"
  env
}

.nonmemFamilyFit <- function(env, ...) {
  .ui <- env$ui
  .control <- .ui$control
  .data <- env$data
  .ret <- new.env(parent=emptyenv())
  .ret$table <- env$table
  .ret$nonmemControl <- .control
  .tmp  <- bblDatToNonmem(.ui, .data, table=env$table, env=.ret)
  .ret$nonmemData <- .nonmemFormatData(.tmp, .ui)
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
  .nmctl <- .ui$nonmemModel
  .contra <- .ui$nonmemContra
  .hashMd5 <- digest::digest(list(.nmctl, .contra, .ret$nonmemData))
  .foundModelName <- FALSE

  .hashFile <- file.path(.ui$nonmemExportPath, .ui$nonmemHashFile)
  while (!.foundModelName) {
    if (!file.exists(.hashFile)) {
      .foundModelName <- TRUE
    } else {
      if (readLines(.hashFile) == .hashMd5) {
        .foundModelName <- TRUE
      } else {
        .num <- rxode2::rxGetControl(.ui, ".modelNumber", 0) + 1
        rxode2::rxAssignControlValue(.ui, ".modelNumber", .num)
        .hashFile <- file.path(.ui$nonmemExportPath, .ui$nonmemHashFile)
      }
    }
  }
  .exportPath <- .ui$nonmemExportPath
  if (!dir.exists(.exportPath)) dir.create(.exportPath)
  .csv <- file.path(.exportPath, .ui$nonmemCsv)
  .nmctlFile <- file.path(.exportPath, .ui$nonmemNmctl)
  .contraFile <- file.path(.exportPath, .ui$nonmemContraName)
  .ccontraFile <- file.path(.exportPath, .ui$nonmemCcontraName)
  .qs <- file.path(.exportPath, .ui$nonmemQs)

  if (file.exists(.qs)) {
    .minfo("load saved nlmixr2 object")
    .ret <- qs::qread(.qs)
    return(.ret)
  } else if (!file.exists(.nmctlFile)) {
    .minfo("writing nonmem files")
    writeLines(text=.nmctl, con=.nmctlFile)
    writeLines(text=.hashMd5, con=.hashFile)
    if (!is.null(.contra)) {
      writeLines(text=.contra, con=.contraFile)
      .ccontra <- .ui$nonmemCcontra
      writeLines(text=.ccontra, con=.ccontraFile)
    }
    write.csv(.ret$nonmemData, file=.csv, na = ".", row.names = FALSE,
              quote=FALSE)
    .minfo("done")
  }
  .cmd <- rxode2::rxGetControl(.ui, "runCommand", "")
  if (.cmd != "") {
    .arg <- paste0(.ui$nonmemNmctl, " ", .ui$nonmemNmlst)
    .minfo(paste0("run NONMEM: ", sprintf(.cmd, .arg)))
    system(sprintf(.cmd, .arg))
  } else {
    .minfo("run NONMEM manually or setup NONMEM's run command")
  }

  stop("need to run and test as next step")
  ## if (!dir.exists(.exportPath)) {
  ##   .minfo("waiting for nonmem output")
  ##   .i <- 0
  ##   while (!dir.exists(.exportPath)) {
  ##     .i <- .i + 1
  ##     message(".", appendLF=FALSE)
  ##     if (.i %% 50 == 0) {
  ##       message(paste0(.i, "\n"), appendLF=TRUE)
  ##     } else if (.i %% 10 == 0) {
  ##       message("|", appendLF=FALSE)
  ##     }
  ##     Sys.sleep(1)
  ##   }
  ##   message("")
  ## }
  ## .ret <- .nonmemFinalizeEnv(.ret, .ui)
  ## if (inherits(.ret, "nlmixr2FitData")) {
  ##   .tmp <- .ret$ui$nonmemParHistory
  ##   if (is.null(.tmp)) {
  ##     .minfo("nonmem parameter history needs expoted charts, please export charts")
  ##   } else {
  ##     assign("parHist", .tmp, .ret$env)
  ##     .minfo("nonmem parameter history integrated into fit object")
  ##     qs::qsave(.ret, .qs)
  ##   }
  ##   qs::qsave(.ret, .qs)
  ## }
  ## return(.ret)
}

#' @export
nlmixr2Est.nonmem <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiTransformNormal(.ui, " for the estimation routine 'nonmem'", .var.name=.ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'nonmem'", .var.name=.ui$modelName)
  rxode2::assertRxUiEstimatedResiduals(.ui, " for the estimation routine 'nonmem'", .var.name=.ui$modelName)
  .nonmemFamilyControl(env, ...)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  }, add=TRUE)
  .nonmemFamilyFit(env, ...)
}
