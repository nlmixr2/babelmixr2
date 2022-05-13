#' Get the saem control statement and install it into the ui
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
    .control <- saemControl()
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
  .n <- names(.ret)
  rxode2::rxAssignControlValue(ui, ".hasRate",
                               ifelse(any(.n == "RATE"), TRUE, ifelse(any(.n == "TINF"), FALSE, NA)))
  rxode2::rxAssignControlValue(ui, ".hasCens", any(.n == "CENS"))
  rxode2::rxAssignControlValue(ui, ".hasLimit", any(.n == "LIMIT"))
  .ret
}

.monolixFamilyFit <- function(env, ...) {
  .ui <- env$ui
  .control <- .ui$control
  .data <- env$data
  .ret <- new.env(parent=emptyenv())
  .ret$table <- env$table
  .tmp  <- nlmixr2extra::nlmixrDataToMonolix(.ui, .data, table=env$table, env=.ret)
  .ret$monolixData <- .monolixFormatData(.tmp$monolix)
  .tmp <- .tmp$adm
  rxode2::rxAssignControlValue(.ui, ".adm", .tmp)
  .tmp$f <- NA_real_
  .tmp$dur <- NA_real_
  .tmp$lag <- NA_real_
  .tmp$rate <- NA_real_

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
    .cmd <- rxode2::rxGetControl(.ui, "runCommand", "")
    if (.cmd != "") {
      writeLines("", .runLock)
      system(sprintf(.cmd, .mlx))
    } else {
      message("run monolix manually or setup monolix's run command")
    }
    return(invisible())
  } else if (!file.exists(.exportPath)) {
    if (file.exists(.runLock)) {
      .minfo(paste0("may still be runnning '", .runLock, "'"))
    } else {
      .minfo(paste0("the export location '", .exportPath, "' doens't have files in it yet"))
    }
    return(invisible())
  } else {

  }
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
