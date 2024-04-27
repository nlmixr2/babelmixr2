#' @export
rxUiGet.monolixModelName <- function(x, ...) {
  .ui <- x[[1]]
  .modelName <- rxode2::rxGetControl(.ui, "modelName", NULL)
  if (is.null(.modelName)) {
    .modelName <- .ui$modelName
  } else {
    assign("modelName", .modelName, .ui)
  }
  if (isTRUE(checkmate::checkCharacter(.modelName, len=1, any.missing=FALSE))) {
    return(.modelName)
  }
  "x"
}
#' @export
rxUiGet.monolixExportPath <- function(x, ...) {
  .ui <- x[[1]]
  # Handle monolix2rx as well
  .mlxtran <- monolix2rx::.monolixGetMlxtran(.ui)
  if (inherits(.mlxtran, "monolix2rxMlxtran")) {
    .wd <- attr(.mlxtran, "dirn")
    if (!checkmate::testDirectoryExists(.wd)) .wd <- getwd()
    withr::with_dir(.wd, {
      .exportPath <- .mlxtran$MONOLIX$SETTINGS$GLOBAL$exportpath
      return(path.expand(file.path(.wd, .exportPath)))
    })
  }
  .extra <- ""
  .num <- rxode2::rxGetControl(.ui, ".modelNumber", 0)
  if (.num > 0) {
    .extra <- sprintf("-%03d", .num)
  }
  .base <- paste0(rxUiGet.monolixModelName(x, ...), .extra, "-monolix")
  if (rxode2::rxGetControl(.ui, "absolutePath", FALSE)) {
    file.path(getwd(), .base)
  } else {
    file.path(.base)
  }
}

#' @export
rxUiGet.monolixModelHashFileName <- function(x, ...) {
  paste0(rxUiGet.monolixExportPath(x, ...), ".md5")
}

#' @export
rxUiGet.monolixModelFileName <- function(x, ...) {
  paste0(rxUiGet.monolixExportPath(x, ...), ".txt")
}
#' @export
rxUiGet.monolixDataFile <- function(x, ...) {
  paste0(rxUiGet.monolixExportPath(x, ...), ".csv")
}

#' @export
rxUiGet.monolixQs <- function(x, ...) {
  file.path(rxUiGet.monolixExportPath(x, ...), "nlmixr.qs")
}

#' @export
rxUiGet.monolixCvParam <- function(x, ...) {
  file.path(rxUiGet.monolixExportPath(x, ...), "ChartsData", "Saem", "CvParam.txt")
}

#' @export
rxUiGet.monolixRunLock <- function(x, ...) {
  paste0(rxUiGet.monolixExportPath(x, ...), ".run-lock")
}

#' @export
rxUiGet.monolixMlxtranFile <- function(x, ...) {
  paste0(rxUiGet.monolixExportPath(x, ...), ".mlxtran")
}

#' @export
rxUiGet.mlxtranModel <- function(x, ...) {
  .ui <- x[[1]]
  # note there is some categorical covariates that are not taken care of here...
  paste0("<MODEL>\n\n",
         rxUiGet.mlxtranModelCovariate(x, ...),
         rxUiGet.mlxtranModelIndividual(x, ...),"\n\n",
         rxUiGet.mlxtranModelLongitudinal(x, ...))
}

#' @export
rxUiGet.mlxtran <- function(x, ...) {
  paste(c(rxUiGet.mlxtranDatafile(x, ...),"",
          rxUiGet.mlxtranModel(x, ...),"",
          rxUiGet.mlxtranFit(x, ...),"",
          rxUiGet.mlxtranParameter(x, ...),"",
          rxUiGet.mlxtranMonolix(x, ...)),
        collapse="\n")
}
