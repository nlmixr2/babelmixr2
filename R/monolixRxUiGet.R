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
attr(rxUiGet.monolixModelName, "rstudio") <- "character"

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
attr(rxUiGet.monolixExportPath, "rstudio") <- "character"

#' @export
rxUiGet.monolixModelHashFileName <- function(x, ...) {
  paste0(rxUiGet.monolixExportPath(x, ...), ".md5")
}
attr(rxUiGet.monolixModelHashFileName, "rstudio") <- "character"

#' @export
rxUiGet.monolixModelFileName <- function(x, ...) {
  paste0(rxUiGet.monolixExportPath(x, ...), ".txt")
}
attr(rxUiGet.monolixModelFileName,"rstudio") <- "character"

#' @export
rxUiGet.monolixDataFile <- function(x, ...) {
  paste0(rxUiGet.monolixExportPath(x, ...), ".csv")
}
attr(rxUiGet.monolixDataFile, "rstudio") <- "character"

#' @export
rxUiGet.monolixQs <- function(x, ...) {
  file.path(rxUiGet.monolixExportPath(x, ...), "nlmixr.qs")
}
attr(rxUiGet.monolixQs, "rstudio") <- "character"

#' @export
rxUiGet.monolixCvParam <- function(x, ...) {
  file.path(rxUiGet.monolixExportPath(x, ...), "ChartsData", "Saem", "CvParam.txt")
}
attr(rxUiGet.monolixCvParam, "rstudio") <- ""

#' @export
rxUiGet.monolixRunLock <- function(x, ...) {
  paste0(rxUiGet.monolixExportPath(x, ...), ".run-lock")
}
attr(rxUiGet.monolixRunLock, "rstudio") <- ""

#' @export
rxUiGet.monolixMlxtranFile <- function(x, ...) {
  paste0(rxUiGet.monolixExportPath(x, ...), ".mlxtran")
}
attr(rxUiGet.monolixMlxtranFile, "rstudio") <- ""

#' @export
rxUiGet.mlxtranModel <- function(x, ...) {
  .ui <- x[[1]]
  # note there is some categorical covariates that are not taken care of here...
  paste0("<MODEL>\n\n",
         rxUiGet.mlxtranModelCovariate(x, ...),
         rxUiGet.mlxtranModelIndividual(x, ...),"\n\n",
         rxUiGet.mlxtranModelLongitudinal(x, ...))
}
attr(rxUiGet.mlxtranModel, "rstudio") <- ""

#' @export
rxUiGet.mlxtran <- function(x, ...) {
  paste(c(rxUiGet.mlxtranDatafile(x, ...),"",
          rxUiGet.mlxtranModel(x, ...),"",
          rxUiGet.mlxtranFit(x, ...),"",
          rxUiGet.mlxtranParameter(x, ...),"",
          rxUiGet.mlxtranMonolix(x, ...)),
        collapse="\n")
}
attr(rxUiGet.mlxtran, "rstudio") <- ""
