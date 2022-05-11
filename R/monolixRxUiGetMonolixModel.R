#' @export
rxUiGet.monolixModel <- function(x, ...) {
  .ui <- x[[1]]
  .split <- rxUiGet.getSplitMuModel(x, ...)
  assignInMyNamespace(".monolixResponses", NULL)
  # first drop the error lines
  .mainModel <- rxode2::rxCombineErrorLines(.ui,
                    errLines=nmGetDistributionMonolixLines(.ui),
                    paramsLine=NA,
                    modelVars=TRUE,
                    cmtLines=FALSE,
                    dvidLine=FALSE,
                    lstExpr=.split$modelWithDrop,
                    useIf=FALSE)
  .norm <- rxode2::rxNorm(eval(.mainModel))
  .mv <- rxode2::rxModelVars(.ui)
  .mod <- rxToMonolix(.norm, ui=.ui)
  .txtFile <- rxUiGet.monolixModelFileName(x, ...)
  .regressors <- ""
  if (length(.ui$allCov) > 0) {
    .regressors <- paste(paste0(.ui$allCov, "= {use=regressor}"), collapse="\n")
  }
  paste0("DESCRIPTION:\n",
         paste0("model translated from `babelmixr2` and `nlmixr2` function ", .ui$modelName, " to ", .txtFile, "\n\n"),
         "[LONGITUDINAL]\n",
         "input={", paste(setNames(c(.split$pureMuRef, .split$taintMuRef, .ui$allCov), NULL), collapse=",") , "}",
         .regressors,
          ifelse(rxode2::rxGetControl(ui, "stiff", FALSE), "\n\nodeType = stiff", ""),
         "\n\nPK:\n; Define compartment(s)\n",
         paste(paste0("compartment(cmt=", seq_along(.mv$state), ", amount=", .mv$state, ")"), collapse="\n"),
         paste("\n\n;Define depot compartment information\n"),
         paste(paste0("depot(type=", seq_along(.mv$state), ", target=", .mv$state, ", Tlag=", monolixTlag(.mv$state), ", p=", monolixP(.mv$state), ")"), collapse="\n"),
         "\n\nEQUATION:\n", .mod,
         "\n\nOUTPUT:\noutput={",
         paste(.monolixResponses, collapse=", "), "}\n")
}
