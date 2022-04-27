#' @export
rxUiGet.getMonolixResidualDistribution <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  #   .varp <- str2lang(paste0("rx_pred_", .pred1[["var"]]))
  paste0("rx_prd_", .predDf$var, "={distribution = normal, prediction = rx_pred_", .predDf$var, ", errorModel=")
}
#' @export
rxUiGet.getMonolixModel0 <- function(x, ...) {
  .ui <- x[[1]]
  .split <- rxUiGet.getSplitMuModel(x, ...)
  .predDf <- .ui$predDf
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
  .def <- paste(paste0(names(.resMod), "_pred= {distribution = normal, prediction = ", names(.resMod), ", errorModel=",
                       monolixGetErr(.resMod, uif, control), "}"), collapse="\n")
  .mod <- rxToMonolix(.norm)
}
