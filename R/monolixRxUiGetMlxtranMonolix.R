#' @export
rxUiGet.mlxtranMonolix <- function(x, ...) {
  .ui <- x[[1]]
  .r <- c(
    "<MONOLIX>",
    "",
    "[TASKS]",
    "populationParameters()",
    "individualParameters(method = {conditionalMean, conditionalMode})",
    "fim(method = Linearization)",
    "logLikelihood(method = Linearization)",
    "plotResult(method = {outputplot, indfits, obspred, residualsscatter, residualsdistribution, parameterdistribution, covariatemodeldiagnosis, randomeffects, covariancemodeldiagnosis, saemresults})",
    "",
    "[SETTINGS]",
    "GLOBAL:",
    paste0("exportpath = '", rxUiGet.monolixExportPath(x, ...), "'"),
    "POPULATION:",
    paste0("exploratoryautostop = ", ifelse(rxode2::rxGetControl(.ui, "exploratoryautostop", FALSE), "yes", "no")),
    paste0("smoothingautostop = ", ifelse(rxode2::rxGetControl(.ui, "smoothingautostop",FALSE), "yes", "no")),
    paste0("burniniterations = ", rxode2::rxGetControl(.ui, "burniniterations", 5)),
    paste0("exploratoryiterations = ", rxode2::rxGetControl(.ui, "exploratoryiterations", 250)),
    paste0("simulatedannealingiterations = ", rxode2::rxGetControl(.ui, "simulatedannealingiterations", 250)),
    paste0("smoothingiterations = ", rxode2::rxGetControl(.ui, "smoothingiterations",200)),
    paste0("exploratoryalpha = ", rxode2::rxGetControl(.ui, "exploratoryalpha",0.0)),
    paste0("exploratoryinterval = ", rxode2::rxGetControl(.ui, "exploratoryinterval", 200)),
    paste0("omegatau = ", rxode2::rxGetControl(.ui, "omegatau", 0.95)),
    paste0("errormodeltau = ", rxode2::rxGetControl(.ui, "errormodeltau", 0.95)))
  .var <- rxode2::rxGetControl(.ui, "variability", "none")
  if (.var != "none") {
    .r <- c(.r, paste0("variability = ", .var))
  }
  paste(.r, collapse="\n")
}
