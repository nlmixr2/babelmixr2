#' @export
rxUiGet.mlxtranMonolix <- function(x, ...) {
  .ui <- x[[1]]
  .r <- c(
    "<MONOLIX>",
    "",
    "[TASKS]",
    "populationParameters()",
    "individualParameters(method = {conditionalMode})",
    paste0("fim(method = ", ifelse(rxode2::rxGetControl(.ui, "useLinearization", TRUE),
                                   "Linearization", "StochasticApproximation"), ")"),
    paste0("logLikelihood(method = ", ifelse(rxode2::rxGetControl(.ui, "useLinearization", TRUE),
                                             "Linearization", "ImportanceSampling"), ")"),
    "plotResult(method = {outputplot, indfits, obspred, residualsscatter, residualsdistribution, parameterdistribution, covariatemodeldiagnosis, randomeffects, covariancemodeldiagnosis, saemresults})",
    "",
    "[SETTINGS]",
    "GLOBAL:",
    paste0("exportpath = '", rxUiGet.monolixExportPath(x, ...), "'"),
    "POPULATION:",
    paste0("exploratoryautostop = ", ifelse(rxode2::rxGetControl(.ui, "exploratoryAutoStop", FALSE), "yes", "no")),
    paste0("smoothingautostop = ", ifelse(rxode2::rxGetControl(.ui, "smoothingAutoStop",FALSE), "yes", "no")),
    paste0("burniniterations = ", rxode2::rxGetControl(.ui, "burnInIterations", 5)),
    paste0("exploratoryiterations = ", rxode2::rxGetControl(.ui, "exploratoryIterations", 250)),
    paste0("simulatedannealingiterations = ", rxode2::rxGetControl(.ui, "simulatedAnnealingIterations", 250)),
    paste0("smoothingiterations = ", rxode2::rxGetControl(.ui, "smoothingIterations",200)),
    paste0("exploratoryalpha = ", rxode2::rxGetControl(.ui, "exploratoryAlpha",0.0)),
    paste0("exploratoryinterval = ", rxode2::rxGetControl(.ui, "exploratoryInterval", 200)),
    paste0("omegatau = ", rxode2::rxGetControl(.ui, "omegaTau", 0.95)),
    paste0("errormodeltau = ", rxode2::rxGetControl(.ui, "errorModelTau", 0.95)))
    #FIXEME: optimizationIterations, tolerance
  .var <- rxode2::rxGetControl(.ui, "variability", "none")
  if (.var != "none") {
    .r <- c(.r, paste0("variability = ", .var))
  }
  paste(.r, collapse="\n")
}
