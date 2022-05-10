#' Based on the name and ui, figure out what type of data is being provided
#'
#' @param name compartment name
#' @param ui User interface function for retrieving covariate/regressor info
#' @return Monolix use type (string)
#' @author Matthew L. Fidler
#' @noRd
.monolixMapDataUse <- function(name, ui) {
  .use <- switch(name,
                 "ID"="id",
                 "TIME"="time",
                 "EVID"="eventidentifier",
                 "AMT"="amount",
                 "II"="ii",
                 "DV"="observation",
                 "CENS"="censored",
                 "LIMIT"="limit",
                 "YTYPE"="observationtype",
                 "ADM"="administration",
                 "SS"="steadystate",
                 NA_character_)
  if (is.na(.use)) {
    # Determine if this is a regressor or a mu-referenced covariate
    .cov <- ui$saemMuRefCovariateDataFrame
    if (name %in% .cov$covariate) return("{use=covariate, type=continuous}")
    if (name %in% ui$allCovs) return("{use=regressor}")
    return("{use=ignore}")
  }
  if (.use == "steadystate") {
    paste0(.use,", nbdoses=", rxode2::rxGetControl(ui, "nbSSDoses", 5))
  } else if (.use == "observation") {
    .predDf <- ui$predDf
  }
  .use <- paste0("use=", .use)
  .use
}





#' @export
rxUiGet.mlxtranDatafile <- function(x, ...) {
  .ui <- x[[1]]
  # This assumes the data are generated from nlmixr2extra
  .hasRate <- rxode2::rxGetControl(ui, ".hasRate", NA)
  if (is.na(.hasRate)) {
    .rateData <- NULL
  } else if (.hasRate) {
    .rateData <- "RATE"
  } else {
    .rateData <- "TINF"
  }
  .hasCens <- rxode2::rxGetControl(ui, ".hasCens", FALSE)
  .censData <- NULL
  if (.hasCens) {
    .censData <- "CENS"
  }

  .hasLimit <- rxode2::rxGetControl(ui, ".hasLimit", FALSE)
  .limitData <- NULL
  if (.hasLimit) {
    .limitData <- "LIMIT"
  }
  .col0 <- c("ID", "TIME", "EVID", "AMT", "II", "DV", "ADM", "YTYPE", "SS", .rateData,
             .censData, .limitData,
             .ui$allCovs, "nlmixrRowNums")

  .use <- vapply(.col0, .monolixMapDataUse,
                 character(1),
                 ui=.ui, USE.NAMES=FALSE)
  .ret <- c("<DATAFILE>", "",
            "[FILEINFO]",
            paste0("file='", rxUiGet.monolixDataFile(x, ...), "'"),
            "delimter = tab",
            paste0("header = {", paste(.col0, collapse=", "), "}"),"",
            "[CONTENT]")
}
