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
                 "CMT"="administration",
                 "SS"="steadystate",
                 NA_character_)
  if (is.na(.use)) {
    # Determine if this is a regressor or a mu-referenced covariate
    if (name == "nlmixrRowNums") return("ignore")
    .cov <- ui$saemMuRefCovariateDataFrame
    if (name %in% .cov$covariate) return("regressor")
    return("covariate")
  }
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
  .col0 <- c("ID", "TIME", "EVID", "AMT", "II", "DV", "CMT", "YTYPE", "SS", .rateData,
             .censData, .limitData,
             .ui$allCovs, "nlmixrRowNums")
  .ret <- c("<DATAFILE>", ""
            "[FILEINFO]",
            paste0("file='", rxUiGet.monolixDataFile(x, ...), "'"),
            "delimter = tab",
            paste0("header = {", paste(.col0, collapse=", "), "}"),"",
            "[CONTENT]")
}
