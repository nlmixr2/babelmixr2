#' PharmML CovariateModel block
#'
#' Every covariate the model references is declared as a continuous covariate.
#' Categorical covariates are not detected yet: `rxode2` carries no type
#' information for them at this level, so they would have to be declared by the
#' caller.  Continuous is the safe default because it is what a numeric column
#' in the dataset means.
#'
#' @param ui rxode2 UI
#'
#' @param indent Indent depth
#'
#' @return character(1) holding the `CovariateModel` block, or `""` when the
#'   model has no covariates
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.pharmmlCovariateModel <- function(ui, data = NULL, indent = 0L) {
  .covs <- ui$allCovs
  if (length(.covs) == 0L) return("")
  .info <- .pharmmlCovariateInfo(ui, data)
  .children <- vapply(.covs, function(.c) {
    .body <- if (identical(.info[[.c]]$type, "categorical")) {
      .pmlNode("mdef:Categorical",
               children = vapply(.info[[.c]]$levels, function(.l) {
                 .pmlNode("mdef:Category", attrs = c(catId = .l))
               }, character(1), USE.NAMES = FALSE))
    } else {
      .pmlNode("mdef:Continuous")
    }
    .pmlNode("mdef:Covariate", attrs = c(symbId = .c), children = .body)
  }, character(1), USE.NAMES = FALSE)
  .pmlNode("mdef:CovariateModel",
           attrs = c(blkId = .pmlBlk[["covariate"]]),
           children = .children,
           indent = indent)
}

#' @export
rxUiGet.pharmmlCovariateModel <- function(x, ...) {
  .pharmmlCovariateModel(x[[1]])
}
attr(rxUiGet.pharmmlCovariateModel, "rstudio") <- "character"

#' Assert that a model can be written as PharmML
#'
#' The same guards the Monolix estimation backend applies, for the same reason:
#' PharmML's structured parameter model assumes mu-referencing, one level of
#' between-subject variability, and estimated residuals.
#'
#' @param ui rxode2 UI
#' @return Nothing, called for the side effect of erroring
#' @noRd
.pharmmlAssertUi <- function(ui) {
  .what <- " for PharmML translation"
  rxode2::assertRxUiMuRefOnly(ui, .what, .var.name = ui$modelName)
  rxode2::assertRxUiTransformNormal(ui, .what, .var.name = ui$modelName)
  rxode2::assertRxUiRandomOnIdOnly(ui, .what, .var.name = ui$modelName)
  rxode2::assertRxUiEstimatedResiduals(ui, .what, .var.name = ui$modelName)
  invisible()
}

#' PharmML ModelDefinition
#'
#' Assembles the five ModelDefinition blocks in the order the schema's
#' `xs:sequence` requires: variability, covariate, parameter, structural,
#' observation.
#'
#' @param ui rxode2 UI
#'
#' @param indent Indent depth
#'
#' @return character(1) holding the `ModelDefinition` element
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.pharmmlModelDefinition <- function(ui, data = NULL, indent = 0L) {
  .pharmmlAssertUi(ui)
  .blocks <- c(.pharmmlVariabilityModel(ui),
               .pharmmlCovariateModel(ui, data),
               .pharmmlParameterModel(ui),
               .pharmmlStructuralModel(ui),
               .pharmmlObservationModel(ui))
  .blocks <- .blocks[nzchar(.blocks)]
  .pmlNode("mdef:ModelDefinition", children = .blocks, indent = indent)
}

#' @export
rxUiGet.pharmmlModelDefinition <- function(x, ...) {
  .pharmmlModelDefinition(x[[1]])
}
attr(rxUiGet.pharmmlModelDefinition, "rstudio") <- "character"
