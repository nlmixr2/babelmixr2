#' Block ids used across the PharmML ModelDefinition
#'
#' These are fixed rather than generated: a PharmML document written by
#' `babelmixr2` always uses the same block names, which keeps `blkIdRef`
#' resolution trivial and the output diffable.
#'
#' @noRd
.pmlBlk <- c(
  variabilityParameter = "vm1",
  variabilityResidual  = "vm2",
  covariate            = "cm1",
  parameter            = "pm1",
  structural           = "sm1"
)

#' Does this model have between-subject variability?
#'
#' @param ui rxode2 UI
#' @return TRUE when the model has at least one eta
#' @noRd
.pharmmlHasEta <- function(ui) {
  .iniDf <- ui$iniDf
  any(!is.na(.iniDf$neta1))
}

#' PharmML VariabilityModel blocks
#'
#' `babelmixr2` only writes models whose random effects are on `id` (the same
#' restriction the Monolix and NONMEM writers apply), so there is at most one
#' parameter-variability level.  The residual level is always written, because
#' every supported model has an estimated residual.
#'
#' @param ui rxode2 UI
#'
#' @param indent Indent depth
#'
#' @return character(1) holding the `VariabilityModel` blocks
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.pharmmlVariabilityModel <- function(ui, indent = 0L) {
  .ret <- character(0)
  if (.pharmmlHasEta(ui)) {
    .ret <- c(.ret,
              .pmlNode("mdef:VariabilityModel",
                       attrs = c(blkId = .pmlBlk[["variabilityParameter"]],
                                 type = "parameterVariability"),
                       children = .pmlNode("mdef:Level",
                                           attrs = c(referenceLevel = "true",
                                                     symbId = "id"),
                                           indent = indent + 1L),
                       indent = indent))
  }
  .ret <- c(.ret,
            .pmlNode("mdef:VariabilityModel",
                     attrs = c(blkId = .pmlBlk[["variabilityResidual"]],
                               type = "residualError"),
                     children = .pmlNode("mdef:Level",
                                         attrs = c(symbId = "residual"),
                                         indent = indent + 1L),
                     indent = indent))
  paste(.ret, collapse = "\n")
}

#' @export
rxUiGet.pharmmlVariabilityModel <- function(x, ...) {
  .pharmmlVariabilityModel(x[[1]])
}
attr(rxUiGet.pharmmlVariabilityModel, "rstudio") <- "character"
