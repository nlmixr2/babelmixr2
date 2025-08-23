#' Monolix get compartment information
#'
#' @param ui rxode2 user interface
#'
#' @return the compartment specification for defined compartments in
#'   adm dataset
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.monolixGetCompartmentInformation <- function(ui) {
  .adm <- .monolixGetAdm(ui)
  .cmts <- sort(unique(.adm$cmt))
  .state <- rxode2::rxState(ui)[.cmts]
  paste(paste0("compartment(cmt=", .cmts, ", amount=", .state, ")"), collapse="\n")
}
#' Monolix get PK macros for I
#'
#' @param i integer for adm to process
#'
#' @param adm adm data frame
#'
#' @param state A character vector of the states
#'
#' @return macro for adm id i
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.monolixGetPkMacrosForI <- function(i, adm, state) {
  .adm <- adm[i, ]
  .type <- paste(.adm$type)
  if (.type %in% c("bolus", "infusion")) {
    paste0("depot(type=", .adm$adm, ", target=", state[.adm$cmt],
           ", Tlag=", ifelse(is.na(.adm$lag), "0", .adm$lag), ", p=",
           ifelse(is.na(.adm$f), "1", .adm$f), ")")
  } else if (.type == "modelRate") {
    paste0("depot(type=", .adm$adm, ", target=", state[.adm$cmt],
           ", Tk0=amtDose/", .adm$rate,
           ", Tlag=", ifelse(is.na(.adm$lag), "0", .adm$lag), ", p=",
           ifelse(is.na(.adm$f), "1", .adm$f), ")")
  } else if (.type == "modelDur") {
    paste0("depot(type=", .adm$adm, ", target=", state[.adm$cmt],
           ", Tk0=", .adm$dur,
           ", Tlag=", ifelse(is.na(.adm$lag), "0", .adm$lag), ", p=",
           ifelse(is.na(.adm$f), "1", .adm$f), ")")
  } else if (.type == "empty") {
    paste0("empty(adm=", .adm$adm, "target=", state[.adm$cmt], ")")
  }
}
#' Get the PK macros
#'
#' @param ui rxode2 ui
#' @return Monolix macros for dosing
#' @author Matthew L. Fidler
#' @noRd
.monolixGetPkMacros <- function(ui) {
  .adm <- .monolixGetAdm(ui)
  .state <- rxode2::rxState(ui)
  paste(vapply(seq_along(.adm$adm),
               .monolixGetPkMacrosForI,
               character(1), adm=.adm, state=.state,
               USE.NAMES=FALSE), collapse="\n")
}

#' @export
rxUiGet.monolixModel <- function(x, ...) {
  .ui <- x[[1]]
  .split <- .ui$getSplitMuModel
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
  .regress <- .ui$allCovs
  .cov <- .ui$saemMuRefCovariateDataFrame
  .regress <- .regress[!(.regress %in% .cov)]
  .regressors <- ""
  if (length(.ui$allCovs) > 0) {
    .regressors <- paste0("\n", paste(paste0(.regress, "= {use=regressor}"), collapse="\n"))
  }
  paste0("DESCRIPTION:\n",
         paste0("model translated from `babelmixr2` and `nlmixr2` function ", .ui$modelName, " to ", .txtFile, "\n\n"),
         "[LONGITUDINAL]\n",
         "input={", paste(setNames(c(.split$pureMuRef, .split$taintMuRef, .regress), NULL), collapse=",") , "}",
         .regressors,
          ifelse(rxode2::rxGetControl(.ui, "stiff", FALSE), "\n\nodeType = stiff", ""),
         "\n\nPK:\n; Define compartments with administrations\n",
         .monolixGetCompartmentInformation(.ui),
         "\n; Define PK macros\n",
         .monolixGetPkMacros(.ui),
         "\n\nEQUATION:\n", .mod,
         "\n\nOUTPUT:\noutput={",
         paste(.monolixResponses, collapse=", "), "}\n")
}
attr(rxUiGet.monolixModel, "rstudio") <- ""
