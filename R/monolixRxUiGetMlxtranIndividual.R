#' mlxtran switch distributions
#'
#' @param curEval Current evaluation
#' @return Monolix distribution type
#' @author Matthew L. Fidler
#' @noRd
.mlxTranCurEvalToDistribution <- function(curEval) {
  .ret <- switch(ifelse(curEval == "", "add", curEval),
                 exp="logNormal",
                 expit="logitNormal",
                 probitInv="probitNormal",
                 add="normal",
                 NA_character_)
  if (is.na(.ret))
    stop(paste0("monolix translation of '", curEval, "' is unknown"),
         call.=FALSE)
  paste0("distribution=", .ret)
}
#' Determine if the estimate is a Population estimate only
#'
#' @param est Estimated value
#' @param muRefTable Mu ref table
#' @return `TRUE` for population only estimates
#' @author Matthew L. Fidler
#' @noRd
.mlxTranIsPopOnly <- function(est, muRefTable) {
  !(est %in% muRefTable$theta)
}


.mlxTranInputForIndividual <- NULL

#' Get the omega variability component name
#'
#' @param var monolix variable name
#' @param est typical value estimation variable name
#' @param muRefTable Mu reference table
#' @return monolix sd=expression or no-variability
#' @author Matthew L. Fidler
#' @noRd
.mlxTranGetVaraibility <- function(var, est, muRefTable) {
  if (.mlxTranIsPopOnly(est, muRefTable)) {
    "no-variability"
  } else {
    assignInMyNamespace(".mlxTranInputForIndividual",
                        c(.mlxTranInputForIndividual, paste0("omega_", var)))
    paste0("sd=omega_", var)
  }
}

.mlxTranGetLimits <- function(curEval, .low, .hi) {
  if (curEval == "expit") {
    c(paste0("min=", ifelse(is.na(.low), "0", .low)),
      paste0("max=", ifelse(is.na(.hi), "1", .hi)))
  } else {
    NULL
  }
}
#' Get mlxtran mu reference covariate definitions
#'
#' @param var Monolix modeled variable
#' @param est Monolix estimated variable
#' @param muRefCovariateDataFrame Mu referenced covariate data frame
#' @return covariate and coefficient expressions for mu-referenced
#'   covariates. NULL otherwise
#' @author Matthew L. Fidler
#' @noRd
.mlxTranGetMuRefCovariate <- function(var, est, muRefCovariateDataFrame) {
  if (is.null(muRefCovariateDataFrame)) return(NULL)
  .covs <- muRefCovariateDataFrame[muRefCovariateDataFrame == est, ]
  if (length(.covs$covariate) == 0) return(NULL)
  .cov <- paste0("covariate = {", paste(.covs$covariate, collapse=", "), "}")
  .coef <- paste0("beta_", var, "_", .covs$covariate)
  assignInMyNamespace(".mlxTranInputForIndividual",
                      c(.mlxTranInputForIndividual, .coef))
  .coef <- paste0("coefficient = {", paste(.coef, collapse=", "), "}")
  c(.cov, .coef)
}
#' Get the mlxtran individual estimate for the mu-referenced variables
#'
#' @param var Monolix modeled variable
#' @param est Theta estimated variable
#' @param muRefCurEval This is the mu reference current evaluation
#'   function
#' @param muRefTable This is the mu-reference table
#' @param muRefCovariateDataFrame mu-ref covariate data frame
#' @return A single line that gives the individual definition
#' @author Matthew L. Fidler
#' @noRd
.mlxtranIndividualDef <- function(var, est, muRefCurEval, muRefTable, muRefCovariateDataFrame) {
  .w <- which(muRefCurEval$parameter == est)
  if (length(.w) != 1) stop("duplicate/missing parameter in `muRefCurEval`", call.=FALSE)
  .curEval <- muRefCurEval$curEval[.w]
  .low <- muRefCurEval$low[.w]
  .hi <- muRefCurEval$hi[.w]
  if (is.na(.low) && !is.na(.hi)) .low <- 0
  assignInMyNamespace(".mlxTranInputForIndividual",
                      c(.mlxTranInputForIndividual, paste0(var, "_pop")))
  paste0(var, " = {", paste(c(.mlxTranCurEvalToDistribution(.curEval),
          .mlxTranGetLimits(.curEval, .low, .hi),
          paste0("typical=", var, "_pop"),
          .mlxTranGetMuRefCovariate(var, est, muRefCovariateDataFrame),
          .mlxTranGetVaraibility(var, est, muRefTable)),
        collapse=", "), "}")
}
#' Get individual mu ref eta monolix names
#'
#' @param ui roxde2 ui
#' @param num eta number
#' @param muRef mu-reference theta to monolix variables
#' @return monolix variable for eta
#' @author Matthew L. Fidler
#' @noRd
.mlxtranGetIndividualMuRefEtaMonolixName <- function(ui, num, muRef) {
  .iniDf <- ui$iniDf
  .etaName <- .iniDf[which(.iniDf$neta1 == num & .iniDf$neta2 == num), "name"]
  .w <- which(ui$muRefTable$eta == .etaName)
  if (length(.w) == 1) {
    return(setNames(muRef[ui$muRefTable$theta[.w]], NULL))
  }
  stop("only should get eta names where the mu referencing is known",
       call.=FALSE)
}
#' Get individual correlation definition for var1 and var2
#'
#' @param ui rxode2 ui
#' @param num1 eta number 1
#' @param num2 eta number 2
#' @param muRef Mu reference number
#' @return Gives the correlation definition, ie r(p1, p2)=corr_p1_p2
#' @author Matthew L. Fidler
#' @noRd
.mlxtranGetIndividualCorDefinition <- function(ui, num1, num2, muRef) {
  .par1 <- .mlxtranGetIndividualMuRefEtaMonolixName(ui, num1, muRef)
  .par2 <- .mlxtranGetIndividualMuRefEtaMonolixName(ui, num2, muRef)
  .cor <- paste0("corr_", .par1, "_", .par2)
  assignInMyNamespace(".mlxTranInputForIndividual", c(.mlxTranInputForIndividual, .cor))
  paste0("r(", .par1, ", ", .par2, ")=", .cor)
}
#' Get individual correlation statement (if needed)
#'
#' @param ui rxode2 ui
#' @param muRefs Mu reference to monolix varaible
#' @return the correlation= statement in monolix
#' @author Matthew L. Fidler
#' @noRd
.mlxtranIndividualCor <- function(ui, muRefs) {
  .eta <- ui$iniDf[!(is.na(ui$iniDf$neta1)),, drop=FALSE]
  if (length(.eta$neta1) == 0L) stop("need eta for monolix model",
                                     call.=FALSE)
  .w <- which(.eta$neta1 != .eta$neta2)
  if (length(.w) > 0) {
    .eta <- .eta[.w, ]
    return(paste0("correlation = {",
                  paste(c("level=id",
                          vapply(seq_along(.eta$neta1), function(i) {
                            .cur <- .eta[i, ]
                            .mlxtranGetIndividualCorDefinition(ui, .cur$neta1, .cur$neta2, muRef)
                          }, character(1), USE.NAMES=TRUE)), collapse=", "), "}"))
  }
  ""
}

#' @export
rxUiGet.mlxtranModelIndividual <- function(x, ...) {
  .ui <- x[[1]]
  .split <- rxUiGet.getSplitMuModel(x, ...)
  .muRef <- c(.split$pureMuRef, .split$taintMuRef)
  assignInMyNamespace(".mlxTranInputForIndividual", NULL)
  .muRefCov <- .ui$saemMuRefCovariateDataFrame
  .def <- vapply(seq_along(.muRef), function(.i){
    .est <- names(.muRef)[.i]
    .var <- setNames(.muRef[.i], NULL)
    .mlxtranIndividualDef(.var, .est, .ui$muRefCurEval, .ui$muRefTable, .muRefCov)
  }, character(1), USE.NAMES=FALSE)
  .def <- paste(.def, collapse="\n")
  .cor <- .mlxtranIndividualCor(.ui, muRefs)
  if (.cor != "") .def <- paste0(.def, "\n", .cor)
  paste0("[INDIVIDUAL]\n",
         "input={", paste(.mlxTranInputForIndividual, collapse=", "), "}\n\n",
         "DEFINITION:\n",
         .def
         )
}
