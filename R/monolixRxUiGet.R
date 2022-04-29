.getMonolixResidual <- function(i, ui, input=FALSE) {
  .pred1 <- ui$predDf[i, ]
  .iniDf <- ui$iniDf
  .cond <- .pred1$cond
  .errType <- as.integer(.pred1$errType)
  .iniDf <- .iniDf[which(.iniDf$condition == .cond), ]
  if (.errType == 1L) { # add
    .tmp <- as.character(.iniDf$name)
    .tmp <- eval(str2lang(paste0("rxToMonolix(", .tmp, ")")))
    if (input) return(.tmp)
    return(paste0("constant(", .tmp, ")"))
  } else if (.errType == 2L) { # prop
    .tmp <- as.character(.iniDf$name)
    .getMonolixResidualAddPar(.tmp)
    .tmp <- eval(str2lang(paste0("rxToMonolix(", .tmp, ")")))
    if (input) return(.tmp)
    return(paste0("proportional(", .tmp, ")"))
  } else if (.errType == 3L) { # pow
    stop("residual type 'pow' is not supported in the nlmixr2 monolix conversion yet",
         call.=FALSE)
  } else if (.errType == 4L) { # add + prop
    .addProp <- .pred1$addProp
    if (.addProp == "default") {
      .addProp <- rxode2::rxGetControl(ui, "addProp", "combined2")
    }
    .add <- .iniDf[.iniDf$err == "add", "name"]
    .add <- eval(str2lang(paste0("rxToMonolix(", .add, ")")))
    .prop <- .iniDf[.iniDf$err == "prop", "name"]
    .prop <- eval(str2lang(paste0("rxToMonolix(", .prop, ")")))
    if (input) return(c(.add, .prop))
    return(paste0(.addProp, "(", .add, ",", .prop, ")"))
  } else if (.errType == 5L) { #  add + pow
    stop("residual type 'add+pow' is not supported in the nlmixr2 monolix conversion yet",
         call.=FALSE)
  } else {
    stop("residual type 'none' is not supported in monolix",
         call.=FALSE)
  }
}

#' @export
rxUiGet.mlxtranModelLongitudinal <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  .err <- vapply(seq_along(.predDf$var),
                 .getMonolixResidual, character(1), ui=.ui, USE.NAMES=FALSE)
  .inputErr <- unlist(lapply(seq_along(.predDf$var),
                             .getMonolixResidual,  ui=.ui, input=TRUE))
  .iniDf <- .ui$iniDf
  paste0("[LONGITUDINAL]\n",
       "input={", paste(.inputErr, collapse=", "), "}\n",
       "file=", rxUiGet.monolixModelFileName(x, ...), "\n\n",
       "DEFINITION:\n",
       paste(paste0("rx_prd_", .predDf$var, "={distribution = normal, prediction = rx_pred_", .predDf$var, ", errorModel=", .err, "}"), collapse='\n'))
}

#' @export
rxUiGet.monolixModelFileName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(.ui$modelName, ".txt")
}

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
  .mod <- rxToMonolix(.norm)
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
         paste(paste0("depot(type=1, target=", .mv$state, ", Tlag=", monolixTlag(.mv$state), ", p=", monolixP(.mv$state), ")"), collapse="\n"),
         "\n\nEQUATION:\n", .mod,
         "\n\nOUTPUT:\noutput={",
         paste(.monolixResponses, collapse=", "), "}\n")
}




