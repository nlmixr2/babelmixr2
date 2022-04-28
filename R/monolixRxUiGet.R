.getMonolixResidual <- function(ui, i) {
  .pred1 <- ui$predDf[i, ]
  .iniDf <- ui$iniDf
  .cond <- .pred1$cond
  .errType <- .pred1$errType
  .iniDf <- .iniDf[which(.iniDf$condition == .cond), ]
  if (.errType == 1L) { # add

  } else if (.errType == 2L) { # prop

  } else if (.errType == 3L) { # pow

  } else if (.errType == 4L) { # add + prop
    .addProp <- .predDf$addProp
    if (.addProp == "default") {
      .addProp <- rxode2::rxGetControl(ui, "addProp", "combined2")
    }
    .add <- .iniDf[.iniDf$err == "add", "name"]
    .add <- eval(str2lang(paste0("rxToMonolix(", .add, ")")))
    .prop <- .iniDf[.iniDf$err == "prop", "name"]
    .prop <- eval(parse(text=paste0("rxToMonolix(", .prop, ")")))
    .lst[[length(.lst) + 1]] <- data.frame(name=.ini[.ini$err == "prop", "name"], errName=.prop)
    assignInMyNamespace(".monolixErrs", c(.monolixErrs, .add, .prop))
    assignInMyNamespace(".monolixGetErr", .lst)
    return(paste0(control$addProp, "(", .add, ",", .prop, ")"))

  } else if (.errType == 5L) { #  add + pow

  } else {
    stop("residual type 'none' is not supported in monolix",
         call.=FALSE)
  }
}

#' @export
rxUiGet.getMonolixResidualDistribution <- function(x, ...) {
  .ui <- x[[1]]
  .predDf <- .ui$predDf
  #   .varp <- str2lang(paste0("rx_pred_", .pred1[["var"]]))
  .iniDf <- .ui$iniDf
  paste0("rx_prd_", .predDf$var, "={distribution = normal, prediction = rx_pred_", .predDf$var, ", errorModel=", "}")
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
  .mod <- rxToMonolix(.norm)
  .mod
}
