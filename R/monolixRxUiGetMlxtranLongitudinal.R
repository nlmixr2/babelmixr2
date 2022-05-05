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