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
  .ret
}

.mlxtranIndividualDef <- function(var, est, muRefCurEval) {
  .w <- which(muRefCurEval$parameter == est)
  if (length(.w) != 1) stop("duplicate/missing parameter in `muRefCurEval`", call.=FALSE)
  .curEval <- muRefCurEval$curEval[.w]
  .low <- muRefCurEval$low[.w]
  .hi <- muRefCurEval$hi[.w]
  if (is.na(.low) && !is.na(.hi)) .low <- 0
}

#' @export
rxUiGet.mlxtranModelIndividual <- function(x, ...) {
  .ui <- x[[1]]
  .split <- rxUiGet.getSplitMuModel(x, ...)
  .muRef <- c(.split$pureMuRef, .split$taintMuRef)
  ## .muRef <- lapply(seq_along(.muRef), function(.i){
  ##   .est <- names(.muRef)[.i]
  ##   .var <- setNames(.muRef[.i], NULL)
  ##   .createMuRefPkBlock(.var, .est, .ui$muRefCurEval)
  ## })

  paste0("[INDIVIDUAL]\n",
         "input={}\n\n",
         "DEFINITION:\n",
         )
}

#' @export
rxUiGet.monolixModelFileName <- function(x, ...) {
  .ui <- x[[1]]
  paste0(.ui$modelName, ".txt")
}



