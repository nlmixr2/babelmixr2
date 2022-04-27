#' This is a S3 method for getting the distribution lines for a base rxode2 saem problem
#'
#' @param line Parsed rxode2 model environment
#' @return Lines for the estimation of monolix
#' @author Matthew Fidler
#' @keywords internal
#' @export
nmGetDistributionMonolixLines <- function(line) {
  UseMethod("nmGetDistributionMonolixLines")
}

#' Creates a monolix line object from a predDf line
#'
#' @param x rxode2 ui object
#' @param line Line number for monolix error line object
#' @param len Number of prediction statements
#' @return nmGetDistributionMonolixLines object
#' @author Matthew L. Fidler
#' @noRd
.createMonolixLineObject <- function(x, line) {
  .predDf <- get("predDf", x)
  if (line > nrow(.predDf)) {
    return(NULL)
  }
  .predLine <- .predDf[line, ]
  .ret <- list(x, .predLine)
  class(.ret) <- c(paste(.predLine$distribution), "nmGetDistributionMonolixLines")
  .ret
}

#' @rdname nmGetDistributionMonolixLines
#' @export
nmGetDistributionMonolixLines.rxUi <- function(line) {
  .predDf <- get("predDf", line)
  lapply(seq_along(.predDf$cond), function(c){
    .mod <- .createMonolixLineObject(line, c)
    nmGetDistributionMonolixLines(.mod)
  })
}

#' @rdname nmGetDistributionMonolixLines
#' @export
nmGetDistributionMonolixLines.norm <- function(line) {
  .rx <- line[[1]]
  .pred1 <- line[[2]]
  if (.pred1[["linCmt"]]) {
    stop("linCmt() translation not supported (yet)",
         call.=FALSE)
  }
  .var <- str2lang(.pred1[["var"]])
  .varp <- str2lang(paste0("rx_pred_", .pred1[["var"]]))
  return(list(bquote(.(.varp) <- .(.var))))
}

#' @export
nmGetDistributionMonolixLines.t <- function(line) {
  stop("t isn't supported yet")
}

#' @export
nmGetDistributionMonolixLines.default  <- function(line) {
  stop("Distribution not supported")
}
