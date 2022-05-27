#' This is a S3 method for getting the distribution lines for a base rxode2 saem problem
#'
#' @param line Parsed rxode2 model environment
#' @return Lines for the estimation of nonmem
#' @author Matthew Fidler
#' @keywords internal
#' @export
nmGetDistributionNonmemLines <- function(line) {
  UseMethod("nmGetDistributionNonmemLines")
}

#' Creates a nonmem line object from a predDf line
#'
#' @param x rxode2 ui object
#' @param line Line number for nonmem error line object
#' @param len Number of prediction statements
#' @return nmGetDistributionNonmemLines object
#' @author Matthew L. Fidler
#' @noRd
.createNonmemLineObject <- function(x, line, extra="") {
  .predDf <- get("predDf", x)
  if (line > nrow(.predDf)) {
    return(NULL)
  }
  .predLine <- .predDf[line, ]
  .ret <- list(x, .predLine)
  class(.ret) <- c(paste0(.predLine$distribution, extra), "nmGetDistributionNonmemLines")
  .ret
}

#' @rdname nmGetDistributionNonmemLines
#' @export
nmGetDistributionNonmemLines.rxUi <- function(line) {
  .predDf <- get("predDf", line)
  lapply(seq_along(.predDf$cond), function(c){
    .mod <- .createNonmemLineObject(line, c)
    nmGetDistributionNonmemLines(.mod)
  })
}

.nonmemResponses <- NULL

#' @rdname nmGetDistributionNonmemLines
#' @export
nmGetDistributionNonmemLines.norm <- function(line) {
  .env <- line[[1]]
  .pred1 <- line[[2]]
  if (.pred1[["linCmt"]]) {
    stop("linCmt() translation not supported (yet)",
         call.=FALSE)
  }
  # Take out transformation; not supported right for multiple endpoint models
  .ret <- vector("list", 1)
  .ret[[1]] <- bquote(rx_pred_f_ ~ .(rxode2::.rxGetPredictionF(.env, .pred1)))
  .ret
}


#' @export
nmGetDistributionNonmemLines.t <- function(line) {
  stop("t isn't supported yet")
}

#' @export
nmGetDistributionNonmemLines.default  <- function(line) {
  stop("Distribution not supported")
}
