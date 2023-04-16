#' Convert an object to a nlmixr2 fit object
#'  
#' @param x Object to convert
#' @param ... Other arguments
#' @return nlmixr2 fit object
#' @export 
#' @author Matthew L. Fidler
#' @examples 
as.nlmixr2 <- function(x, ...) {
  UseMethod("as.nlmixr2")
}
#' @export
as.nlmixr <- as.nlmixr2

as.nlmixr2.default <- function(x, ...) {
  stop("cannot figure out how to create an nlmixr2 object from the input",
       call.=FALSE)
}
