#' Control options for PharmML export
#'
#' @param version PharmML version to write.  See `pharmmlVersions()`.
#'
#' @param validate When `TRUE` (default) the emitted document is validated
#'   against the PharmML schema before it is returned, and an invalid document
#'   is an error.  Validation is offline and takes a fraction of a second, so
#'   there is rarely a reason to turn it off.
#'
#' @param dataFile Name the PharmML document should use when referring to the
#'   exported dataset.  When `NULL` this is derived from the model name.
#'
#' @param description Free text written into the document's `ct:Description`.
#'   When `NULL` a short provenance line naming `babelmixr2` is written.
#'
#' @param writeData When `TRUE` (default) `as.pharmml(file=)` also writes the
#'   NONMEM-format dataset next to the model file.  A PharmML document refers
#'   to its data by path rather than embedding it, so a document written
#'   without its dataset is incomplete.
#'
#' @return a list of PharmML export options
#'
#' @author Matthew L. Fidler
#'
#' @export
#' @examples
#' pharmmlControl(validate = FALSE)
pharmmlControl <- function(version = "0.9",
                           validate = TRUE,
                           dataFile = NULL,
                           description = NULL,
                           writeData = TRUE) {
  checkmate::assertCharacter(version, len = 1, any.missing = FALSE)
  checkmate::assertLogical(validate, len = 1, any.missing = FALSE)
  checkmate::assertLogical(writeData, len = 1, any.missing = FALSE)
  if (!is.null(dataFile)) {
    checkmate::assertCharacter(dataFile, len = 1, any.missing = FALSE)
  }
  if (!is.null(description)) {
    checkmate::assertCharacter(description, len = 1, any.missing = FALSE)
  }
  .pmlNs(version) # errors on an unsupported version
  .ret <- list(version = version,
               validate = validate,
               dataFile = dataFile,
               description = description,
               writeData = writeData)
  class(.ret) <- "pharmmlControl"
  .ret
}
