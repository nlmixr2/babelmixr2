#' PharmML namespace URIs, by supported version
#'
#' @noRd
.pmlNsTable <- list(
  "0.9" = c(
    pharmml = "http://www.pharmml.org/pharmml/0.9/PharmML",
    ct      = "http://www.pharmml.org/pharmml/0.9/CommonTypes",
    math    = "http://www.pharmml.org/pharmml/0.9/Maths",
    mdef    = "http://www.pharmml.org/pharmml/0.9/ModelDefinition",
    mstep   = "http://www.pharmml.org/pharmml/0.9/ModellingSteps",
    ds      = "http://www.pharmml.org/pharmml/0.9/Dataset",
    design  = "http://www.pharmml.org/pharmml/0.9/TrialDesign",
    po      = "http://www.pharmml.org/probonto/ProbOnto"
  )
)

#' Supported PharmML versions
#'
#' `babelmixr2` vendors the PharmML schemas, so this is the set of versions
#' that can be written and validated without network access.
#'
#' @return character vector of supported PharmML versions
#'
#' @author Matthew L. Fidler
#'
#' @export
#' @examples
#' pharmmlVersions()
pharmmlVersions <- function() {
  names(.pmlNsTable)
}

#' Namespace URIs for a PharmML version
#'
#' @param version PharmML version string
#' @return named character vector of namespace URIs
#' @noRd
.pmlNs <- function(version = "0.9") {
  .ret <- .pmlNsTable[[version]]
  if (is.null(.ret)) {
    stop("unsupported PharmML version '", version, "'; supported: ",
         paste(pharmmlVersions(), collapse = ", "),
         call. = FALSE)
  }
  .ret
}

#' Path to the vendored schema root for a PharmML version
#'
#' @param version PharmML version string
#' @return path to the directory holding pharmml.xsd and its siblings
#' @noRd
.pmlSchemaPath <- function(version = "0.9") {
  .pmlNs(version) # validates the version
  .ret <- system.file("pharmml", version, package = "babelmixr2")
  if (!nzchar(.ret)) {
    stop("cannot find vendored PharmML ", version, " schemas", call. = FALSE)
  }
  .ret
}

#' Coerce validation input to an xml_document
#'
#' @param x path, document text, or xml_document
#' @return xml_document
#' @noRd
.pmlAsXmlDocument <- function(x) {
  if (inherits(x, "xml_document")) return(x)
  checkmate::assertCharacter(x, min.len = 1, any.missing = FALSE)
  if (length(x) == 1L && !grepl("<", x, fixed = TRUE)) {
    checkmate::assertFileExists(x)
    return(xml2::read_xml(x))
  }
  xml2::read_xml(paste(x, collapse = "\n"))
}

#' Validate a PharmML document against its schema
#'
#' Validation is entirely offline: the schemas are vendored in the package and
#' their `schemaLocation` references have been rewritten to relative paths.
#' The upstream schema host (`pharmml.org`) no longer serves them.
#'
#' @param x Either a path to a PharmML file, or a length-one character vector
#'   holding the document itself, or an `xml_document`.
#'
#' @param version PharmML version to validate against.  See
#'   `pharmmlVersions()`.
#'
#' @param error When `TRUE` (default) an invalid document raises an error
#'   listing the schema violations.  When `FALSE` the function returns `FALSE`
#'   with the violations in the `"errors"` attribute.
#'
#' @return `TRUE` when the document validates.  When `error = FALSE` and the
#'   document does not validate, `FALSE` with an `"errors"` attribute.
#'
#' @author Matthew L. Fidler
#'
#' @export
#' @examples
#' \dontrun{
#' pharmmlValidate("model.xml")
#' }
pharmmlValidate <- function(x, version = "0.9", error = TRUE) {
  checkmate::assertLogical(error, len = 1, any.missing = FALSE)
  .schemaDir <- .pmlSchemaPath(version)
  .doc <- .pmlAsXmlDocument(x)
  .schema <- xml2::read_xml(file.path(.schemaDir, "pharmml.xsd"))
  .ret <- xml2::xml_validate(.doc, .schema)
  if (!isTRUE(as.logical(.ret)) && error) {
    stop("PharmML document does not validate against the ", version,
         " schema:\n  ",
         paste(attr(.ret, "errors"), collapse = "\n  "),
         call. = FALSE)
  }
  .ret
}
