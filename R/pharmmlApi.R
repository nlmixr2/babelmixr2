#' The namespace declarations for a PharmML root element
#'
#' @param version PharmML version
#' @return character(1) of `xmlns:` attributes, one per line
#' @noRd
.pharmmlNsAttrs <- function(version = "0.9") {
  .ns <- .pmlNs(version)
  .default <- paste0('xmlns="', .ns[["pharmml"]], '"')
  .rest <- vapply(setdiff(names(.ns), "pharmml"), function(.p) {
    paste0('    xmlns:', .p, '="', .ns[[.p]], '"')
  }, character(1), USE.NAMES = FALSE)
  paste(c(.default, .rest), collapse = "\n")
}

#' Assemble a complete PharmML document
#'
#' @param ui rxode2 UI
#' @param data The dataset, or `NULL` for a model-only document
#' @param control `pharmmlControl()` options
#' @return character(1) holding the whole document
#' @noRd
.pharmmlDocument <- function(ui, data = NULL, control = pharmmlControl()) {
  .version <- control$version
  .name <- ui$modelName
  if (is.null(.name) || !nzchar(.name)) .name <- "model"
  .desc <- control$description
  if (is.null(.desc)) {
    .desc <- paste0("Model '", .name,
                    "' translated to PharmML ", .version,
                    " by babelmixr2 ",
                    utils::packageVersion("babelmixr2"))
  }
  .dataFile <- control$dataFile
  if (is.null(.dataFile)) .dataFile <- paste0(.name, ".csv")

  .body <- c(.pmlText("ct:Name", .name, 1L),
             .pmlText("ct:Description", .desc, 1L),
             .pmlNode("IndependentVariable", attrs = c(symbId = "t"), indent = 1L),
             .pharmmlModelDefinition(ui, data = data, indent = 1L))
  if (!is.null(data)) {
    .body <- c(.body,
               .pharmmlTrialDesign(ui, data, dataFile = .dataFile, indent = 1L),
               .pharmmlModellingSteps(ui, indent = 1L))
  }

  .ret <- paste0('<?xml version="1.0" encoding="UTF-8"?>\n',
                 "<PharmML ", .pharmmlNsAttrs(.version), "\n",
                 '    writtenVersion="', .version, '" id="', .pharmmlId(.name), '">\n',
                 paste(.body, collapse = "\n"), "\n",
                 "</PharmML>\n")
  if (isTRUE(control$validate)) {
    pharmmlValidate(.ret, version = .version)
  }
  .ret
}

#' A document id derived from the model name
#'
#' `ct:IdType` is an NCName, so a name starting with a digit or containing a
#' space would be rejected.
#'
#' @param name model name
#' @return an NCName-safe id
#' @noRd
.pharmmlId <- function(name) {
  .ret <- gsub("[^A-Za-z0-9._-]", "_", name)
  if (!grepl("^[A-Za-z_]", .ret)) .ret <- paste0("i", .ret)
  .ret
}

#' Convert an rxode2 expression to PharmML
#'
#' This is the expression-level translator, the PharmML counterpart of
#' [rxToNonmem()] and [rxToMonolix()].  To translate a whole model, use
#' [as.pharmml()].
#'
#' @param x Expression, or a character vector holding one
#'
#' @param ui rxode2 ui.  Optional: when supplied, symbols are qualified with
#'   the PharmML block that owns them.
#'
#' @return PharmML syntax, as a character string
#'
#' @author Matthew L. Fidler
#'
#' @export
#' @examples
#' rxToPharmml("ka * depot")
rxToPharmml <- function(x, ui = NULL) {
  if (!is.null(ui)) {
    ui <- rxode2::assertRxUi(ui)
    ui <- rxode2::rxUiDecompress(ui)
  }
  if (is(substitute(x), "character")) {
    force(x)
  } else if (is(substitute(x), "{")) {
    x <- deparse1(substitute(x))
  } else if (!is.character(x)) {
    x <- deparse1(substitute(x))
  }
  .lst <- lapply(as.list(str2lang(paste0("{", paste(x, collapse = "\n"), "}")))[-1],
                 function(.e) .rxToPharmml(.e, ui))
  paste(unlist(.lst), collapse = "\n")
}

#' Translate an nlmixr2 model to PharmML
#'
#' PharmML is a standard XML description of a pharmacometric model, so this
#' makes an `nlmixr2` model archivable and exchangeable without the receiving
#' tool needing to understand `rxode2` syntax.
#'
#' The emitted document is validated against the PharmML schema before it is
#' returned.  Constructs PharmML cannot express -- non-normal residuals,
#' inter-occasion variability, mixture models -- raise an error naming the
#' construct rather than producing a silently wrong document.
#'
#' A PharmML document refers to its dataset by path rather than embedding it,
#' so when `file` is given the NONMEM-format dataset is written alongside the
#' model unless `pharmmlControl(writeData = FALSE)`.
#'
#' @param model An `rxode2`/`nlmixr2` model function or ui
#'
#' @param data The dataset the model is to be estimated with.  When `NULL`
#'   only the `ModelDefinition` is written: a valid PharmML model description
#'   with no trial design or estimation step.
#'
#' @param file Path to write the document to.  When `NULL` (default) nothing
#'   is written and the document is returned.
#'
#' @param control `pharmmlControl()` options
#'
#' @return the PharmML document, invisibly when `file` is given
#'
#' @author Matthew L. Fidler
#'
#' @export
#' @examples
#' \dontrun{
#' one.cmt <- function() {
#'   ini({
#'     tka <- log(1.57); tcl <- log(2.72); tv <- log(31.5)
#'     eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl + eta.cl)
#'     vc <- exp(tv + eta.v)
#'     d/dt(depot) <- -ka * depot
#'     d/dt(center) <- ka * depot - cl / vc * center
#'     cp <- center / vc
#'     cp ~ add(add.sd)
#'   })
#' }
#' as.pharmml(one.cmt, nlmixr2data::theo_sd, file = "theo.xml")
#' }
as.pharmml <- function(model, data = NULL, file = NULL,
                       control = pharmmlControl()) {
  if (!inherits(control, "pharmmlControl")) {
    stop("'control' must come from pharmmlControl()", call. = FALSE)
  }
  .ui <- rxode2::assertRxUi(model, " to translate to PharmML")
  .ui <- rxode2::rxUiDecompress(.ui)
  if (!is.null(data)) {
    checkmate::assertDataFrame(data, min.rows = 1)
  }
  .ret <- .pharmmlDocument(.ui, data, control)
  class(.ret) <- "babelmixr2pharmml"
  if (is.null(file)) return(.ret)

  writeLines(.ret, file)
  if (!is.null(data) && isTRUE(control$writeData)) {
    .dataFile <- control$dataFile
    if (is.null(.dataFile)) {
      .name <- .ui$modelName
      if (is.null(.name) || !nzchar(.name)) .name <- "model"
      .dataFile <- paste0(.name, ".csv")
    }
    .pharmmlWriteData(.ui, data, file.path(dirname(file), basename(.dataFile)))
  }
  invisible(.ret)
}

#' @export
print.babelmixr2pharmml <- function(x, ...) {
  cat(paste(x, collapse = "\n"))
  invisible(x)
}
