#' Error-model parameter names for one endpoint
#'
#' Error parameters live in `iniDf`, keyed by the endpoint's `condition` with a
#' non-`NA` `err` column saying which role each one plays.
#'
#' @param ui rxode2 UI
#' @param cond endpoint condition (the `predDf$cond` value)
#' @return named character vector, `err` type -> parameter name
#' @noRd
.pharmmlErrParams <- function(ui, cond) {
  .iniDf <- ui$iniDf
  .w <- which(!is.na(.iniDf$err) & .iniDf$condition == cond)
  setNames(.iniDf$name[.w], .iniDf$err[.w])
}

#' Resolve the `addProp` combination form
#'
#' Mirrors the Monolix writer, which resolves `"default"` through
#' `rxGetControl(ui, "addProp", "combined2")`.
#'
#' @param ui rxode2 UI
#' @param predLine one row of `predDf`
#' @return "combined1" or "combined2"
#' @noRd
.pharmmlAddProp <- function(ui, predLine) {
  .ret <- paste(predLine[["addProp"]])
  if (.ret == "default") {
    .ret <- rxode2::rxGetControl(ui, "addProp", "combined2")
  }
  .ret
}

#' The error-model expression for one endpoint, as an R language object
#'
#' PharmML's `Standard` observation model is `u(y) = u(f) + g(...)*eps`, so this
#' returns the `g(...)` term.
#'
#' @param ui rxode2 UI
#' @param predLine one row of `predDf`
#' @return R language object
#' @noRd
.pharmmlErrorModelExpr <- function(ui, predLine) {
  .cond <- paste(predLine[["cond"]])
  .type <- paste(predLine[["errType"]])
  .par <- .pharmmlErrParams(ui, .cond)
  .f <- str2lang(paste(predLine[["var"]]))

  if (.type == "add") {
    return(str2lang(.par[["add"]]))
  }
  if (.type == "prop") {
    return(bquote(.(str2lang(.par[["prop"]])) * .(.f)))
  }
  if (.type == "add + prop") {
    .a <- str2lang(.par[["add"]])
    .b <- str2lang(.par[["prop"]])
    if (.pharmmlAddProp(ui, predLine) == "combined1") {
      return(bquote(.(.a) + .(.b) * .(.f)))
    }
    return(bquote(sqrt(.(.a)^2 + (.(.b) * .(.f))^2)))
  }
  stop("PharmML translation of the residual error type '", .type,
       "' is not supported", call. = FALSE)
}

#' PharmML ObservationModel blocks, one per endpoint
#'
#' @param ui rxode2 UI
#'
#' @param indent Indent depth
#'
#' @return character(1) holding the `ObservationModel` blocks
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.pharmmlObservationModel <- function(ui, indent = 0L) {
  .predDf <- ui$predDf
  .ret <- vapply(seq_len(nrow(.predDf)), function(.i) {
    .line <- .predDf[.i, ]
    .cond <- paste(.line[["cond"]])
    .dist <- paste(.line[["distribution"]])
    if (.dist != "norm") {
      stop("PharmML translation of the '", .dist,
           "' residual distribution is not supported", call. = FALSE)
    }
    if (isTRUE(.line[["linCmt"]])) {
      stop("linCmt() translation is not supported yet", call. = FALSE)
    }
    .eps <- paste0("eps_", .cond)
    .popPars <- vapply(setNames(.pharmmlErrParams(ui, .cond), NULL),
                       function(.p) {
                         .pmlNode("mdef:PopulationParameter", attrs = c(symbId = .p))
                       }, character(1), USE.NAMES = FALSE)

    .pmlNode(
      "mdef:ObservationModel",
      attrs = c(blkId = paste0("om", .i)),
      children = .pmlNode(
        "mdef:ContinuousData",
        children = c(
          .popPars,
          .pmlNode("mdef:RandomVariable", attrs = c(symbId = .eps),
                   children = c(
                     .pmlNode("ct:VariabilityReference",
                              children = .pmlNode(
                                "ct:SymbRef",
                                attrs = c(blkIdRef = .pmlBlk[["variabilityResidual"]],
                                          symbIdRef = "residual"))),
                     .pmlNode("mdef:Distribution",
                              children = .pmlNode("po:ProbOnto",
                                                  attrs = c(name = "StandardNormal1"))))),
          .pmlNode("mdef:Standard", attrs = c(symbId = paste0(.cond, "_obs")),
                   children = c(
                     .pmlNode("mdef:Output",
                              children = .pmlNode(
                                "ct:SymbRef",
                                attrs = c(blkIdRef = .pmlBlk[["structural"]],
                                          symbIdRef = paste(.line[["var"]])))),
                     .pmlNode("mdef:ErrorModel",
                              children = .pharmmlAssign(
                                .rxToPharmml(.pharmmlErrorModelExpr(ui, .line)))),
                     .pmlNode("mdef:ResidualError",
                              children = .pmlNode("ct:SymbRef",
                                                  attrs = c(symbIdRef = .eps))))))),
      indent = indent)
  }, character(1), USE.NAMES = FALSE)
  paste(.ret, collapse = "\n")
}

#' @export
rxUiGet.pharmmlObservationModel <- function(x, ...) {
  .pharmmlObservationModel(x[[1]])
}
attr(rxUiGet.pharmmlObservationModel, "rstudio") <- "character"
