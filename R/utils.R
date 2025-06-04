#' Get compartment number from name
#'
#' @param nm Name of the compartment
#' @param ui rxode2 ui to query against
#' @param error Should this error if not found (if `FALSE` return `NA_integer`)
#' @return compartment number (or possibily NA)
#' @author Bill Denney and Matthew L. Fidler
#' @noRd
.rxGetCmtNumber <- function(nm, ui, error=TRUE) {
  if (is.numeric(nm)) {
    nm
  } else {
    if (is.name(nm)) {
      nm <- as.character(nm)
    } else if (is.character(nm)) {
      checkmate::assertCharacter(nm, len=1, any.missing=FALSE)
    } else if (error) {
      stop("invalid attempt to find a compartment number: ", as.character(nm))
    } else {
      return(NA_integer_)
    }
    .cmt <- which(rxode2::rxModelVars(ui)$state %in% nm)
    if (length(.cmt) == 0) {
      if (error) {
        stop("cannot find compartment '", nm, "' in the model",
             call.=FALSE)
      } else {
        return(NA_integer_)
      }
    }
    .cmt
  }
}

#' Determine if the expression is d/dt()
#'
#' @param expr expression
#' @return boolean determining if expression is a d/dt() expression
#' @author Matthew L. Fidler
#' @noRd
.rxIsDdt <- function(expr) {
  if (length(expr) != 3) return(FALSE)
  if (!identical(expr[[1]], quote(`/`))) return(FALSE)
  if (!identical(expr[[2]], quote(`d`))) return(FALSE)
  if (length(expr[[3]]) != 2) return(FALSE)
  if (!identical(expr[[3]][[1]], quote(`dt`))) return(FALSE)
  TRUE
}
#' Is a variable known to be non-zero
#'
#' @param variable Varible name to check
#' @param ui rxode2 ui
#' @return logical to say if the variable is known to be non-zero
#' @author Matthew L. Fidler
#' @noRd
.rxIsKnownNonZeroVariable <- function(variable, ui) {
  if (!is.null(names(variable))) {
    if (.rxIsKnownNonZeroVariable(names(variable), ui)) return(TRUE)
  }
  variable <- as.character(variable)
  if (variable %in% toupper(ui$allCovs)) return(TRUE)
  .split <- ui$getSplitMuModel
  .mus <- c(.split$pureMuRef, .split$taintMuRef)
  .w <- which(variable == toupper(.mus))
  if (length(.w) != 1) return(FALSE)
  .tv <- names(.w)
  .w <- which(ui$muRefCurEval$parameter == .tv)
  if (length(.w) != 1) return(FALSE)
  .curEval <- ui$muRefCurEval$curEval[.w]
  if (.curEval == "exp") return(TRUE)
  if (any(.curEval == c("expit", "probitInv"))) {
    .low <- ui$muRefCurEval$parameter$low[.w]
    if (.low >= 0) return(TRUE)
    .hi <- ui$muRefCurEval$parameter$hi[.w]
    if (.hi <= 0) return(TRUE)
  }
  FALSE
}
#' Should the variable be protected from being zero?
#'
#' @param variable name of the variable
#' @param ui rxode2 ui
#' @return boolean of if the variable needs protection
#' @author Matthew L. Fidler
#' @noRd
.rxShouldProtectZeros <- function(variable, ui) {
  # should protect zeros if requested, not in an if/else block
  # and if the variable is known to be something non-zero
  if (!rxode2::rxGetControl(ui, "protectZeros", getOption("babelmixr2.protectZeros", TRUE))) return(FALSE)
  if (rxode2::rxGetControl(ui, ".ifelse", getOption("babelmixr2.ifelse", FALSE))) return(FALSE)
  if (.rxIsKnownNonZeroVariable(variable, ui)) return(FALSE)
  return(TRUE)
}


#' @import nlmixr2data
#' @import nlmixr2plot
#' @importFrom stats predict logLik na.fail pchisq approxfun cov cov2cor dlnorm median na.omit qchisq qnorm
#' @noRd
.genHardReExport <- function(fun) { # nocov start
  message("Writing hard reexport: ", fun)
  .args <- deparse(eval(str2lang(paste0("args(", fun, ")"))))
  if (fun == "nlmixr2est::nlmixr2") {
    .args <- gsub("data *=*[^,]*,", "data = NULL,", .args)
  }
  .args <- .args[-length(.args)]
  .formalArgs <- as.character(eval(str2lang(paste0("formalArgs(", fun, ")"))))
  .w <- which(.formalArgs == "...")
  .formalArgs <- paste0(.formalArgs, "=", .formalArgs)
  .has3 <- FALSE
  if (length(.w) > 0) {
    .formalArgs[.w] <- "..."
    .has3 <- TRUE
  }
  .formalArgs <- paste(.formalArgs, collapse=", ")
  .newFun <- strsplit(fun, "::")[[1]][2]
  .has3text <- NULL
  if (.has3) {
    .has3text <- paste0("#' @param ... Additional arguments passed to [", fun, "()].")
  }
  ret <-
    paste(
      c(paste("#' @inherit", fun),
        .has3text,
        "#' @export",
        trimws(
          deparse(str2lang(paste0(
            c(paste0(.newFun, " <- ", paste0(.args, collapse="\n"), " {"),
              paste0(fun, "(", .formalArgs, ")"),
              "}"
              ),
            collapse="\n"
          ))),
          which = "right"
        )
        ),
      collapse="\n"
    )
  ret <- gsub(x = ret, pattern = "{", replacement = "{ # nocov start", fixed = TRUE)
  ret <- gsub(x = ret, pattern = "}", replacement = "} # nocov end", fixed = TRUE)
  ret
} # nocov end

.genSoftReExport <- function(fun, alias=NULL) { # nocov start
  message("Writing soft reexport: ", fun)
  .newFun <- strsplit(fun, "::")[[1]]
  .pkg <- .newFun[1]
  .fun <- .newFun[2]
  if (is.na(.fun)) return("")
  .aliasText <- NULL
  if (!is.null(alias)) {
    .aliasText <- paste0("#' @rdname ", .fun)
  }
  paste(
    c(paste("#' @importFrom", .pkg, .fun),
      .aliasText,
      "#' @export",
      fun
      ),
    collapse="\n"
  )
} # nocov end

.genReexports <- function(soft=c(
  "lotri::lotri",
  "magrittr::`%>%`",
  "monolix2rx::mlxtran",
  "monolix2rx::monolix2rx",
  "monolix2rx::monolix2rx",
  "nlmixr2est::.nlmixrNlmeFun",
  "nlmixr2est::ACF",
  "nlmixr2est::VarCorr",
  "nlmixr2est::augPred",
  "nlmixr2est::fixed.effects",
  "nlmixr2est::fixef",
  "nlmixr2est::getData",
  "nlmixr2est::getValidNlmixrCtl",
  "nlmixr2est::getVarCov",
  "nlmixr2est::groupedData",
  "nlmixr2est::nlme",
  "nlmixr2est::nlmixr",
  "nlmixr2est::nlmixr2",
  "nlmixr2est::nlmixr2AllEst",
  "nlmixr2est::nlmixr2Est",
  "nlmixr2est::nlmixr2NlmeControl",
  "nlmixr2est::nmObjGetControl",
  "nlmixr2est::nmObjGetFoceiControl",
  "nlmixr2est::nmObjHandleControlObject",
  "nlmixr2est::pdBlocked",
  "nlmixr2est::pdCompSymm",
  "nlmixr2est::pdConstruct",
  "nlmixr2est::pdDiag",
  "nlmixr2est::pdFactor",
  "nlmixr2est::pdIdent",
  "nlmixr2est::pdLogChol",
  "nlmixr2est::pdMat",
  "nlmixr2est::pdMatrix",
  "nlmixr2est::pdNatural",
  "nlmixr2est::pdSymm",
  "nlmixr2est::random.effects",
  "nlmixr2est::ranef",
  "nlmixr2est::reStruct",
  "nlmixr2est::varComb",
  "nlmixr2est::varConstPower",
  "nlmixr2est::varExp",
  "nlmixr2est::varFixed",
  "nlmixr2est::varFunc",
  "nlmixr2est::varIdent",
  "nlmixr2est::varPower",
  "nlmixr2est::varWeights",
  "nonmem2rx::as.nonmem2rx",
  "nonmem2rx::nmcov",
  "nonmem2rx::nmcov",
  "nonmem2rx::nmext",
  "nonmem2rx::nminfo",
  "nonmem2rx::nmtab",
  "nonmem2rx::nmxml",
  "nonmem2rx::nonmem2rx",
  "rxode2::.minfo",
  "rxode2::RxODE",
  "rxode2::`ini<-`",
  "rxode2::`model<-`",
  "rxode2::add.dosing",
  "rxode2::add.sampling",
  "rxode2::et",
  "rxode2::etExpand",
  "rxode2::eventTable",
  "rxode2::expit",
  "rxode2::geom_amt",
  "rxode2::geom_cens",
  "rxode2::ini",
  "rxode2::logit",
  "rxode2::lotri",
  "rxode2::model",
  "rxode2::probit",
  "rxode2::probitInv",
  "rxode2::rxCat",
  "rxode2::rxClean",
  "rxode2::rxControl",
  "rxode2::rxDerived",
  "rxode2::rxFun",
  "rxode2::rxModelVars",
  "rxode2::rxParam",
  "rxode2::rxParams",
  "rxode2::rxSetPipingAuto",
  "rxode2::rxSetSeed",
  "rxode2::rxSolve",
  "rxode2::rxUiGet",
  "rxode2::rxode",
  "rxode2::rxode2",
  "rxode2::stat_amt",
  "rxode2::stat_cens"
),
hard=c(
  "nlmixr2plot::traceplot",
  "nlmixr2est::vpcSim",
  "nlmixr2plot::vpcPlot",
  "nlmixr2plot::vpcPlotTad",
  "nlmixr2plot::vpcCens",
  "nlmixr2plot::vpcCensTad",
  "nlmixr2est::saemControl",
  "nlmixr2est::foceiControl",
  "nlmixr2est::nlmeControl",
  "nlmixr2est::tableControl",
  "nlmixr2est::bobyqaControl",
  "nlmixr2est::lbfgsb3cControl",
  "nlmixr2est::n1qn1Control",
  "nlmixr2est::newuoaControl",
  "nlmixr2est::nlmControl",
  "nlmixr2est::nlminbControl",
  "nlmixr2est::nlsControl",
  "nlmixr2est::optimControl",
  "nlmixr2est::uobyqaControl",
  "nlmixr2est::addCwres",
  "nlmixr2est::addNpde",
  "nlmixr2est::addTable",
  "nlmixr2est::setOfv",
  "nlmixr2extra::profileFixed",
  "nlmixr2extra::profileFixedSingle",
  "nlmixr2extra::profileLlp",
  "nlmixr2extra::preconditionFit",
  "nlmixr2extra::bootstrapFit",
  "nlmixr2extra::bootplot"
)
                          ) { # nocov start
  writeLines(c("# Generated from .genReexports()\n",
               paste(vapply(soft, .genSoftReExport, character(1), USE.NAMES=FALSE),
                     collapse="\n\n")),
             devtools::package_file("R/reexports.R"))
  writeLines(c("# Generated from .genReexports()\n",
               paste(vapply(hard, .genSoftReExport, character(1), USE.NAMES=FALSE),
                     collapse="\n\n")),
               ## paste(vapply(hard, .genHardReExport, character(1), USE.NAMES=FALSE),
               ##       collapse="\n\n")),
             devtools::package_file("R/hardReexports.R"))
  ""
} # nocov end
