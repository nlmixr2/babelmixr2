#' Monolix get compartment information
#'
#' @param ui rxode2 user interface
#'
#' @return the compartment specification for defined compartments in
#'   adm dataset
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.monolixGetCompartmentInformation <- function(ui) {
  .adm <- .monolixGetAdm(ui)
  .cmts <- sort(unique(.adm$cmt))
  .state <- rxode2::rxState(ui)[.cmts]
  paste(paste0("compartment(cmt=", .cmts, ", amount=", .state, ")"), collapse="\n")
}
#' Monolix get PK macros for I
#'
#' @param i integer for adm to process
#'
#' @param adm adm data frame
#'
#' @param state A character vector of the states
#'
#' @return macro for adm id i
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.monolixGetPkMacrosForI <- function(i, adm, state) {
  .adm <- adm[i, ]
  .type <- paste(.adm$type)
  if (.type %in% c("bolus", "infusion")) {
    paste0("depot(type=", .adm$adm, ", target=", state[.adm$cmt],
           ", Tlag=", ifelse(is.na(.adm$lag), "0", .adm$lag), ", p=",
           ifelse(is.na(.adm$f), "1", .adm$f), ")")
  } else if (.type == "modelRate") {
    paste0("depot(type=", .adm$adm, ", target=", state[.adm$cmt],
           ", Tk0=amtDose/", .adm$rate,
           ", Tlag=", ifelse(is.na(.adm$lag), "0", .adm$lag), ", p=",
           ifelse(is.na(.adm$f), "1", .adm$f), ")")
  } else if (.type == "modelDur") {
    paste0("depot(type=", .adm$adm, ", target=", state[.adm$cmt],
           ", Tk0=", .adm$dur,
           ", Tlag=", ifelse(is.na(.adm$lag), "0", .adm$lag), ", p=",
           ifelse(is.na(.adm$f), "1", .adm$f), ")")
  } else if (.type == "empty") {
    paste0("empty(adm=", .adm$adm, ", target=", state[.adm$cmt], ")")
  }
}
#' Get the PK macros
#'
#' @param ui rxode2 ui
#' @return Monolix macros for dosing
#' @author Matthew L. Fidler
#' @noRd
.monolixGetPkMacros <- function(ui) {
  .adm <- .monolixGetAdm(ui)
  .state <- rxode2::rxState(ui)
  paste(vapply(seq_along(.adm$adm),
               .monolixGetPkMacrosForI,
               character(1), adm=.adm, state=.state,
               USE.NAMES=FALSE), collapse="\n")
}

#' @export
rxUiGet.monolixModel <- function(x, ...) {
  .ui <- x[[1]]
  .bblLinCmtAssertOde(.ui, "monolix")
  .split <- .ui$getSplitMuModel
  assignInMyNamespace(".monolixResponses", NULL)
  .lstExpr <- .split$modelWithDrop
  .pkmodel <- .monolixPkModel(.ui, .lstExpr)
  if (!is.null(.pkmodel)) .lstExpr <- .pkmodel$lstExpr
  # first drop the error lines
  .mainModel <- rxode2::rxCombineErrorLines(.ui,
                    errLines=nmGetDistributionMonolixLines(.ui),
                    paramsLine=NA,
                    modelVars=TRUE,
                    cmtLines=FALSE,
                    dvidLine=FALSE,
                    lstExpr=.lstExpr,
                    useIf=FALSE)
  .norm <- rxode2::rxNorm(eval(.mainModel))
  .mv <- rxode2::rxModelVars(.ui)
  .mod <- rxToMonolix(.norm, ui=.ui)
  if (is.null(.pkmodel)) {
    .pk <- paste0("\n\nPK:\n; Define compartments with administrations\n",
                  .monolixGetCompartmentInformation(.ui),
                  "\n; Define PK macros\n",
                  .monolixGetPkMacros(.ui))
  } else {
    .pk <- ""
    .mod <- .monolixPkModelInsert(.mod, .pkmodel, .ui)
  }
  .txtFile <- rxUiGet.monolixModelFileName(x, ...)
  .regress <- .ui$allCovs
  .cov <- .ui$saemMuRefCovariateDataFrame
  # mu-referenced covariates enter through [INDIVIDUAL]; only the rest
  # are time-varying regressors of the structural model
  .regress <- .regress[!(.regress %in% .cov$covariate)]
  .regressors <- ""
  if (length(.regress) > 0) {
    .regressors <- paste0("\n", paste(paste0(.regress, "= {use=regressor}"), collapse="\n"))
  }
  paste0("DESCRIPTION:\n",
         paste0("model translated from `babelmixr2` and `nlmixr2` function ", .ui$modelName, " to ", .txtFile, "\n\n"),
         "[LONGITUDINAL]\n",
         "input={",
         paste(setNames(c(.monolixMuRef(.ui), .regress), NULL), collapse = ","),
         "}",
         .regressors,
          ifelse(rxode2::rxGetControl(.ui, "stiff", FALSE), "\n\nodeType = stiff", ""),
         .pk,
         "\n\nEQUATION:\n", .mod,
         "\n\nOUTPUT:\noutput={",
         paste(.monolixResponses, collapse=", "), "}\n")
}
attr(rxUiGet.monolixModel, "rstudio") <- ""

.monolixPkModelVar <- "rx_cc"

#' Is this line a compartment property of a state?
#'
#' @param expr model expression
#' @return list with `fun` (the property, like `"f"` or `"alag"`, or
#'   `"ode"` for an ODE) and `state`, or `NULL` when the line is not a
#'   compartment property
#' @noRd
#' @author Matthew L. Fidler
.monolixPkModelProp <- function(expr) {
  if (!.bblLinCmtIsCmtLine(expr)) return(NULL)
  .lhs <- expr[[2]]
  if (identical(.lhs[[1]], quote(`/`))) return(list(fun="ode", state=""))
  .state <- if (length(.lhs) >= 2L) deparse1(.lhs[[2]]) else ""
  list(fun=as.character(.lhs[[1]]), state=.state)
}

#' Write a linCmt() model with Monolix's pkmodel()
#'
#' The model is the ODE version of a pure `linCmt()` model (from
#' `rxode2::linToOde()`).  The ODEs are dropped, the concentration
#' `central/(v)` (also inside an expression) becomes the output of
#' `pkmodel()` and the lag time
#' and bioavailability of the dosing compartment become `pkmodel()`'s
#' `Tlag` and `p`.  The model keeps the same number of lines (dropped
#' lines become `_drop`) so the error lines stay in place.
#'
#' @param ui rxode2 ui with the micro-constants in the `.linCmtMicro`
#'   control value
#' @param lstExpr model lines
#' @return `NULL` when the model needs ODEs (for instance a modeled
#'   rate, doses into more than one compartment, or the amounts are
#'   used), otherwise a list with the model lines (`lstExpr`), the
#'   lines to calculate before `pkmodel()` (`pre`) and the `pkmodel()`
#'   call (`pkmodel`)
#' @noRd
#' @author Matthew L. Fidler
.monolixPkModel <- function(ui, lstExpr) {
  .micro <- rxode2::rxGetControl(ui, ".linCmtMicro", NULL)
  if (is.null(.micro)) return(NULL)
  .states <- .bblLinCmtStates(.micro)
  .doseState <- .states[1]
  .oral <- .micro$oral0 == 1L
  # pkmodel() has a single administration: all doses go to the depot
  # (oral) or central compartment (bolus or infusion)
  .adm <- .monolixGetAdm(ui)
  if (!all(.adm$cmt == 1L) || length(unique(.adm$adm)) > 1L) return(NULL)
  .types <- if (.oral) "bolus" else c("bolus", "infusion")
  if (!all(paste(.adm$type) %in% .types)) return(NULL)
  .conc <- call("/", quote(central), call("(", .micro$v))
  .pre <- list()
  .args <- NULL
  .found <- FALSE
  .add <- function(name, value) {
    .pre[[length(.pre) + 1L]] <<- call("<-", str2lang(name), value)
  }
  for (.i in seq_along(lstExpr)) {
    .e <- lstExpr[[.i]]
    .prop <- .monolixPkModelProp(.e)
    if (!is.null(.prop)) {
      if (.prop$fun != "ode") {
        if (.prop$state != .doseState) return(NULL)
        if (.prop$fun %in% c("alag", "lag")) {
          .add("rx_tlag", .e[[3]])
          .args <- c(.args, "Tlag=rx_tlag")
        } else if (.prop$fun %in% c("f", "F")) {
          .add("rx_p", .e[[3]])
          .args <- c(.args, "p=rx_p")
        } else {
          # modeled rate/duration and initial conditions need ODEs
          return(NULL)
        }
      }
      lstExpr[[.i]] <- quote(`_drop`)
      next
    }
    # the concentration central/(v) is the output of pkmodel()
    .e2 <- .monolixPkModelReplace(.e, .conc, str2lang(.monolixPkModelVar))
    if (!identical(.e2, .e)) {
      .found <- TRUE
      lstExpr[[.i]] <- .e2
      .e <- .e2
    }
    # anything else that uses the amounts needs ODEs
    if (any(all.vars(.e) %in% .states)) return(NULL)
  }
  if (!.found) return(NULL)
  .add("rx_v", .micro$v)
  .add("rx_k", .micro$k)
  .pars <- c("V=rx_v", "k=rx_k")
  for (.p in c("k12", "k21", "k13", "k31")) {
    if (!is.null(.micro[[.p]])) {
      .add(paste0("rx_", .p), .micro[[.p]])
      .pars <- c(.pars, paste0(.p, "=rx_", .p))
    }
  }
  if (.oral) {
    .add("rx_ka", .micro$ka)
    .pars <- c(.pars, "ka=rx_ka")
  }
  list(lstExpr=lstExpr, pre=.pre,
       pkmodel=paste0("pkmodel(", paste(c(.pars, .args), collapse=", "), ")"))
}

#' Replace a call inside an expression
#'
#' @param expr expression
#' @param what call to replace
#' @param by replacement
#' @return expression with `what` replaced by `by`
#' @noRd
#' @author Matthew L. Fidler
.monolixPkModelReplace <- function(expr, what, by) {
  if (identical(expr, what)) return(by)
  if (!is.call(expr)) return(expr)
  as.call(lapply(as.list(expr), .monolixPkModelReplace, what=what, by=by))
}

#' Insert the pkmodel() lines into the translated Monolix equations
#'
#' @param mod translated Monolix `EQUATION:` lines
#' @param pkmodel list from `.monolixPkModel()`
#' @param ui rxode2 ui
#' @return Monolix equations with the `pkmodel()` call before its
#'   concentration is used
#' @noRd
#' @author Matthew L. Fidler
.monolixPkModelInsert <- function(mod, pkmodel, ui) {
  .pre <- vapply(pkmodel$pre, function(e) .rxToMonolix(e, ui=ui),
                 character(1), USE.NAMES=FALSE)
  .pre <- c(paste0("   ", trimws(.pre)),
            "   ; closed-form linear compartment model",
            paste0("   ", .monolixPkModelVar, " = ", pkmodel$pkmodel))
  .mod <- strsplit(mod, "\n")[[1]]
  .w <- which(grepl(paste0("\\b", .monolixPkModelVar, "\\b"), .mod))[1]
  if (is.na(.w)) {
    stop("could not place Monolix's pkmodel() in the model; use monolixControl(linCmt=\"ode\")",
         call.=FALSE)
  }
  paste(c(.mod[seq_len(.w - 1L)], .pre, .mod[seq(.w, length(.mod))]),
        collapse="\n")
}
