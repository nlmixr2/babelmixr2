#' Does this model use linCmt()?
#'
#' @param ui rxode2 ui function
#' @return TRUE when the model has a solved linear compartment
#' @noRd
#' @author Matthew L. Fidler
.bblHasLinCmt <- function(ui) {
  .flags <- rxode2::rxModelVars(ui)$flags
  if (!any(names(.flags) == "linCmtFlg")) return(FALSE)
  !identical(as.integer(.flags[["linCmtFlg"]]), 0L)
}

#' Get the linCmt() micro-constants from rxode2
#'
#' `rxode2::linCmtMicro()` only exists in newer rxode2; without it the
#' closed-form translations are not available and every `linCmt()`
#' model is translated to ODEs instead.
#'
#' @param ui rxode2 ui
#' @return list of micro-constants (see `rxode2::linCmtMicro()`) or
#'   `NULL` when rxode2 cannot supply them
#' @noRd
#' @author Matthew L. Fidler
.bblLinCmtMicro <- function(ui) {
  if (!("linCmtMicro" %in% getNamespaceExports("rxode2"))) return(NULL)
  .fun <- getExportedValue("rxode2", "linCmtMicro")
  .ret <- try(.fun(ui), silent=TRUE)
  if (inherits(.ret, "try-error")) return(NULL)
  .ret
}

#' The compartment names that linToOde() gives a linCmt() model
#'
#' @param micro one micro-constant list from `rxode2::linCmtMicro()`
#' @return character vector of the ODE state names, in order
#' @noRd
#' @author Matthew L. Fidler
.bblLinCmtStates <- function(micro) {
  c(if (micro$oral0 == 1L) "depot", "central",
    if (micro$ncmt > 1L) paste0("peripheral", seq_len(micro$ncmt - 1L)))
}

#' Variables assigned by a model line
#'
#' @param expr model expression
#' @return character vector of variables assigned in `expr` (including
#'   inside `if`/`else` blocks); compartment properties like `f(depot)`
#'   and ODEs return nothing
#' @noRd
#' @author Matthew L. Fidler
.bblLinCmtLhs <- function(expr) {
  if (!is.call(expr)) return(character(0))
  .op <- as.character(expr[[1]])
  if (.op %in% c("<-", "=", "~")) {
    if (is.name(expr[[2]])) return(as.character(expr[[2]]))
    return(character(0))
  }
  unique(unlist(lapply(as.list(expr)[-1], .bblLinCmtLhs)))
}

#' Is this line an ODE, initial condition or compartment property?
#'
#' @param expr model expression
#' @return logical
#' @noRd
#' @author Matthew L. Fidler
.bblLinCmtIsCmtLine <- function(expr) {
  if (!is.call(expr)) return(FALSE)
  .op <- as.character(expr[[1]])
  if (!(.op %in% c("<-", "="))) return(FALSE)
  is.call(expr[[2]])
}

#' Find which model lines depend on the ODE states or time
#'
#' Lines that do not depend on the states (or on time) can be
#' calculated before the closed-form solution (NONMEM's `$PK` or
#' Monolix's `PK:` block); the others need the solved amounts.
#'
#' @param exprs list of model expressions
#' @param states state names
#' @return list with `dep` (logical, line depends on a state or time),
#'   `cmt` (logical, ODE/property line) and `vars` (variables that
#'   depend on a state or time)
#' @noRd
#' @author Matthew L. Fidler
.bblLinCmtStateDep <- function(exprs, states) {
  # time and the dose related variables change between records
  .timeVars <- c("t", "time", "tad", "tafd", "tlast", "tfirst", "podo",
                 "dosenum")
  .depVars <- c(states, .timeVars)
  .dep <- logical(length(exprs))
  .cmt <- logical(length(exprs))
  for (.i in seq_along(exprs)) {
    .e <- exprs[[.i]]
    if (.bblLinCmtIsCmtLine(.e)) {
      .cmt[.i] <- TRUE
      next
    }
    .lhs <- .bblLinCmtLhs(.e)
    .used <- setdiff(all.vars(.e), .lhs)
    if (any(.used %in% .depVars)) {
      .dep[.i] <- TRUE
      .depVars <- unique(c(.depVars, .lhs))
    }
  }
  list(dep=.dep, cmt=.cmt, vars=setdiff(.depVars, .timeVars))
}

#' Can this linCmt() model use the closed-form solution?
#'
#' @param ode ODE version of the model (from `rxode2::linToOde()`)
#' @param micro micro-constants from `rxode2::linCmtMicro()`
#' @return the single micro-constant list when the closed-form solution
#'   can be used, otherwise `NULL`
#' @noRd
#' @author Matthew L. Fidler
.bblLinCmtNativeMicro <- function(ode, micro) {
  if (length(micro) != 1L) return(NULL)
  .micro <- micro[[1]]
  .states <- .bblLinCmtStates(.micro)
  # other ODEs (like an effect compartment) need a general ODE solver
  if (!identical(rxode2::rxState(ode), .states)) return(NULL)
  .dep <- .bblLinCmtStateDep(ode$lstExpr, .states)
  .microVars <- unique(unlist(lapply(.micro[c("ka", "v", "k", "k12", "k21", "k13", "k31")],
                                     function(e) if (is.null(e)) character(0) else all.vars(e))))
  # closed-form solutions need parameters that are constant between records
  if (any(.microVars %in% c(.dep$vars, "t", "time", "tad", "tafd", "tlast",
                            "tfirst", "podo", "dosenum"))) return(NULL)
  .micro
}

#' Translate linCmt() models so an external program can fit them
#'
#' A `linCmt()` model is changed to its ODE equivalent with
#' `rxode2::linToOde()`.  The translation for the external program is
#' then based on this ODE model; when `native` is `TRUE` and the model
#' is a pure linear compartment model, the micro-constants are returned
#' so the translator can replace the ODEs with the program's
#' closed-form solutions.
#'
#' @param env nlmixr2 estimation environment; `env$ui` is replaced by
#'   the ODE model
#' @param software name of the software for the message
#' @param native should the closed-form micro-constants be returned
#'   when possible
#' @return micro-constant list for the closed-form solution, or `NULL`
#'   (for ODE solving or when the model does not have `linCmt()`)
#' @noRd
#' @author Matthew L. Fidler
.bblLinCmtToOde <- function(env, software, native=TRUE) {
  .ui <- rxode2::rxUiDecompress(env$ui)
  if (!.bblHasLinCmt(.ui)) return(NULL)
  .micro <- NULL
  if (native) .micro <- .bblLinCmtMicro(.ui)
  .ode <- rxode2::rxUiDecompress(suppressMessages(rxode2::linToOde(.ui)))
  # older rxode2 cannot translate every linCmt() model (like a
  # `linCmt() ~ add(add.sd)` endpoint) to ODEs
  if (!isFALSE(try(.bblHasLinCmt(.ode), silent=TRUE))) {
    stop("rxode2 could not translate this linCmt() model to ODEs for ", software,
         "; update rxode2, or write the endpoint as `cp <- linCmt()` and `cp ~ ...`",
         call.=FALSE)
  }
  if (!is.null(.micro)) .micro <- .bblLinCmtNativeMicro(.ode, .micro)
  # keep what nlmixr2 attached to the ui (like the model name used for
  # the output files)
  for (.v in c("modelName", "boundedTransforms")) {
    if (exists(.v, envir=.ui, inherits=FALSE)) {
      assign(.v, get(.v, envir=.ui, inherits=FALSE), envir=.ode)
    }
  }
  env$ui <- .ode
  if (is.null(.micro)) {
    .minfo(paste0("translating linCmt() to ODEs for ", software))
  } else {
    .minfo(paste0("translating linCmt() to ", software, "'s closed-form solution"))
  }
  .micro
}

#' Get the linCmt() option from a (possibly incomplete) control
#'
#' @param control control object or list (may be `NULL`)
#' @param default default value
#' @return option value
#' @noRd
#' @author Matthew L. Fidler
.bblLinCmtControl <- function(control, default) {
  if (is.list(control) && !is.null(control$linCmt)) {
    return(control$linCmt[1])
  }
  default
}

#' Stop when a linCmt() model is translated directly
#'
#' `nlmixr2(..., est="nonmem")` and `nlmixr2(..., est="monolix")`
#' translate a `linCmt()` model to ODEs (or the closed form) first; the
#' ui properties like `$nonmemModel` need that translated model.
#'
#' @param ui rxode2 ui
#' @param est estimation method for the message
#' @return nothing, called for the error
#' @noRd
#' @author Matthew L. Fidler
.bblLinCmtAssertOde <- function(ui, est) {
  if (.bblHasLinCmt(ui)) {
    stop("linCmt() models are translated when fit with nlmixr2(..., est=\"", est,
         "\"); to translate the model directly, first convert it with rxode2::linToOde()",
         call.=FALSE)
  }
  invisible()
}
