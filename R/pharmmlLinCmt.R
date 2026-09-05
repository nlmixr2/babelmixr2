# linCmt() -> PharmML PKmacros.
#
# The design note in pharmml-plan.md section 4.2 explains why this cannot be
# keyed on `linCmtFlg`: at the rxUi level rxode2 has not resolved the
# compartment count, so 1-, 2- and 3-compartment models all report the same
# flag.  Only the units digit (the depot flag) is usable.  The compartment
# count and the parameterisation come from the linCmt() parameter names, which
# is how rxode2 itself resolves them.

# Accepted names for the central volume, by volume "style".  rxode2 refuses to
# mix styles (see its linCmtVStyle()), and so does this.
.pharmmlLinCmtVCentral <- list(
  plain = c("v", "vc"),
  numbered = "v1"
)

# Accepted names for the elimination term.
.pharmmlLinCmtCl <- c("cl")
.pharmmlLinCmtK <- c("k", "ke", "kel")

# Peripheral parameter names, by compartment index and style.
.pharmmlLinCmtPeripheral <- list(
  # clearance style: inter-compartmental clearance + peripheral volume
  list(q = "q",  vPlain = "vp",  vNumbered = "v2", k1i = "k12", ki1 = "k21"),
  list(q = "q2", vPlain = "vp2", vNumbered = "v3", k1i = "k13", ki1 = "k31")
)

# Names that mean a parameterisation this writer deliberately does not map.
.pharmmlLinCmtUnsupported <- c("vss", "vm", "km", "alpha", "beta", "gamma",
                               "aob", "a", "b", "c", "ktr", "mtt")

#' Which linCmt() parameters a model defines
#'
#' `linCmt()` takes its arguments implicitly, from the model's own defined
#' variables.  Under the mu-referencing this writer requires, those are exactly
#' the individual parameters.
#'
#' @param ui rxode2 UI
#' @return character vector of candidate parameter names, lower-cased
#' @noRd
.pharmmlLinCmtPars <- function(ui) {
  tolower(.pharmmlParameterSymbols(ui))
}

#' Find which of `cand` is present, erroring when more than one is
#'
#' @param pars available parameter names
#' @param cand candidate names
#' @param what human-readable role, for the error message
#' @return the matched name, or NA_character_
#' @noRd
.pharmmlLinCmtPick <- function(pars, cand, what) {
  .w <- cand[cand %in% pars]
  if (length(.w) == 0L) return(NA_character_)
  if (length(.w) > 1L) {
    stop("linCmt() model defines more than one ", what, " parameter ('",
         paste(.w, collapse = "', '"), "'); PharmML translation is ambiguous",
         call. = FALSE)
  }
  .w
}

#' Resolve a linCmt() model's structure and parameterisation
#'
#' @param ui rxode2 UI
#'
#' @return a list with `ncmt` (1-3), `depot` (logical), `elimination`
#'   (`"cl"` or `"k"`), `vStyle`, and the resolved parameter names.  The names
#'   are the *original* (case-preserved) model symbols, so they can be emitted
#'   as `SymbRef`s.
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.pharmmlLinCmtInfo <- function(ui) {
  .orig <- .pharmmlParameterSymbols(ui)
  .pars <- tolower(.orig)
  .asIs <- function(n) {
    if (is.na(n)) return(NA_character_)
    .orig[match(n, .pars)]
  }

  .bad <- .pharmmlLinCmtUnsupported[.pharmmlLinCmtUnsupported %in% .pars]
  if (length(.bad) > 0L) {
    stop("linCmt() parameter '", .bad[1],
         "' is not supported in PharmML translation yet", call. = FALSE)
  }

  # rxode2 already refuses to mix volume styles ("cannot mix 'Vp' and 'V#'
  # volume styles"), and it does so while building the UI -- before this
  # function can be reached.  No duplicate guard here: the style is only
  # detected, to know whether the peripheral volume is `vp` or `v2`.
  .vPlain <- .pharmmlLinCmtPick(.pars, .pharmmlLinCmtVCentral$plain, "central volume")
  .vNum <- .pharmmlLinCmtPick(.pars, .pharmmlLinCmtVCentral$numbered, "central volume")
  .vStyle <- if (!is.na(.vNum)) "numbered" else "plain"
  .v <- if (!is.na(.vNum)) .vNum else .vPlain
  if (is.na(.v)) {
    stop("cannot find the central volume of the linCmt() model; PharmML ",
         "translation needs one of '",
         paste(c(.pharmmlLinCmtVCentral$plain, .pharmmlLinCmtVCentral$numbered),
               collapse = "', '"), "'", call. = FALSE)
  }

  .cl <- .pharmmlLinCmtPick(.pars, .pharmmlLinCmtCl, "clearance")
  .k <- .pharmmlLinCmtPick(.pars, .pharmmlLinCmtK, "elimination rate")
  if (!is.na(.cl) && !is.na(.k)) {
    stop("linCmt() model defines both a clearance ('", .cl,
         "') and an elimination rate ('", .k, "')", call. = FALSE)
  }
  if (is.na(.cl) && is.na(.k)) {
    stop("cannot find the elimination term of the linCmt() model", call. = FALSE)
  }
  .elim <- if (!is.na(.cl)) "cl" else "k"

  .periph <- list()
  for (.i in seq_along(.pharmmlLinCmtPeripheral)) {
    .p <- .pharmmlLinCmtPeripheral[[.i]]
    .vName <- if (.vStyle == "numbered") .p$vNumbered else .p$vPlain
    .q <- .pharmmlLinCmtPick(.pars, .p$q, "inter-compartmental clearance")
    .pv <- .pharmmlLinCmtPick(.pars, .vName, "peripheral volume")
    .k1i <- .pharmmlLinCmtPick(.pars, .p$k1i, "peripheral rate")
    .ki1 <- .pharmmlLinCmtPick(.pars, .p$ki1, "peripheral rate")
    if (!is.na(.q) && !is.na(.pv)) {
      .periph[[length(.periph) + 1L]] <-
        list(style = "q", q = .asIs(.q), v = .asIs(.pv))
    } else if (!is.na(.k1i) && !is.na(.ki1)) {
      .periph[[length(.periph) + 1L]] <-
        list(style = "k", k1i = .asIs(.k1i), ki1 = .asIs(.ki1))
    } else if (!is.na(.q) || !is.na(.pv) || !is.na(.k1i) || !is.na(.ki1)) {
      stop("peripheral compartment ", .i,
           " of the linCmt() model is only partly specified", call. = FALSE)
    } else {
      break
    }
  }

  .ka <- .pharmmlLinCmtPick(.pars, "ka", "absorption rate")

  list(ncmt = 1L + length(.periph),
       depot = !is.na(.ka),
       elimination = .elim,
       vStyle = .vStyle,
       v = .asIs(.v),
       cl = .asIs(.cl),
       k = .asIs(.k),
       ka = .asIs(.ka),
       peripheral = .periph)
}

#' Does this model use linCmt()?
#'
#' @param ui rxode2 UI
#' @return TRUE when any endpoint is a solved-system prediction
#' @noRd
.pharmmlIsLinCmt <- function(ui) {
  .predDf <- ui$predDf
  !is.null(.predDf) && any(isTRUE(.predDf$linCmt) | .predDf$linCmt %in% TRUE)
}

#' State names rxode2 allocates for a solved model
#'
#' These match what `rxSolve()` reports for a `linCmt()` model, so the PharmML
#' compartment amounts line up with rxode2's own naming.
#'
#' @param i compartment index (1 = central)
#' @return state name
#' @noRd
.pharmmlLinCmtState <- function(i) {
  if (i == 1L) return("central")
  paste0("peripheral", i - 1L)
}

#' A PK macro `Value` element
#'
#' @param arg argument name, or `NULL` for a positional value
#' @param child already-emitted child node
#' @return character(1)
#' @noRd
.pharmmlMacroValue <- function(arg, child) {
  .attrs <- if (is.null(arg)) NULL else c(argument = arg)
  .pmlNode("mdef:Value", attrs = .attrs, children = child)
}

#' A PK macro `Value` holding a reference to a model parameter
#'
#' @param arg argument name, or `NULL` for a positional value
#' @param sym parameter name
#' @param ui rxode2 UI, used to resolve the owning block
#' @return character(1)
#' @noRd
.pharmmlMacroPar <- function(arg, sym, ui) {
  .pharmmlMacroValue(arg, .rxToPharmml(str2lang(sym), ui))
}

#' PharmML PKmacros for a linCmt() model
#'
#' Peripheral compartments are always emitted as micro-constants.  PharmML 0.9
#' gives the `Peripheral` macro no named argument for the transfer terms -- only
#' `amount`/`volume`/`concentration` are named, and the upstream `advan3`/
#' `advan4` reference files pass the rates positionally -- so a `Q`/`Vp` model
#' is converted to `k12 = Q/V`, `k21 = Q/Vp` rather than emitted as an
#' ambiguous pair of unnamed values.
#'
#' @param ui rxode2 UI
#'
#' @param indent Indent depth
#'
#' @return character(1) holding the `PKmacros` element
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.pharmmlPkMacros <- function(ui, indent = 0L) {
  .info <- .pharmmlLinCmtInfo(ui)
  .var <- paste(ui$predDf$var[1])

  .macros <- .pmlNode(
    "mdef:Compartment",
    children = c(
      .pharmmlMacroValue("cmt", .pmlText("ct:Int", 1L)),
      .pharmmlMacroValue("amount",
                         .pmlNode("ct:SymbRef",
                                  attrs = c(symbIdRef = .pharmmlLinCmtState(1L)))),
      .pharmmlMacroPar("volume", .info$v, ui),
      .pharmmlMacroValue("concentration",
                         .pmlNode("ct:SymbRef", attrs = c(symbIdRef = .var)))))

  for (.i in seq_along(.info$peripheral)) {
    .p <- .info$peripheral[[.i]]
    if (.p$style == "q") {
      # A computed rate is an expression, and a macro Value accepts only a
      # SymbRef, a scalar or a ct:Assign -- so the expression is wrapped.
      .k1i <- .pharmmlAssign(
        .rxToPharmml(bquote(.(str2lang(.p$q)) / .(str2lang(.info$v))), ui))
      .ki1 <- .pharmmlAssign(
        .rxToPharmml(bquote(.(str2lang(.p$q)) / .(str2lang(.p$v))), ui))
    } else {
      .k1i <- .rxToPharmml(str2lang(.p$k1i), ui)
      .ki1 <- .rxToPharmml(str2lang(.p$ki1), ui)
    }
    .macros <- c(
      .macros,
      .pmlNode("mdef:Peripheral",
               children = c(
                 .pharmmlMacroValue(NULL, .k1i),
                 .pharmmlMacroValue(NULL, .ki1),
                 .pharmmlMacroValue("amount",
                                    .pmlNode("ct:SymbRef",
                                             attrs = c(symbIdRef = .pharmmlLinCmtState(.i + 1L)))))))
  }

  if (.info$depot) {
    .macros <- c(.macros,
                 .pmlNode("mdef:Oral",
                          children = c(
                            .pharmmlMacroValue("adm", .pmlText("ct:Int", 1L)),
                            .pharmmlMacroValue("cmt", .pmlText("ct:Int", 1L)),
                            .pharmmlMacroPar("ka", .info$ka, ui))))
  } else {
    .macros <- c(.macros,
                 .pmlNode("mdef:IV",
                          children = c(
                            .pharmmlMacroValue("adm", .pmlText("ct:Int", 1L)),
                            .pharmmlMacroValue("cmt", .pmlText("ct:Int", 1L)))))
  }

  .elim <- .pharmmlMacroValue("cmt", .pmlText("ct:Int", 1L))
  if (.info$elimination == "cl") {
    .elim <- c(.elim,
               .pharmmlMacroPar("V", .info$v, ui),
               .pharmmlMacroPar("CL", .info$cl, ui))
  } else {
    .elim <- c(.elim, .pharmmlMacroPar("k", .info$k, ui))
  }
  .macros <- c(.macros, .pmlNode("mdef:Elimination", children = .elim))

  .pmlNode("mdef:PKmacros", children = .macros, indent = indent)
}
