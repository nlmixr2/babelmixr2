#' Convert an nlmixr2 model to a NONMEM control stream or Monolix project
#'
#' These write the files that `nlmixr2(model, data, est="nonmem")` or
#' `nlmixr2(model, data, est="monolix")` would write, without running
#' NONMEM or Monolix.  The file names are the same as the ones from
#' `nonmemControl(modelName=)` and `monolixControl(modelName=)`:
#'
#' - NONMEM: `<modelName>-nonmem/<modelName>.nmctl` and
#'   `<modelName>-nonmem/<modelName>.csv`
#'
#' - Monolix: `<modelName>-monolix.mlxtran`, `<modelName>-monolix.txt`
#'   and `<modelName>-monolix.csv`
#'
#' Part of the translation depends on the data: which covariates are
#' time-varying (these are not mu-referenced), which compartments are
#' dosed and how (the Monolix administrations, `RATE`), and if the data
#' has steady state (`SS`/`II`) or censoring (`CENS`/`LIMIT`) columns.
#' When `data` is given, all of these come from the data.  Otherwise a
#' one-subject template dataset is made from the model and the
#' `timeVarying`, `dose`, `doseType`, `ss`, `cens` and `limit`
#' arguments.  This template dataset is written as the `.csv` file so
#' you can see the columns the model expects; replace it with your
#' data before running NONMEM or Monolix.
#'
#' @param x nlmixr2 model function, `rxode2` ui or fit
#' @param data optional data to translate with.  When `NULL` and `x` is
#'   a fit, the fit's data are used; otherwise a template dataset is
#'   made from the arguments below.
#' @param modelName the base name of the files.  When `NULL`, the
#'   control's `modelName` is used, then the name of the model
#'   variable.
#' @param timeVarying character vector of the covariates that change
#'   with time (by default none do).  Time-varying covariates are not
#'   mu-referenced (and are Monolix regressors).
#' @param dose the compartments receiving doses (names or numbers); by
#'   default the first compartment of the model
#' @param doseType the type of dose for each `dose` compartment
#'   (recycled): `"bolus"`, `"infusion"`, `"modelRate"` (a modeled
#'   `rate()`) or `"modelDur"` (a modeled `dur()`)
#' @param ss logical; does the data have steady state doses (adds `SS`
#'   and `II`)
#' @param cens logical; does the data have censored observations (adds
#'   `CENS`)
#' @param limit logical; does the data have censoring limits (adds
#'   `LIMIT`, and implies `cens=TRUE`)
#' @param control the `nonmemControl()` or `monolixControl()` options
#'   for the translation; `runCommand` is ignored since the program is
#'   not run
#' @param ... other arguments (currently ignored)
#' @return invisibly, a named character vector of the written files
#' @export
#' @author Matthew L. Fidler
#' @examples
#' \donttest{
#' one.cmt <- function() {
#'   ini({
#'     tka <- 0.45
#'     tcl <- log(c(0, 2.7, 100))
#'     tv <- 3.45
#'     cl.wt <- 0.1
#'     eta.ka ~ 0.6
#'     eta.cl ~ 0.3
#'     eta.v ~ 0.1
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl + eta.cl + cl.wt * WT)
#'     v <- exp(tv + eta.v)
#'     linCmt() ~ add(add.sd)
#'   })
#' }
#'
#' withr::with_tempdir({
#'   print(as.nonmem(one.cmt))
#'   print(as.monolix(one.cmt, timeVarying = "WT"))
#' })
#' }
as.nonmem <- function(x, data = NULL, modelName = NULL, timeVarying = character(0),
                      dose = NULL, doseType = "bolus", ss = FALSE, cens = FALSE,
                      limit = FALSE, control = nonmemControl(), ...) {
  if (is.null(modelName) && is.null(control$modelName)) modelName <- .bblAsToolName(substitute(x))
  checkmate::assertClass(control, "nonmemControl")
  .bblAsToolRun(
    x, data = data, modelName = modelName, timeVarying = timeVarying,
    dose = dose, doseType = doseType, ss = ss, cens = cens, limit = limit,
    control = control, est = "nonmem"
  )
}

#' @rdname as.nonmem
#' @export
as.monolix <- function(x, data = NULL, modelName = NULL, timeVarying = character(0),
                       dose = NULL, doseType = "bolus", ss = FALSE, cens = FALSE,
                       limit = FALSE, control = monolixControl(), ...) {
  if (is.null(modelName) && is.null(control$modelName)) modelName <- .bblAsToolName(substitute(x))
  checkmate::assertClass(control, "monolixControl")
  .bblAsToolRun(
    x, data = data, modelName = modelName, timeVarying = timeVarying,
    dose = dose, doseType = doseType, ss = ss, cens = cens, limit = limit,
    control = control, est = "monolix"
  )
}

#' The model name from the variable passed to as.nonmem()/as.monolix()
#'
#' @param x the unevaluated `x` argument
#' @return the variable name, or `NULL` when `x` is not a variable
#' @noRd
.bblAsToolName <- function(x) {
  if (is.name(x)) {
    return(as.character(x))
  }
  NULL
}

#' Translate a model with nlmixr2() without running the external program
#'
#' @inheritParams as.nonmem
#' @param est the estimation method ("nonmem" or "monolix")
#' @return invisibly, the named paths of the written files
#' @noRd
.bblAsToolRun <- function(x, data, modelName, timeVarying, dose, doseType,
                          ss, cens, limit, control, est) {
  if (!is.null(modelName)) {
    checkmate::assertCharacter(modelName, len = 1, any.missing = FALSE, min.chars = 1)
    control$modelName <- modelName
  }
  if (is.null(data) && inherits(x, "nlmixr2FitCore")) {
    data <- x$origData
  }
  .ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(x))
  if (is.null(data)) {
    data <- .bblAsToolData(.ui, timeVarying = timeVarying, dose = dose,
                           doseType = doseType, ss = ss, cens = cens, limit = limit)
  } else if (length(timeVarying) > 0 || !is.null(dose) || !identical(doseType, "bolus") ||
               !isFALSE(ss) || !isFALSE(cens) || !isFALSE(limit)) {
    warning("'data' is given, so 'timeVarying', 'dose', 'doseType', 'ss', 'cens' and 'limit' are ignored",
            call. = FALSE)
  }
  control$runCommand <- NA
  .ret <- nlmixr2est::nlmixr2(.ui, data, est = est, control = control)
  # the translation records the files it wrote (their names include
  # the number given to a changed model)
  .files <- NULL
  if (inherits(.ret, "rxUi")) {
    .ret <- rxode2::rxUiDecompress(.ret)
    if (exists(".bblExport", envir = .ret)) {
      .files <- get(".bblExport", envir = .ret)
      rm(list = ".bblExport", envir = .ret)
    }
  }
  if (is.null(.files) && exists(".bblExport", envir = .ui)) {
    .files <- get(".bblExport", envir = .ui)
    rm(list = ".bblExport", envir = .ui)
  }
  if (is.null(.files)) {
    stop("the ", est, " translation did not finish", call. = FALSE)
  }
  invisible(.files)
}

#' Make a one-subject template dataset for translating a model
#'
#' @inheritParams as.nonmem
#' @param ui rxode2 ui
#' @return data frame with one subject that has the doses,
#'   observations and covariates the translation needs
#' @noRd
.bblAsToolData <- function(ui, timeVarying = character(0), dose = NULL, doseType = "bolus",
                           ss = FALSE, cens = FALSE, limit = FALSE) {
  checkmate::assertLogical(ss, len = 1, any.missing = FALSE)
  checkmate::assertLogical(cens, len = 1, any.missing = FALSE)
  checkmate::assertLogical(limit, len = 1, any.missing = FALSE)
  .covs <- ui$allCovs
  checkmate::assertCharacter(timeVarying, any.missing = FALSE, unique = TRUE)
  # covariates are matched ignoring case, like rxode2
  .badTv <- timeVarying[!(tolower(timeVarying) %in% tolower(.covs))]
  if (length(.badTv) > 0) {
    stop("'timeVarying' has variables that are not covariates in the model: ",
         paste(.badTv, collapse = ", "),
         if (length(.covs) > 0) paste0("\n  the covariates are: ", paste(.covs, collapse = ", ")),
         call. = FALSE)
  }
  .state <- ui$state
  if (is.null(dose)) {
    if (length(.state) == 0) {
      stop("the model has no compartments to dose", call. = FALSE)
    }
    dose <- .state[1]
  }
  if (is.numeric(dose)) {
    checkmate::assertIntegerish(dose, lower = 1, upper = length(.state),
                                any.missing = FALSE, min.len = 1, unique = TRUE)
    dose <- .state[dose]
  }
  checkmate::assertCharacter(dose, any.missing = FALSE, min.len = 1, unique = TRUE)
  .badDose <- dose[!(dose %in% .state)]
  if (length(.badDose) > 0) {
    stop("'dose' has compartments that are not in the model: ",
         paste(.badDose, collapse = ", "),
         "\n  the compartments are: ", paste(.state, collapse = ", "),
         call. = FALSE)
  }
  checkmate::assertCharacter(doseType, any.missing = FALSE, min.len = 1)
  doseType <- vapply(doseType, match.arg, character(1),
                     choices = c("bolus", "infusion", "modelRate", "modelDur"),
                     USE.NAMES = FALSE)
  if (length(doseType) != 1L && length(doseType) != length(dose)) {
    stop("'doseType' must be length 1 or the same length as 'dose'", call. = FALSE)
  }
  doseType <- rep(doseType, length.out = length(dose))
  .rate <- c(bolus = 0, infusion = 10, modelRate = -1, modelDur = -2)[doseType]
  .doses <- data.frame(id = 1L, time = 0, amt = 100, rate = unname(.rate),
                       evid = 1L, cmt = dose, dv = NA_real_)
  if (ss) {
    .doses$ss <- 1L
    .doses$ii <- 24
  }
  .cond <- ui$predDf$cond
  .obs <- data.frame(id = 1L, time = rep(c(1, 2), each = length(.cond)),
                     amt = 0, rate = 0, evid = 0L,
                     cmt = rep(.cond, 2), dv = 1)
  if (ss) {
    .obs$ss <- 0L
    .obs$ii <- 0
  }
  .ret <- rbind(.doses, .obs)
  if (cens || limit) {
    # the first observation is below the limit of quantitation
    .ret$cens <- 0L
    .ret$cens[nrow(.doses) + 1L] <- 1L
  }
  if (limit) {
    .ret$limit <- -Inf
    .ret$limit[nrow(.doses) + 1L] <- 0
  }
  for (.cov in .covs) {
    if (tolower(.cov) %in% tolower(timeVarying)) {
      .ret[[.cov]] <- seq_len(nrow(.ret))
    } else {
      .ret[[.cov]] <- 1
    }
  }
  .ret
}
