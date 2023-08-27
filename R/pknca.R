#' Estimate starting parameters using PKNCA
#'
#' @details Parameters are estimated as follows:
#'
#' \itemize{
#' \item{ka}{4 half-lives to Tmax but not higher than 3:  \code{log(2)/(tmax/4)}}
#' \item{vc}{Inverse of dose-normalized Cmax}
#' \item{cl}{Estimated as the median clearance}
#' \item{vp,vp2}{2- and 4-fold the \code{vc}, respectively by default,
#'   controlled by the \code{vpMult} and \code{vp2Mult} arguments to
#'   \code{pkncaControl}}
#' \item{q,q2}{0.5- and 0.25-fold the \code{cl}, respectively by default,
#'   controlled by the \code{qMult} and \code{q2Mult} arguments to
#'   \code{pkncaControl}}
#' }
#'
#' The bounds for the parameter estimates are set to 10% of the first percentile
#' and 10 times the 99th percentile.  (For ka, the lower bound is set to the
#' lower of 10% of the first percentile or 0.03 and the upper bound is not
#' modified from 10 times the 99th percentile.)
#'
#' Parameter estimation methods may be changed in a future version.
#'
#' @inheritParams nlmixr2est::nlmixr2Est
#'
#' @return A model with updated starting parameters.  In the model a new element
#'   named "nca" will be available which includes the PKNCA results used for the
#'   calculation.
#' @export
nlmixr2Est.pknca <- function(env, ...) {
  rxode2::rxReq("PKNCA")
  rxode2::rxReq("units")
  control <- env$control[[1]]
  # Get the units from the basic units (before unit conversion)
  dUnitsData <-
    PKNCA::pknca_units_table(
      concu = control$concu,
      doseu = control$doseu,
      timeu = control$timeu
    )

  if (is.null(control$ncaResults)) {
    oNCA <- calcPknca(env, pkncaUnits = dUnitsData)
  } else {
    oNCA <- control$ncaResults
  }

  unitsSetup <-
    stats::setNames(
      object = dUnitsData$PPORRESU,
      nm = dUnitsData$PPTESTCD
    )
  conversionFactors <-
    modelUnitConversion(
      dvu = control$concu,
      amtu = control$doseu,
      timeu = control$timeu,
      volumeu = control$volumeu
    )
  modelDataConversions <-
    data.frame(
      PPORRESU=
        c(
          unitsSetup[["cmax"]],
          unitsSetup[["cmax.dn"]],
          unitsSetup[["cl.last"]],
          unitsSetup[["vss.last"]]
        ),
      PPSTRESU=
        c(
          conversionFactors$cmtu,
          conversionFactors$volumeu,
          conversionFactors$clearanceu,
          conversionFactors$volumeu
        )
    )
  ncaUnitsToModelUnits <-
    merge(
      dUnitsData[dUnitsData$PPTESTCD %in% c("vss.last", "cmax", "cl.last"), ],
      modelDataConversions,
      all.x = TRUE
    )
  ncaUnitsToModelUnits$conversionFactor <- NA_real_
  for (idx in seq_len(nrow(ncaUnitsToModelUnits))) {
    if (is.na(ncaUnitsToModelUnits$PPORRESU[idx]) | is.na(ncaUnitsToModelUnits$PPSTRESU[idx])) {
      ncaUnitsToModelUnits$conversionFactor[idx] <- 1
    } else {
      ncaUnitsToModelUnits$conversionFactor[idx] <-
        units::set_units(
          units::set_units(
            1,
            ncaUnitsToModelUnits$PPORRESU[idx],
            mode = "standard"
          ),
          ncaUnitsToModelUnits$PPSTRESU[idx],
          mode = "standard"
        )
    }
  }

  unitConversions <-
    stats::setNames(
      ncaUnitsToModelUnits$conversionFactor,
      nm = ncaUnitsToModelUnits$PPTESTCD
    )

  pkncaEst <- calcPkncaEst(objectPknca = oNCA)
  paramEstimates <-
    ncaToEst(
      tmax = pkncaEst$tmax,
      cmaxdn = pkncaEst$cmaxdn,
      cl = pkncaEst$cllast,
      control = control,
      unitConversions = unitConversions
    )
  # What parameters should be modified?  And then modify them.
  murefNames <- env$ui$getSplitMuModel$pureMuRef
  updateNames <- intersect(murefNames, names(paramEstimates))
  newEnv <- do.call(ini_transform, append(list(x=env$ui), paramEstimates[updateNames]))

  if (unitConversions[["cmax"]] != 1) {
    # No need to bother with modifications if the unit conversion is unity
    dvParam <- control$dvParam
    if (is.null(dvParam)) {
      dvResid <- getDvLines(env$ui$fun)
      if (length(dvResid) == 1) {
        dvParam <- dvResid[[1]][[2]]
        if (is.name(dvParam)) {
          cli::cli_abort("Could not detect the dependent variable (not a name), use pkncaControl(dvParam) to fix")
        }
      } else {
        cli::cli_abort("Could not detect the dependent variable (no specific line found), use pkncaControl(dvParam) to fix")
      }
    }
    dvAssign <- getDvLines(modelfun = env$ui$fun, dvAssign = dvParam)
    if (length(dvAssign) != 1) {
      cli::cli_abort("Could not detect DV assignment for unit conversion")
    }
    newEnv <-
      eval(str2lang(
        sprintf(
          "rxode2::model(newEnv, %s <- %g*%s)",
          dvParam,
          1/unitConversions[["cmax"]],
          deparse1(dvAssign[[1]][[3]])
        )
      ))
  }

  env$ui <- newEnv
  env$nca <- oNCA

  class(env) <- "pkncaEst"
  env
}

#' Get the lines of the model which assign the dependent variable(s)
#' @return A list of dependent variable lines (or NULL if there are no dependent
#'   variables)
#' @noRd
getDvLines <- function(modelfun, inModel = FALSE, dvAssign = NULL) {
  if (!is.null(dvAssign)) {
    if (is.character(dvAssign)) {
      dvAssign <- lapply(X = dvAssign, FUN = as.name)
    } else if (is.name(dvAssign)) {
      dvAssign <- list(dvAssign)
    }
    allAreNames <- vapply(X = dvAssign, FUN = is.name, FUN.VALUE = TRUE)
    if (!allAreNames) {
      cli::cli_abort("dvAssign must be a name, a character string, or a list of names")
    }
  }
  if (is.function(modelfun)) {
    ret <- getDvLines(methods::functionBody(modelfun), inModel = inModel, dvAssign = dvAssign)
  } else if (inherits(modelfun, "{")) {
    ret <- lapply(X = modelfun, FUN = getDvLines, inModel = inModel, dvAssign = dvAssign)
  } else if (is.name(modelfun)) {
    ret <- NULL
  } else if (inherits(modelfun, "<-") | inherits(modelfun, "=")) {
    if (inModel & !is.null(dvAssign)) {
      if (any(vapply(X = dvAssign, FUN = identical, FUN.VALUE = TRUE, y = modelfun[[2]]))) {
        # Return the DV assignment line(s) for the DV of interest
        ret <- list(modelfun)
      } else {
        ret <- NULL
      }
    } else {
      ret <- NULL
    }
  } else if (is.call(modelfun)) {
    if (identical(modelfun[[1]], as.name("~")) & inModel & is.null(dvAssign)) {
      # Return the residual error lines
      ret <- list(modelfun)
    } else if (identical(modelfun[[1]], as.name("ini"))) {
      # The ini block doesn't have the DV lines
      ret <- NULL
    } else if (identical(modelfun[[1]], as.name("model"))) {
      # This is what we want
      ret <- lapply(X = modelfun, FUN = getDvLines, inModel = TRUE, dvAssign = dvAssign)
    } else {
      # No other call has information that we want (I think), but recurse in
      # case the model is within the other call.
      ret <- lapply(X = modelfun, FUN = getDvLines, inModel = inModel, dvAssign = dvAssign)
    }
  } else {
    cli::cli_abort("Error finding DV lines, please report a bug") # nocov
  }
  if (is.list(ret)) {
    # The list of calls will still be a list, but it will not be nested.
    ret <- unlist(ret, use.names = FALSE)
  }
  ret
}

#' Perform NCA calculations with PKNCA
#'
#' @inheritParams nlmixr2est::nlmixr2Est
#' @return A PKNCAresults object
#' @noRd
calcPknca <- function(env, pkncaUnits) {
  # Normalize column names
  rxControl <- env$control[[1]]$rxControl
  control <- env$control[[1]]
  rawData <- env$data

  if (!is.null(control$ncaData)) {
    # as.data.frame() due to https://github.com/nlmixr2/nlmixr2est/pull/262
    rawData <- as.data.frame(control$ncaData)
  }
  cleanData <- bblDatToPknca(model = env$ui, data = rawData, rxControl=rxControl)
  cleanColNames <- getStandardColNames(cleanData$obs)
  oConcFormula <-
    stats::as.formula(sprintf(
      "%s~%s|%s",
      cleanColNames[["dv"]], cleanColNames[["time"]],
      paste(c(control$groups, cleanColNames[["id"]]), collapse="+")
    ))
  oDoseFormula <-
    stats::as.formula(sprintf(
      "%s~%s|%s",
      cleanColNames[["amt"]], cleanColNames[["time"]],
      paste(c(control$groups, cleanColNames[["id"]]), collapse="+")
    ))
  # Determine route of administration
  obsCmt <- unique(cleanData$obs[[cleanColNames[["cmt"]]]])
  doseCmt <- unique(cleanData$dose[[cleanColNames[["cmt"]]]])
  doseRoute <- ifelse(obsCmt == doseCmt, yes = "intravascular", no = "extravascular")

  oConc <- PKNCA::PKNCAconc(data = cleanData$obs, oConcFormula, sparse = control$sparse)
  oDose <- PKNCA::PKNCAdose(data = cleanData$dose, oDoseFormula, route = doseRoute)

  oData <- PKNCA::PKNCAdata(oConc, oDose, units = pkncaUnits)
  intervals <- oData$intervals
  intervals$cl.last <- intervals$auclast
  intervals$cmax.dn <- intervals$cmax
  intervals$vss.last <- intervals$auclast
  oData$intervals <- intervals
  oNCA <- PKNCA::pk.nca(oData)
  oNCA
}

#' Extract desired PKNCAresults to a list and set bounds
#' @param objectPknca A PKNCAresults object with at least tmax, cmax.dn, and cl.last calculated
#' @return A list with named values for tmax, cmaxdn, cllast
#' @noRd
calcPkncaEst <- function(objectPknca) {
  ncaParams <- as.data.frame(objectPknca)
  # one compartment parameters including unit conversion
  tmaxValues <-
    c(0.1, 1, 10)*
    stats::quantile(ncaParams$PPORRES[ncaParams$PPTESTCD == "tmax"], probs = c(0.01, 0.5, 0.99), na.rm=TRUE, names = FALSE)
  cmaxdnValues <-
    c(0.1, 1, 10)*
    stats::quantile(ncaParams$PPORRES[ncaParams$PPTESTCD == "cmax.dn"], probs = c(0.01, 0.5, 0.99), na.rm=TRUE, names = FALSE)
  cllastValues <-
    c(0.1, 1, 10)*
    stats::quantile(ncaParams$PPORRES[ncaParams$PPTESTCD == "cl.last"], probs = c(0.01, 0.5, 0.99), na.rm=TRUE, names = FALSE)

  # Ensure that NCA as sufficiently successful
  naValues <- character()
  if (any(is.na(tmaxValues))) {
    naValues <- c(naValues, "tmax")
  }
  if (any(is.na(cmaxdnValues))) {
    naValues <- c(naValues, "cmax.dn")
  }
  if (any(is.na(cllastValues))) {
    naValues <- c(naValues, "cl.last")
  }
  if (length(naValues) > 0) {
    cli::cli_abort(paste(
      "All",
      paste(naValues, collapse = ", "),
      "values were NA for NCA, cannot proceed with PKNCA estimation of parameters"
    ))
  }
  list(
    tmax = tmaxValues,
    cmaxdn = cmaxdnValues,
    cllast = cllastValues
  )
}

#' Convert NCA parameters to compartmental parameter values
#' @noRd
ncaToEst <- function(tmax, cmaxdn, cl, control, unitConversions) {
  ncaEstimates <-
    list(
      ka=
        # 4 absorption half-lives
        sort(log(2)/(tmax/4)),
      vc=
        sort(unitConversions[["vss.last"]] / cmaxdn),
      cl=
        unitConversions[["cl.last"]] * cl
    )
  ncaEstimates$ka <-
    pmin(
      c(0.03, 3, Inf),
      ncaEstimates$ka
    )
  # two compartment parameters
  ncaEstimates$vp <- ncaEstimates$vc*control$vpMult
  ncaEstimates$q <-  ncaEstimates$cl*control$qMult
  # three compartment parameters
  ncaEstimates$vp2 <-  ncaEstimates$vc*control$vp2Mult
  ncaEstimates$q2 <-  ncaEstimates$cl*control$q2Mult
  ncaEstimates
}

# Update the ini() with parameters given by ...; transformations to the
# estimation scale are automatically applied
ini_transform <- function(x, ..., envir = parent.frame()) {
  changeArgs <- list(...)
  # This only works for fixed effects, so formula are not allowed
  checkmate::assert_names(names(changeArgs))
  murefNames <- x$getSplitMuModel$pureMuRef
  murefTrans <- x$muRefCurEval
  inverseTrans <-
    list(
      exp=log,
      logit=rxode2::expit
      # TODO: add all of the other transforms here
    )

  for (nm in names(changeArgs)) {
    if (nm %in% names(murefNames)) {
      # It is already the transformed parameter, no modification required
      x <- do.call(rxode2::ini, append(list(x=x), changeArgs[nm]))
    } else if (nm %in% murefNames) {
      iniName <- names(murefNames)[murefNames == nm]
      currentTrans <- murefTrans$curEval[murefTrans$parameter == iniName]
      if (currentTrans == "") {
        # No transformation
        transFun <- identity
      } else {
        transFun <- inverseTrans[[currentTrans]]
      }
      if (is.null(transFun)) {
        cli::cli_abort(paste("cannot invert the transform (please report a bug):", transFun)) # nocov
      }
      x <-
        do.call(
          rxode2::ini,
          append(
            list(x=x),
            stats::setNames(
              list(transFun(changeArgs[[nm]])),
              iniName
            )
          )
        )
      newValue <- changeArgs[[nm]]
    }
  }
  x
}

#' PKNCA estimation control
#'
#' @inheritParams PKNCA::PKNCAconc
#' @param concu,doseu,timeu concentration, dose, and time units from the source
#'   data (passed to \code{PKNCA::pknca_units_table()}).
#' @param volumeu compartment volume for the model (if \code{NULL}, simplified
#'   units from source data will be used)
#' @param vpMult,qMult,vp2Mult,q2Mult Multipliers for vc and cl to provide
#'   initial estimates for vp, q, vp2, and q2
#' @param dvParam The parameter name in the model that should be modified for
#'   concentration unit conversions.  It must be assigned on a line by itself,
#'   separate from the residual error model line.
#' @param groups Grouping columns for NCA summaries by group (required if
#'   \code{sparse = TRUE})
#' @param ncaData Data to use for calculating NCA parameters.  Typical use is
#'   when a subset of the original data are informative for NCA.
#' @param ncaResults Already computed NCA results (a PKNCAresults object) to
#'   bypass automatic calculations.  At least the following parameters must be
#'   calculated in the NCA: tmax, cmax.dn, cl.last
#' @return A list of parameters
#' @export
pkncaControl <- function(concu = NA_character_, doseu = NA_character_, timeu = NA_character_,
                         volumeu = NA_character_,
                         vpMult=2, qMult=1/2,
                         vp2Mult=4, q2Mult=1/4,
                         dvParam = "cp",
                         groups = character(),
                         sparse = FALSE,
                         ncaData = NULL,
                         ncaResults = NULL,
                         rxControl=rxode2::rxControl()) {
  getValidNlmixrCtl.pknca(
    list(
      concu = concu,
      doseu = doseu,
      timeu = timeu,
      volumeu = volumeu,
      vpMult = vpMult,
      qMult = qMult,
      vp2Mult = vp2Mult,
      q2Mult = q2Mult,
      dvParam = dvParam,
      groups = groups,
      sparse = sparse,
      ncaData = ncaData,
      ncaResults = ncaResults,
      rxControl=rxControl
    )
  )
}

#' @export
getValidNlmixrCtl.pknca <- function(control) {
  orig <- control
  if (inherits(control, "getValidNlmixrControl")) {
    if (is.null(orig[[1]])) {
      # Use default values
      orig[[1]] <- pkncaControl()
    } else if (is.list(orig[[1]]) && length(orig[[1]]) == 0) {
      # Use default values
      orig[[1]] <- pkncaControl()
    }
    control <- orig[[1]]
  }

  checkmate::assert_names(
    x = names(control),
    permutation.of = names(formals(pkncaControl))
  )
  # verify units look like units
  for (unitNm in c("concu", "doseu", "timeu", "volumeu")) {
    checkmate::assert_character(control[[unitNm]], .var.name = unitNm, null.ok = FALSE, len = 1, min.chars = 1)
  }
  # Verify that multipliers are numbers
  for (multNm in c("vpMult", "qMult", "vp2Mult", "q2Mult")) {
    checkmate::assert_number(control[[multNm]], .var.name = multNm, na.ok = FALSE, finite = TRUE)
  }

  checkmate::assert_data_frame(control$ncaData, min.rows = 1, null.ok = TRUE)
  checkmate::assert_class(control$ncaResults, classes = "PKNCAresults", null.ok = TRUE)
  checkmate::assert_character(control$dvParam, min.chars = 1, len = 1, null.ok = FALSE)
  checkmate::assert_character(control$groups, min.chars = 1, min.len = as.numeric(control$sparse))
  checkmate::assert_logical(control$sparse, len = 1, any.missing = FALSE)
  orig
}

nlmixr2.pkncaEst <- function(object, data, est = NULL,
                             control = list(), table = nlmixr2est::tableControl(),
                             ..., save = NULL, envir = parent.frame()) {
  # Estimate using the ui part of the object
  nlmixr2(
    object = object$ui,
    data = data,
    est = est,
    control = control,
    table = table,
    ...,
    save = save,
    envir = envir
  )
}

#' @export
print.pkncaEst <- function(x, ...) {
  cat("x$ui:\n")
  print(x$ui)
  cat("x$nca:\n")
  print(summary(x$nca))
  invisible(x)
}
