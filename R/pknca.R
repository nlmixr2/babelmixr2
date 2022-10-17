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
  browser()
  control <- env$control[[1]]
  # Get the units from the basic units (before unit conversion)
  dUnitsData <-
    PKNCA::pknca_units_table(
      concu = control$concu,
      doseu = control$doseu,
      timeu = control$timeu
    )

  oNCA <- calcPknca(env, pkncaUnits = dUnitsData)

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
  # Drop unit conversion when units are not specified
  modelDataConversions <-
    modelDataConversions[
      !is.na(modelDataConversions$PPORRESU) &
        !is.na(modelDataConversions$PPSTRESU),
    ]
  dUnitsModel <-
    PKNCA::pknca_units_table(
      concu = control$concu,
      doseu = control$doseu,
      timeu = control$timeu,
      conversions = modelDataConversions
    )
  # set any missing conversion factors to 1
  if (!("conversion_factor" %in% names(dUnitsModel))) {
    dUnitsModel$conversion_factor <- 1
  } else {
    dUnitsModel$conversion_factor[is.na(dUnitsModel$conversion_factor)] <- 1
  }
  unitConversions <-
    stats::setNames(
      dUnitsModel$conversion_factor,
      nm = dUnitsModel$PPTESTCD
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

  # Add unit conversion to the estimation
  # TODO: This will only work when cp is assigned to center/vc.  Need detection
  # of the correct unit conversion assignment.
  browser()
  newEnv <-
    eval(str2lang(
      sprintf("rxode2::model(newEnv, cp <- %g*center/vc)", 1/unitConversions[["cmax"]])
    ))

  env$ui <- newEnv
  env$nca <- oNCA

  env
}

#' Perform NCA calculations with PKNCA
#'
#' @inheritParams nlmixr2est::nlmixr2Est
#' @return A PKNCAresults object
#' @noRd
calcPknca <- function(env, pkncaUnits) {
  # Normalize column names
  cleanData <- bblDatToPknca(model = env$ui, data = env$data)
  control <- env$control[[1]]
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
        sort(log(2)/(tmaxValues/4)),
      vc=
        sort(unitConversions[["vss.last"]] / cmaxdnValues),
      cl=
        unitConversions[["cl.last"]] * cllastValues
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
  checkmate::expect_names(names(changeArgs))
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
#' @return A list of parameters
#' @export
pkncaControl <- function(concu = NA_character_, doseu = NA_character_, timeu = NA_character_,
                         volumeu = NA_character_,
                         vpMult=2, qMult=1/2,
                         vp2Mult=4, q2Mult=1/4,
                         dvParam = "cp",
                         groups = character(),
                         sparse = FALSE) {
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
      sparse = sparse
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
    }
    control <- orig[[1]]
  }
  checkmate::expect_names(
    x = names(control),
    permutation.of = names(formals(pkncaControl))
  )
  # verify units look like units
  for (unitNm in c("concu", "doseu", "timeu", "volumeu")) {
    checkmate::expect_character(control[[unitNm]], label = unitNm, null.ok = FALSE, len = 1, min.chars = 1)
  }
  # Verify that multipliers are numbers
  for (multNm in c("vpMult", "qMult", "vp2Mult", "q2Mult")) {
    checkmate::expect_number(control[[multNm]], label = multNm, na.ok = FALSE, finite = TRUE)
  }

  checkmate::expect_character(control$dvParam, min.chars = 1, len = 1, null.ok = FALSE)
  checkmate::expect_character(control$groups, min.chars = 1, min.len = as.numeric(control$sparse))
  checkmate::expect_logical(control$sparse, len = 1, any.missing = FALSE)
  orig
}
