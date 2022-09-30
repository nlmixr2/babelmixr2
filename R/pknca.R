#' Estimate starting parameters using PKNCA
#'
#' @details Parameters are estimated as follows:
#'
#' \itemize{
#' \item{ka}{4 half-lives to Tmax but not higher than 3:  \code{log(2)/(tmax/4)}}
#' \item{vc}{Inverse of dose-normalized Cmax}
#' \item{cl}{Estimated as the median clearance}
#' \item{vp,vp2}{2- and 4-fold the \code{vc}, respectively}
#' \item{q,q2}{0.5- and 0.25-fold the \code{cl}, respectively}
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
  # TODO: Detect column names for EVID, CMT, ID
  obsData <- env$data[env$data$EVID == 0, ]
  doseData <- env$data[env$data$EVID %in% c(1, 4), ]
  obsCmt <- unique(obsData$CMT)
  doseCmt <- unique(doseData$CMT)
  stopifnot(length(obsCmt) == 1)
  stopifnot(length(doseCmt) == 1)
  control <- env$control[[1]]
  oConcFormula <-
    stats::as.formula(paste0(
      "DV~TIME|", paste(c(control$groups, "ID"), collapse="+")
    ))
  oDoseFormula <-
    stats::as.formula(paste0(
      "AMT~TIME|", paste(c(control$groups, "ID"), collapse="+")
    ))
  oConc <- PKNCA::PKNCAconc(data = obsData, oConcFormula, sparse = control$sparse)
  doseRoute <- ifelse(obsCmt == doseCmt, yes = "intravascular", no = "extravascular")
  oDose <- PKNCA::PKNCAdose(data = doseData, oDoseFormula, route = doseRoute)

  # Get the units from the basic units (before unit conversion)
  dUnitsData <-
    PKNCA::pknca_units_table(
      concu = control$concu,
      doseu = control$doseu,
      timeu = control$timeu
    )
  # Setup unit conversions for human-usable units
  unitsSetup <- setNames(dUnitsData$PPORRESU, nm = dUnitsData$PPTESTCD)
  if (is.null(control$volumeu)) {
    # simplify the volume units
    volumeuRaw <-
      attr(
        units::set_units(1, unitsSetup[["vss.last"]], mode = "standard"),
        "units"
      )
    for (idx in rev(seq_along(volumeuRaw$numerator))) {
      foundDenom <- which(volumeuRaw$denominator == volumeuRaw$numerator[idx])
      if (length(foundDenom) > 0) {
        # Cancel out the units
        volumeuRaw$numerator <- volumeuRaw$numerator[-idx]
        volumeuRaw$denominator <- volumeuRaw$denominator[-foundDenom[1]]
      }
    }
    numerator <- paste(volumeuRaw$numerator, collapse = "*")
    if (length(volumeuRaw$numerator) > 1) {
      numerator <- sprintf("(%s)", numerator)
    }
    denominator <- paste(volumeuRaw$denominator, collapse = "*")
    if (length(volumeuRaw$denominator) > 1) {
      denominator <- sprintf("(%s)", denominator)
    }
    control$volumeu <- paste0(numerator, "/", denominator)
  }

  clearanceu <- paste0(control$volumeu, "/", control$timeu)
  modelConcu <- paste0("(", control$doseu, ")/(", control$volumeu, ")")
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
          modelConcu,
          paste0("1/(", control$volumeu, ")"),
          clearanceu,
          control$volumeu
        )
    )
  dUnitsModel <-
    PKNCA::pknca_units_table(
      concu = control$concu,
      doseu = control$doseu,
      timeu = control$timeu,
      conversions = modelDataConversions
    )
  unitConversions <-
    stats::setNames(
      dUnitsModel$conversion_factor,
      nm = dUnitsModel$PPTESTCD
    )

  oData <- PKNCA::PKNCAdata(oConc, oDose, units = dUnitsData)
  intervals <- oData$intervals
  intervals$cl.last <- intervals$auclast
  intervals$cmax.dn <- intervals$cmax
  intervals$vss.last <- intervals$auclast
  oData$intervals <- intervals
  oNCA <- PKNCA::pk.nca(oData)

  ncaParams <- as.data.frame(oNCA)
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

  # What parameters should be modified?  And then modify them.
  murefNames <- env$ui$getSplitMuModel$pureMuRef
  updateNames <- intersect(murefNames, names(ncaEstimates))
  newEnv <- do.call(ini_transform, append(list(x=env$ui), ncaEstimates[updateNames]))

  # Add unit conversion to the estimation
  # TODO: This will only work when cp is assigned to center/vc.  Need detection
  # of the correct unit conversion assignment.
  newEnv <-
    eval(str2lang(
      sprintf("rxode2::model(newEnv, cp <- %g*center/vc)", 1/unitConversions[["cmax"]])
    ))

  env$ui <- newEnv
  env$nca <- oNCA

  env
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
      transFun <- inverseTrans[[currentTrans]]
      if (is.null(transFun)) {
        cli::cli_abort(paste("cannot invert the transform:", transFun))
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
#' @param concu,doseu,amountu,timeu concentration, dose, dosing amount, and time
#'   units from the source data (passed to \code{PKNCA::pknca_units_table()}).
#' @param volumeu compartment volume for the model (if \code{NULL}, simplified
#'   units from source data will be used)
#' @param vpMult,qMult,vp2Mult,q2Mult Multipliers for vc and cl to provide
#'   initial estimates for vp, q, vp2, and q2
#' @param groups Grouping columns for NCA summaries by group (required if
#'   \code{sparse = TRUE})
#' @return A list of parameters
#' @export
pkncaControl <- function(concu = NULL, doseu = NULL, amountu = NULL, timeu = NULL,
                         volumeu = NULL,
                         vpMult=2, qMult=1/2,
                         vp2Mult=4, q2Mult=1/4,
                         groups = character(),
                         sparse = FALSE) {
  getValidNlmixrCtl.pknca(
    list(
      concu = concu,
      doseu = doseu,
      timeu = timeu,
      amountu = amountu,
      volumeu = volumeu,
      groups = groups,
      sparse = sparse
    )
  )
}

#' @export
getValidNlmixrCtl.pknca <- function(control) {
  orig <- control
  if (inherits(control, "getValidNlmixrControl")) {
    control <- orig[[1]]
  }
  checkmate::expect_names(
    x = names(control),
    permutation.of = names(formals(pkncaControl))
  )
  # verify units look like units
  for (unitNm in c("concu", "doseu", "amountu", "timeu",
                   "volumeu")) {
    checkmate::expect_character(control[[unitNm]], label = unitNm, null.ok = TRUE, len = 1, min.chars = 1)
  }
  # Verify that multipliers are numbers
  for (multNm in c("vpMult", "qMult", "vp2Mult", "q2Mult")) {
    checkmate::expect_number(control[[multNm]], label = multNm, na.ok = FALSE, finite = TRUE)
  }

  checkmate::expect_logical(control$sparse, len = 1, any.missing = FALSE)
  checkmate::expect_character(control$groups, min.chars = 1, min.len = as.numeric(control$sparse))
  orig
}
