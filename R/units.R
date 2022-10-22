#' Simplify units by removing repeated units from the numerator and denominator
#'
#' @details \code{NA} or \code{""} for \code{numerator} and \code{denominator}
#'   are considered unitless.
#'
#' @param numerator The numerator of the units (or the whole unit specification)
#' @param denominator The denominator of the units (or NULL if \code{numerator}
#'   is the whole unit specification)
#' @return The units specified with units that are in both the numerator and
#'   denominator cancelled.
#' @family Unit conversion
#' @examples
#' simplifyUnit("kg", "kg/mL")
#' # units that don't match exactly are not cancelled
#' simplifyUnit("kg", "g/mL")
#' @export
simplifyUnit <- function(numerator="", denominator="") {
  checkmate::expect_character(numerator, len = 1)
  checkmate::expect_character(denominator, len = 1)
  if (is.na(numerator)) numerator <- ""
  if (is.na(denominator)) denominator <- ""
  hasNumeratorInput <- nchar(numerator) > 0
  hasDenominatorInput <- nchar(denominator) > 0
  if (!hasNumeratorInput & !hasDenominatorInput) {
    x <- ""
  } else if (!hasNumeratorInput) {
    x <- sprintf("1/(%s)", denominator)
  } else if (!hasDenominatorInput) {
    x <- numerator
  } else {
    x <- sprintf("(%s)/(%s)", numerator, denominator)
  }

  unitToSimplify <- attr(units::as_units(x), "unit")
  for (idx in rev(seq_along(unitToSimplify$numerator))) {
    foundDenom <- which(unitToSimplify$denominator == unitToSimplify$numerator[idx])
    if (length(foundDenom) > 0) {
      # Cancel out the units
      unitToSimplify$numerator <- unitToSimplify$numerator[-idx]
      unitToSimplify$denominator <- unitToSimplify$denominator[-foundDenom[1]]
    }
  }
  hasNumerator <- length(unitToSimplify$numerator) > 0
  hasDenominator <- length(unitToSimplify$denominator) > 0
  hasLongDenominator <- length(unitToSimplify$denominator) > 1
  numeratorText <- paste(unitToSimplify$numerator, collapse = "*")
  denominatorText <-
    if (hasLongDenominator) {
      sprintf("(%s)", paste(unitToSimplify$denominator, collapse="*"))
    } else {
      unitToSimplify$denominator
    }
  if (!hasNumerator & !hasDenominator) {
    # everything cancelled
    ret <- ""
  } else if (hasNumerator & !hasDenominator) {
    # has numerator with no denominator
    ret <- numeratorText
  } else if (!hasNumerator & hasDenominator) {
    # has denominator with no numerator
    ret <- sprintf("1/%s", denominatorText)
  } else {
    # has both numerator and denominator
    ret <- sprintf("%s/%s", numeratorText, denominatorText)
  }
  ret
}

#' Unit conversion for pharmacokinetic models
#'
#' @param dvu,amtu,timeu The units for the DV, AMT, and TIME columns in the data
#' @param volumeu The units for the volume parameters in the model
#' @return A list with names for the units associated with each parameter
#'   ("amtu", "clearanceu", "volumeu", "timeu", "dvu") and the numeric value to
#'   multiply the modeled estimate (for example, \code{cp}) so that the model is
#'   consistent with the data units.
#' @family Unit conversion
#' @examples
#' modelUnitConversion(dvu = "ng/mL", amtu = "mg", timeu = "hr", volumeu = "L")
#' @export
modelUnitConversion <- function(dvu = NA_character_, amtu = NA_character_, timeu = NA_character_, volumeu = NA_character_) {
  checkmate::assert_character(dvu, min.chars = 1, len = 1, null.ok = FALSE, any.missing = TRUE)
  checkmate::assert_character(amtu, min.chars = 1, len = 1, null.ok = FALSE, any.missing = TRUE)
  checkmate::assert_character(timeu, min.chars = 1, len = 1, null.ok = FALSE, any.missing = TRUE)
  checkmate::assert_character(volumeu, min.chars = 1, len = 1, null.ok = FALSE, any.missing = TRUE)

  dvuConversion <- 1
  dvuBase <- NA_character_
  if (!is.na(dvu) & !is.na(amtu)) {
    # calculate volumeu and dvuConversion
    volumeuBase <- simplifyUnit(numerator=amtu, denominator=dvu)
    if (is.na(volumeu)) {
      # auto-generate the volumeu
      if (units::ud_are_convertible(x = volumeuBase, y = "L")) {
        volumeu <- "L"
      } else if (units::ud_are_convertible(x = volumeuBase, y = "mL/kg")) {
        volumeu <- "mL/kg"
      } else {
        volumeu <- volumeuBase
      }
      cli::cli_alert_info(paste("volumeu detected from amtu and dvu as:", volumeu))
    }
    # default units for dvu with no conversion
    dvuBase <- simplifyUnit(numerator = amtu, denominator = volumeu)
    if (units::ud_are_convertible(x = dvuBase, y = dvu)) {
      dvuConversion <-
        units::set_units(
          units::as_units(x = dvuBase),
          value = dvu,
          mode = "standard"
        )
    } else {
      dvuOrig <- dvu
      dvu <- simplifyUnit(numerator = dvu, denominator = dvuBase)
      cli::cli_warn(sprintf(
        "dvu cannot be converted from %s to %s, actual dvu for the model may be hard to interpret (%s)",
        dvuBase, dvuOrig, dvu
      ))
    }
  }

  clearanceu <- NA_character_
  if (!is.na(volumeu) & !is.na(timeu)) {
    clearanceu <- simplifyUnit(numerator = volumeu, denominator = timeu)
  }

  list(
    amtu = amtu,
    clearanceu = clearanceu,
    volumeu = volumeu,
    timeu = timeu,
    dvu = dvu,
    cmtu = dvuBase,
    dvConversion = as.numeric(dvuConversion)
  )
}
