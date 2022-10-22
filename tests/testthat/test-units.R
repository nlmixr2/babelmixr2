test_that("simplifyUnit", {
  expect_equal(simplifyUnit(""), "")
  expect_equal(simplifyUnit("", ""), "")
  expect_equal(simplifyUnit("kg", ""), "kg")
  expect_equal(simplifyUnit("", "kg"), "1/kg")
  expect_equal(simplifyUnit("kg", NA), "kg")
  expect_equal(simplifyUnit(NA, "kg"), "1/kg")
  expect_equal(simplifyUnit("kg^2", ""), "kg*kg")
  expect_equal(simplifyUnit("", "kg^2"), "1/(kg*kg)")
  expect_equal(simplifyUnit("kg"), "kg")
  expect_equal(simplifyUnit("kg^2"), "kg*kg")
  expect_equal(simplifyUnit("kg", "kg/mL"), "mL")
  expect_equal(simplifyUnit("kg/(kg/mL)"), "mL")

  expect_equal(simplifyUnit(denominator = "mg"), "1/mg")
})

test_that("simplifyUnit expected errors", {
  # units must be recognized by udunits
  expect_error(simplifyUnit("foo"))
  expect_error(simplifyUnit(numerator = c("mg", "mL")))
  expect_error(simplifyUnit(numerator = "mg", denominator = c("mg", "mL")))
})

test_that("modelUnitConversion", {
  # no input gives effectively no output
  expect_equal(
    modelUnitConversion(),
    list(
      amtu = NA_character_,
      clearanceu = NA_character_,
      volumeu = NA_character_,
      timeu = NA_character_,
      dvu = NA_character_,
      cmtu = NA_character_,
      dvConversion = 1
    )
  )
  # Unit conversion is correct
  expect_equal(
    modelUnitConversion(dvu = "ng/mL", amtu = "mg", timeu = "hr", volumeu = "L"),
    list(
      amtu = "mg",
      clearanceu = "L/h",
      volumeu = "L",
      timeu = "hr",
      dvu = "ng/mL",
      cmtu = "mg/L",
      dvConversion = 1000
    )
  )
  # Volume detection works for "L" units
  expect_message(
    conversion <- modelUnitConversion(dvu = "ng/mL", amtu = "mg", timeu = "hr"),
    regexp = "volumeu detected from amtu and dvu as: L"
  )
  expect_equal(
    conversion,
    list(
      amtu = "mg",
      clearanceu = "L/h",
      volumeu = "L",
      timeu = "hr",
      dvu = "ng/mL",
      cmtu = "mg/L",
      dvConversion = 1000
    )
  )
  # Volume detection works for "mL/kg" units
  expect_message(
    conversion <- modelUnitConversion(dvu = "ng/mL", amtu = "mg/kg", timeu = "hr"),
    regexp = "volumeu detected from amtu and dvu as: mL/kg"
  )
  expect_equal(
    conversion,
    list(
      amtu = "mg/kg",
      clearanceu = "mL/(h*kg)",
      volumeu = "mL/kg",
      timeu = "hr",
      dvu = "ng/mL",
      cmtu = "mg/mL",
      dvConversion = 1e6
    )
  )
  # Weird volume detection works for "mL/kg" units
  expect_message(
    conversion <- modelUnitConversion(dvu = "ng/mL", amtu = "mg/g", timeu = "hr"),
    regexp = "volumeu detected from amtu and dvu as: mL/kg"
  )
  expect_equal(
    conversion,
    list(
      amtu = "mg/g",
      clearanceu = "mL/(h*kg)",
      volumeu = "mL/kg",
      timeu = "hr",
      dvu = "ng/mL",
      cmtu = "kg*mg/(g*mL)",
      dvConversion = 1e9
    )
  )
  # Unconvertible volume detection works
  expect_message(
    conversion <- modelUnitConversion(dvu = "ng/mL", amtu = "mol", timeu = "hr"),
    regexp = "volumeu detected from amtu and dvu as: mL*mol/ng",
    fixed = TRUE
  )
  expect_equal(
    conversion,
    list(
      amtu = "mol",
      clearanceu = "mL*mol/(h*ng)",
      volumeu = "mL*mol/ng",
      timeu = "hr",
      dvu = "ng/mL",
      cmtu = "ng/mL",
      dvConversion = 1
    )
  )
  # No time means no clearance
  expect_equal(
    modelUnitConversion(dvu = "ng/mL", amtu = "mg", volumeu = "L"),
    list(
      amtu = "mg",
      clearanceu = NA_character_,
      volumeu = "L",
      timeu = NA_character_,
      dvu = "ng/mL",
      cmtu = "mg/L",
      dvConversion = 1000
    )
  )
})
