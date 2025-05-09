test_that("est='pknca'", {
  modelGood <- function() {
    ini({
      tvka <- 0.45 ; label("Absorption rate (Ka)")
      lcl <- 1 ; label("Clearance (CL)")
      lvc  <- 3.45 ; label("Central volume of distribution (V)")
      prop.err <- 0.5 ; label("Proportional residual error (fraction)")
    })
    model({
      ka <- tvka
      cl <- exp(lcl)
      vc  <- exp(lvc)

      cp <- linCmt()
      cp ~ prop(prop.err)
    })
  }

  # It works with no `control` argument
  suppressMessages(expect_s3_class(
    nlmixr(object = modelGood, data = nlmixr2data::theo_sd, est = "pknca"),
    "pkncaEst"
  ))

  onlyevid0 <- nlmixr2data::theo_sd[nlmixr2data::theo_sd$EVID == 0, ]
  expect_error(
    nlmixr(object = modelGood, data = onlyevid0, est = "pknca"),
    regexp = "no dosing rows (EVID = 1 or 4) detected",
    fixed = TRUE
  )
  noevid0 <- nlmixr2data::theo_sd[nlmixr2data::theo_sd$EVID != 0, ]
  expect_error(
    nlmixr(object = modelGood, data = noevid0, est = "pknca"),
    regexp = "no rows in event table or input data",
    fixed = TRUE
  )
})

test_that("pkncaControl", {
  expect_type(pkncaControl(), "list")
  # All good arguments work
  expect_equal(
    pkncaControl(
      concu = "ng/mL",
      doseu = "mg",
      timeu = "hr",
      volumeu = "L",
      vpMult = 3,
      qMult = 1/3,
      vp2Mult = 6,
      q2Mult = 1/6,
      dvParam = "cp",
      groups = "foo",
      sparse = FALSE
    ),
    list(
      concu = "ng/mL",
      doseu = "mg",
      timeu = "hr",
      volumeu = "L",
      vpMult = 3,
      qMult = 1/3,
      vp2Mult = 6,
      q2Mult = 1/6,
      dvParam = "cp",
      groups = "foo",
      sparse = FALSE,
      ncaData = NULL,
      ncaResults = NULL,
      rxControl= rxode2::rxControl()
    )
  )

  # Confirm some degree of error checking on all arguments
  expect_error(pkncaControl(concu = 1))
  expect_error(pkncaControl(doseu = 1))
  expect_error(pkncaControl(timeu = 1))
  expect_error(pkncaControl(volumeu = 1))
  expect_error(pkncaControl(vpMult = "A"))
  expect_error(pkncaControl(qMult = "A"))
  expect_error(pkncaControl(vp2Mult = "A"))
  expect_error(pkncaControl(q2Mult = "A"))
  expect_error(pkncaControl(dvParam = 1))
  expect_error(pkncaControl(groups = 1))
  expect_error(pkncaControl(sparse = NA))
  expect_error(pkncaControl(ncaData = 1))
  expect_error(pkncaControl(ncaResults = 1))
})

test_that("ini_transform", {
  model <- function() {
    ini({
      tvka <- 0.45 ; label("Absorption rate (Ka)")
      lcl <- 1 ; label("Clearance (CL)")
      lvc  <- 3.45 ; label("Central volume of distribution (V)")
      prop.err <- 0.5 ; label("Proportional residual error (fraction)")
    })
    model({
      ka <- tvka
      cl <- exp(lcl)
      vc  <- exp(lvc)

      linCmt() ~ prop(prop.err)
    })
  }
  suppressMessages(newmod <- ini_transform(rxode2::rxode(model), ka=1.5, cl=2, lvc=3))
  expect_equal(fixef(newmod)[["tvka"]], 1.5)
  expect_equal(fixef(newmod)[["lcl"]], log(2))
  expect_equal(fixef(newmod)[["lvc"]], 3)
})

test_that("dvParam", {
  modelBad <- function() {
    ini({
      tvka <- 0.45 ; label("Absorption rate (Ka)")
      lcl <- 1 ; label("Clearance (CL)")
      lvc  <- 3.45 ; label("Central volume of distribution (V)")
      prop.err <- 0.5 ; label("Proportional residual error (fraction)")
    })
    model({
      ka <- tvka
      cl <- exp(lcl)
      vc  <- exp(lvc)

      linCmt() ~ prop(prop.err)
    })
  }
  modelGood <- function() {
    ini({
      tvka <- 0.45 ; label("Absorption rate (Ka)")
      lcl <- 1 ; label("Clearance (CL)")
      lvc  <- 3.45 ; label("Central volume of distribution (V)")
      prop.err <- 0.5 ; label("Proportional residual error (fraction)")
    })
    model({
      ka <- tvka
      cl <- exp(lcl)
      vc  <- exp(lvc)

      cp <- linCmt()
      cp ~ prop(prop.err)
    })
  }

  suppressMessages(expect_s3_class(
    nlmixr(object = modelGood, data = nlmixr2data::theo_sd, est = "pknca"),
    "pkncaEst"
  ))
  # modelBad is okay if unit conversion is not required
  suppressMessages(expect_s3_class(
    nlmixr(object = modelBad, data = nlmixr2data::theo_sd, est = "pknca"),
    "pkncaEst"
  ))
  # modelBad is not okay if unit conversion is required
  skip_if_not_installed("PKNCA", "0.10.0.9000") # this test will fail due to https://github.com/humanpred/pknca/pull/191
  suppressMessages(expect_error(
    nlmixr(
      object = modelBad, data = nlmixr2data::theo_sd,
      est = "pknca",
      control = pkncaControl(concu = "ng/mL", doseu = "mg", timeu = "hr", volumeu = "L")
    ),
    regexp = "Could not detect DV assignment for unit conversion"
  ))
})

test_that("getDvLines", {
  modelBad <- function() {
    ini({
      tvka <- 0.45 ; label("Absorption rate (Ka)")
      lcl <- 1 ; label("Clearance (CL)")
      lvc  <- 3.45 ; label("Central volume of distribution (V)")
      prop.err <- 0.5 ; label("Proportional residual error (fraction)")
    })
    model({
      ka <- tvka
      cl <- exp(lcl)
      vc  <- exp(lvc)

      linCmt() ~ prop(prop.err)
    })
  }
  modelGood <- function() {
    ini({
      tvka <- 0.45 ; label("Absorption rate (Ka)")
      lcl <- 1 ; label("Clearance (CL)")
      lvc  <- 3.45 ; label("Central volume of distribution (V)")
      prop.err <- 0.5 ; label("Proportional residual error (fraction)")
    })
    model({
      ka <- tvka
      cl <- exp(lcl)
      vc  <- exp(lvc)

      cp <- linCmt()
      cp ~ prop(prop.err)
    })
  }

  expect_equal(
    getDvLines(modelBad),
    list(str2lang("linCmt() ~ prop(prop.err)"))
  )
  expect_equal(
    getDvLines(modelGood),
    list(str2lang("cp ~ prop(prop.err)"))
  )
  expect_equal(
    getDvLines(modelGood, dvAssign = "cp"),
    list(str2lang("cp <- linCmt()"))
  )
})
