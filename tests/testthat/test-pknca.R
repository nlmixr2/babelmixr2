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
    nlmixr2est::nlmixr(object = modelGood, data = nlmixr2data::theo_sd, est = "pknca"),
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

test_that("est='pknca' with non-mu-referenced models (#101)", {
  nonmumod <- function() {
    ini({
      tka <- 0.45
      tcl <- 0.009
      tv  <- 0.003
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      prop.sd <- 0.7
    })
    model({
      ka <- tka * exp(eta.ka)
      cl <- tcl * exp(eta.cl)
      v <- tv * exp(eta.v)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl / v * center
      cp = center / v
      cp ~ prop(prop.sd)
    })
  }
  mumod <- function() {
    ini({
      tka <- 0.45
      tcl <- -7
      tv  <- -8
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      prop.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl / v * center
      cp = center / v
      cp ~ prop(prop.sd)
    })
  }
  dMod <- nlmixr2data::theo_sd
  dModNoZero <- dMod[(dMod$DV != 0 & dMod$EVID == 0) | (dMod$EVID == 101), ]
  ctl <- pkncaControl(ncaData = dMod, concu = "mg/L", doseu = "mg/kg", timeu = "hr", volumeu = "L/kg")

  suppressMessages(
    fitNonMu <- nlmixr(nonmumod, data = dModNoZero, est = "pknca", control = ctl)
  )
  suppressMessages(
    fitMu <- nlmixr(mumod, data = dModNoZero, est = "pknca", control = ctl)
  )
  expect_s3_class(fitNonMu, "pkncaEst")
  expect_s3_class(fitMu, "pkncaEst")

  feNonMu <- fitNonMu$ui$theta
  feMu <- fitMu$ui$theta
  # Both parameterizations give the same parameter values
  expect_equal(feNonMu[["tka"]], exp(feMu[["tka"]]))
  expect_equal(feNonMu[["tcl"]], exp(feMu[["tcl"]]))
  expect_equal(feNonMu[["tv"]], exp(feMu[["tv"]]))
  # and they were updated from the initial values
  expect_false(isTRUE(all.equal(feNonMu[["tka"]], 0.45)))
  expect_false(isTRUE(all.equal(feNonMu[["tcl"]], 0.009)))
  expect_false(isTRUE(all.equal(feNonMu[["tv"]], 0.003)))
  # bounds are on the parameter scale for the non-mu-referenced model
  iniDf <- fitNonMu$ui$iniDf
  expect_true(all(iniDf$lower[iniDf$name %in% c("tka", "tcl", "tv")] > 0))
})

test_that("pkncaParamMap", {
  model <- function() {
    ini({
      tka <- 0.45
      lcl <- 1
      tvc <- 3
      tq <- 2
      lvp <- 3
      tf <- 0.5
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.q ~ 0.1
      prop.err <- 0.5
    })
    model({
      ka <- (tka) * exp(eta.ka + 0)
      cl <- exp(lcl + eta.cl)
      vc <- tvc
      q <- tq * exp(eta.q)
      vp <- exp(lvp) * WT / 70
      kx <- 3
      kx <- kx + 1
      fdepot <- expit(tf)
      cp <- linCmt()
      cp ~ prop(prop.err)
    })
  }
  suppressMessages(ui <- rxode2::rxode(model))
  paramMap <- pkncaParamMap(ui)
  paramMap <- paramMap[order(paramMap$param), c("theta", "param", "curEval")]
  rownames(paramMap) <- NULL
  expect_equal(paramMap$param, c("cl", "fdepot", "ka", "q", "vc"))
  expect_equal(paramMap$theta, c("lcl", "tf", "tka", "tq", "tvc"))
  expect_equal(paramMap$curEval, c("exp", "expit", "", "", ""))

  # Nothing to change returns the model unchanged
  expect_identical(ini_transform(ui), ui)

  suppressMessages(newmod <- ini_transform(ui, ka = 1.5, cl = 2, fdepot = 0.25, vp = 99))
  expect_equal(newmod$theta[["tka"]], 1.5)
  expect_equal(newmod$theta[["lcl"]], log(2))
  expect_equal(newmod$theta[["tf"]], rxode2::logit(0.25))
  # vp is not a simple function of a single theta, so it is unchanged
  expect_equal(newmod$theta[["lvp"]], 3)
})

test_that("pkncaSimplifyZeroEta", {
  expect_equal(pkncaSimplifyZeroEta(quote(tka * exp(eta.ka)), "eta.ka"), quote(tka))
  expect_equal(pkncaSimplifyZeroEta(quote(exp(eta.ka) * tka), "eta.ka"), quote(tka))
  expect_equal(pkncaSimplifyZeroEta(quote(exp(tka + eta.ka)), "eta.ka"), quote(exp(tka)))
  expect_equal(pkncaSimplifyZeroEta(quote(exp(tka - eta.ka)), "eta.ka"), quote(exp(tka)))
  expect_equal(pkncaSimplifyZeroEta(quote(tka / exp(eta.ka)), "eta.ka"), quote(tka))
  expect_equal(pkncaSimplifyZeroEta(quote(eta.ka + tka), "eta.ka"), quote(tka))
  expect_equal(pkncaSimplifyZeroEta(quote(tka * WT), "eta.ka"), quote(tka * WT))
  expect_equal(pkncaSimplifyZeroEta(3, "eta.ka"), 3)
})
