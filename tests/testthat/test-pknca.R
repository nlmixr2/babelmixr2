test_that("est='pknca'", {
  modelGood <- function() {
    ini({
      tvka <- 0.45
      label("Absorption rate (Ka)")
      lcl <- 1
      label("Clearance (CL)")
      lvc <- 3.45
      label("Central volume of distribution (V)")
      prop.err <- 0.5
      label("Proportional residual error (fraction)")
    })
    model({
      ka <- tvka
      cl <- exp(lcl)
      vc <- exp(lvc)

      cp <- linCmt()
      cp ~ prop(prop.err)
    })
  }

  # It works with no `control` argument
  suppressMessages(expect_s3_class(
    nlmixr2est::nlmixr(
      object = modelGood,
      data = nlmixr2data::theo_sd,
      est = "pknca"
    ),
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
      qMult = 1 / 3,
      vp2Mult = 6,
      q2Mult = 1 / 6,
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
      qMult = 1 / 3,
      vp2Mult = 6,
      q2Mult = 1 / 6,
      dvParam = "cp",
      groups = "foo",
      sparse = FALSE,
      ncaData = NULL,
      ncaResults = NULL,
      rxControl = rxode2::rxControl()
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
      tvka <- 0.45
      label("Absorption rate (Ka)")
      lcl <- 1
      label("Clearance (CL)")
      lvc <- 3.45
      label("Central volume of distribution (V)")
      prop.err <- 0.5
      label("Proportional residual error (fraction)")
    })
    model({
      ka <- tvka
      cl <- exp(lcl)
      vc <- exp(lvc)

      linCmt() ~ prop(prop.err)
    })
  }
  suppressMessages(
    newmod <- ini_transform(rxode2::rxode(model), ka = 1.5, cl = 2, lvc = 3)
  )
  expect_equal(fixef(newmod)[["tvka"]], 1.5)
  expect_equal(fixef(newmod)[["lcl"]], log(2))
  expect_equal(fixef(newmod)[["lvc"]], 3)
})

test_that("dvParam", {
  modelBad <- function() {
    ini({
      tvka <- 0.45
      label("Absorption rate (Ka)")
      lcl <- 1
      label("Clearance (CL)")
      lvc <- 3.45
      label("Central volume of distribution (V)")
      prop.err <- 0.5
      label("Proportional residual error (fraction)")
    })
    model({
      ka <- tvka
      cl <- exp(lcl)
      vc <- exp(lvc)

      linCmt() ~ prop(prop.err)
    })
  }
  modelGood <- function() {
    ini({
      tvka <- 0.45
      label("Absorption rate (Ka)")
      lcl <- 1
      label("Clearance (CL)")
      lvc <- 3.45
      label("Central volume of distribution (V)")
      prop.err <- 0.5
      label("Proportional residual error (fraction)")
    })
    model({
      ka <- tvka
      cl <- exp(lcl)
      vc <- exp(lvc)

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
      object = modelBad,
      data = nlmixr2data::theo_sd,
      est = "pknca",
      control = pkncaControl(
        concu = "ng/mL",
        doseu = "mg",
        timeu = "hr",
        volumeu = "L"
      )
    ),
    regexp = "Could not detect DV assignment for unit conversion"
  ))
})

test_that("getDvLines", {
  modelBad <- function() {
    ini({
      tvka <- 0.45
      label("Absorption rate (Ka)")
      lcl <- 1
      label("Clearance (CL)")
      lvc <- 3.45
      label("Central volume of distribution (V)")
      prop.err <- 0.5
      label("Proportional residual error (fraction)")
    })
    model({
      ka <- tvka
      cl <- exp(lcl)
      vc <- exp(lvc)

      linCmt() ~ prop(prop.err)
    })
  }
  modelGood <- function() {
    ini({
      tvka <- 0.45
      label("Absorption rate (Ka)")
      lcl <- 1
      label("Clearance (CL)")
      lvc <- 3.45
      label("Central volume of distribution (V)")
      prop.err <- 0.5
      label("Proportional residual error (fraction)")
    })
    model({
      ka <- tvka
      cl <- exp(lcl)
      vc <- exp(lvc)

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
      tv <- 0.003
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      prop.sd <- 0.7
    })
    model({
      ka <- tka * exp(eta.ka)
      cl <- tcl * exp(eta.cl)
      v <- tv * exp(eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ prop(prop.sd)
    })
  }
  mumod <- function() {
    ini({
      tka <- 0.45
      tcl <- -7
      tv <- -8
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      prop.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ prop(prop.sd)
    })
  }
  dMod <- nlmixr2data::theo_sd
  dModNoZero <- dMod[(dMod$DV != 0 & dMod$EVID == 0) | (dMod$EVID == 101), ]
  ctl <- pkncaControl(
    ncaData = dMod,
    concu = "mg/L",
    doseu = "mg/kg",
    timeu = "hr",
    volumeu = "L/kg"
  )

  suppressMessages(
    fitNonMu <- nlmixr(
      nonmumod,
      data = dModNoZero,
      est = "pknca",
      control = ctl
    )
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

  # When vc is in the model, v is not given the vc estimate
  vcmod <- function() {
    ini({
      tka <- 0.45
      tcl <- 0.009
      tvc <- 0.004
      tv <- 0.003
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      prop.sd <- 0.7
    })
    model({
      ka <- tka * exp(eta.ka)
      cl <- tcl * exp(eta.cl)
      vc <- tvc * exp(eta.v)
      v <- tv * exp(eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / vc * center
      cp <- center / vc + 0 * v
      cp ~ prop(prop.sd)
    })
  }
  suppressMessages(
    fitVc <- nlmixr(vcmod, data = dModNoZero, est = "pknca", control = ctl)
  )
  expect_equal(fitVc$ui$theta[["tv"]], 0.003)
  expect_equal(fitVc$ui$theta[["tvc"]], feNonMu[["tv"]])

  # ... even when vc cannot be updated (not a simple function of one theta)
  suppressMessages(
    vcCovMod <- rxode2::model(vcmod, vc <- tvc * 2 * exp(eta.v))
  )
  suppressMessages(
    fitVcCov <- nlmixr(
      vcCovMod,
      data = dModNoZero,
      est = "pknca",
      control = ctl
    )
  )
  expect_equal(fitVcCov$ui$theta[["tv"]], 0.003)
  expect_equal(fitVcCov$ui$theta[["tvc"]], 0.004)
  expect_message(
    nlmixr(vcCovMod, data = dModNoZero, est = "pknca", control = ctl),
    regexp = "NCA initial estimates not applied to `vc`"
  )

  # ... or when vc is a theta used directly
  vcThetaMod <- function() {
    ini({
      tka <- 0.45
      tcl <- 0.009
      vc <- 0.004
      tv <- 0.003
      eta.v ~ 0.1
      prop.sd <- 0.7
    })
    model({
      ka <- tka
      cl <- tcl
      v <- tv * exp(eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / vc * center
      cp <- center / vc + 0 * v
      cp ~ prop(prop.sd)
    })
  }
  suppressMessages(
    fitVcTheta <- nlmixr(
      vcThetaMod,
      data = dModNoZero,
      est = "pknca",
      control = ctl
    )
  )
  expect_equal(fitVcTheta$ui$theta[["tv"]], 0.003)
  expect_equal(fitVcTheta$ui$theta[["vc"]], feNonMu[["tv"]])

  # A peripheral v is not given the central volume when the central volume
  # is named V
  vPeriphMod <- function() {
    ini({
      tka <- 0.45
      tcl <- 0.009
      tV <- 0.004
      tv <- 0.003
      tq <- 0.001
      eta.v ~ 0.1
      prop.sd <- 0.7
    })
    model({
      ka <- tka
      cl <- tcl
      V <- tV * exp(eta.v)
      v <- tv
      q <- tq
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka *
        depot -
        cl / V * center -
        q / V * center +
        q / v * periph
      d / dt(periph) <- q / V * center - q / v * periph
      cp <- center / V
      cp ~ prop(prop.sd)
    })
  }
  # With more than one central volume name, it is ambiguous
  expect_message(
    fitVPeriph <- nlmixr(
      vPeriphMod,
      data = dModNoZero,
      est = "pknca",
      control = ctl
    ),
    regexp = "the central volume could be any of"
  )
  expect_equal(fitVPeriph$ui$theta[["tv"]], 0.003)
  expect_equal(fitVPeriph$ui$theta[["tV"]], 0.004)

  # Another central volume name gets the central volume estimate
  suppressMessages(vcMod <- rxode2::rxRename(nonmumod, Vc = v, tVc = tv))
  suppressMessages(
    fitVc2 <- nlmixr(vcMod, data = dModNoZero, est = "pknca", control = ctl)
  )
  expect_equal(fitVc2$ui$theta[["tVc"]], feNonMu[["tv"]])

  # A v that cannot be updated (and no vc) is reported
  suppressMessages(
    vComplexMod <- rxode2::model(nonmumod, v <- tv * 2 * exp(eta.v))
  )
  expect_message(
    nlmixr(vComplexMod, data = dModNoZero, est = "pknca", control = ctl),
    regexp = "NCA initial estimates not applied to `v`"
  )
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
      if (WT > 70) {
        q <- tq * 2 * exp(eta.q)
      }
      fdepot <- expit(tf)
      cp <- linCmt()
      cp ~ prop(prop.err)
    })
  }
  suppressMessages(ui <- rxode2::rxode(model))
  paramMap <- pkncaParamMap(ui)
  paramMap <- paramMap[
    paramMap$param %in% c("cl", "fdepot", "ka", "q", "vc", "vp"),
  ]
  paramMap <- paramMap[order(paramMap$param), c("theta", "param", "curEval")]
  rownames(paramMap) <- NULL
  # q is also assigned in an if block, so it is ambiguous and not mapped; vp
  # is not a simple function of a single theta
  expect_equal(paramMap$param, c("cl", "fdepot", "ka", "vc"))
  expect_equal(paramMap$theta, c("lcl", "tf", "tka", "tvc"))
  expect_equal(paramMap$curEval, c("exp", "expit", "", ""))

  # Nothing to change returns the model unchanged
  expect_identical(ini_transform(ui), ui)

  suppressMessages(
    newmod <- ini_transform(ui, ka = 1.5, cl = 2, fdepot = 0.25, vp = 99)
  )
  expect_equal(newmod$theta[["tka"]], 1.5)
  expect_equal(newmod$theta[["lcl"]], log(2))
  expect_equal(newmod$theta[["tf"]], rxode2::logit(0.25))
  # vp is not a simple function of a single theta, so it is unchanged
  expect_equal(newmod$theta[["lvp"]], 3)
})

test_that("ini_transform with the same theta and parameter name", {
  model <- function() {
    ini({
      cl <- 1
      tv <- 3
      tf <- 0
      eta.cl ~ 0.1
      eta.f ~ 0.1
      prop.err <- 0.5
    })
    model({
      cl <- exp(cl + eta.cl)
      v <- tv * exp(-eta.cl)
      fdepot <- expit(tf) * exp(eta.f)
      cp <- linCmt()
      cp ~ prop(prop.err)
    })
  }
  suppressMessages(ui <- rxode2::rxode(model))
  suppressMessages(newmod <- ini_transform(ui, cl = 2, v = 5, fdepot = 0.25))
  expect_equal(newmod$theta[["cl"]], log(2))
  expect_equal(newmod$theta[["tv"]], 5)
  expect_equal(newmod$theta[["tf"]], rxode2::logit(0.25))

  # expit() bounds are respected
  model2 <- rxode2::model(ui, fdepot <- expit(tf, -1, 2) * exp(eta.f))
  suppressMessages(newmod <- ini_transform(model2, fdepot = 0.25))
  expect_equal(newmod$theta[["tf"]], rxode2::logit(0.25, -1, 2))

  # An estimate outside of the expit() bounds leaves the theta unchanged
  expect_warning(
    suppressMessages(newmod <- ini_transform(model2, fdepot = 3)),
    regexp = "cannot transform the estimate for `fdepot`"
  )
  expect_equal(newmod$theta[["tf"]], 0)
})

test_that("pkncaParamMap skips thetas shared by more than one parameter", {
  model <- function() {
    ini({
      tpop <- 1
      tvc <- 3
      tv <- 4
      eta.ka ~ 0.1
      prop.err <- 0.5
    })
    model({
      ka <- tpop * exp(eta.ka)
      cl <- exp(tpop)
      vc <- tvc * exp(eta.ka)
      v <- tv * exp(eta.ka)
      cp <- linCmt()
      cp ~ prop(prop.err)
    })
  }
  suppressMessages(ui <- rxode2::rxode(model))
  paramMap <- pkncaParamMap(ui)
  expect_false("tpop" %in% paramMap$theta)
  expect_true(all(c("vc", "v") %in% paramMap$param))
})

test_that("est='pknca' with parameters defined in ini()", {
  model <- function() {
    ini({
      ka <- 0.45
      cl <- 1
      vc <- 3.45
      prop.err <- 0.5
    })
    model({
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / vc * center
      cp <- center / vc
      cp ~ prop(prop.err)
    })
  }
  suppressMessages(ui <- rxode2::rxode(model))
  paramMap <- pkncaParamMap(ui)
  expect_true(all(c("ka", "cl", "vc") %in% paramMap$param))
  suppressMessages(
    fit <- nlmixr(model, data = nlmixr2data::theo_sd, est = "pknca")
  )
  expect_false(isTRUE(all.equal(fit$ui$theta[["ka"]], 0.45)))
  expect_false(isTRUE(all.equal(fit$ui$theta[["cl"]], 1)))
  expect_false(isTRUE(all.equal(fit$ui$theta[["vc"]], 3.45)))
  expect_equal(fit$ui$theta[["prop.err"]], 0.5)
})

test_that("pkncaParamMap skips direct thetas used in transformations", {
  model <- function() {
    ini({
      tka <- 0.45
      cl <- log(4)
      tv <- 3
      eta.cl ~ 0.3
      prop.err <- 0.5
    })
    model({
      ka <- tka
      CL <- exp(cl + 0.75 * log(WT / 70) + eta.cl)
      v <- tv
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - CL / v * center
      cp <- center / v
      cp ~ prop(prop.err)
    })
  }
  suppressMessages(ui <- rxode2::rxode(model))
  paramMap <- pkncaParamMap(ui)
  expect_false("cl" %in% paramMap$param)
  suppressMessages(newmod <- ini_transform(ui, cl = c(0.3, 3, 30)))
  expect_equal(newmod$theta[["cl"]], log(4))
})

test_that("ini_transform allows infinite bounds", {
  model <- function() {
    ini({
      tka <- 0
      eta.ka ~ 0.1
      tcl <- 1
      tvc <- 3
      prop.err <- 0.5
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl)
      vc <- exp(tvc)
      cp <- linCmt()
      cp ~ prop(prop.err)
    })
  }
  suppressMessages(ui <- rxode2::rxode(model))
  suppressMessages(newmod <- ini_transform(ui, ka = c(0.01, 1, Inf)))
  expect_equal(newmod$theta[["tka"]], 0)
  iniDf <- newmod$iniDf
  expect_equal(iniDf$lower[iniDf$name == "tka"], log(0.01))
  expect_equal(iniDf$upper[iniDf$name == "tka"], Inf)
})

test_that("pkncaParamMap: covariates and linCmt() arguments", {
  covModel <- function() {
    ini({
      ka <- 0.45
      cl <- 1
      vc <- 3
      prop.err <- 0.5
    })
    model({
      CL <- cl * WT
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - CL / vc * center
      cp <- center / vc
      cp ~ prop(prop.err)
    })
  }
  suppressMessages(ui <- rxode2::rxode(covModel))
  paramMap <- pkncaParamMap(ui)
  # vc is divided into the derived CL, so it is not used plainly either
  expect_equal(intersect(c("ka", "cl", "vc"), paramMap$param), "ka")

  linModel <- function() {
    ini({
      ka <- 0.45
      cl <- 1
      vc <- 3.45
      prop.err <- 0.5
    })
    model({
      cp <- linCmt(ka, cl, vc)
      cp ~ prop(prop.err)
    })
  }
  suppressMessages(ui <- rxode2::rxode(linModel))
  expect_true(all(c("ka", "cl", "vc") %in% pkncaParamMap(ui)$param))
  suppressMessages(
    fit <- nlmixr(linModel, data = nlmixr2data::theo_sd, est = "pknca")
  )
  expect_false(isTRUE(all.equal(fit$ui$theta[["cl"]], 1)))
  expect_false(isTRUE(all.equal(fit$ui$theta[["vc"]], 3.45)))

  # Scaled by a constant or another parameter: left unchanged with a message
  scaledModel <- function() {
    ini({
      ka <- 0.45
      cl <- 1
      vc <- 3.45
      f1 <- 0.8
      prop.err <- 0.5
    })
    model({
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / f1 / vc * center
      cp <- center / vc
      cp ~ prop(prop.err)
    })
  }
  suppressMessages(ui <- rxode2::rxode(scaledModel))
  expect_equal(intersect(c("ka", "cl", "vc"), pkncaParamMap(ui)$param), "ka")
  expect_message(
    fit <- nlmixr(scaledModel, data = nlmixr2data::theo_sd, est = "pknca"),
    regexp = "NCA initial estimates not applied"
  )
  expect_equal(fit$ui$theta[["cl"]], 1)
  expect_equal(fit$ui$theta[["vc"]], 3.45)

  vModel <- function() {
    ini({
      ka <- 0.45
      cl <- 1
      v <- 3.45
      prop.err <- 0.5
    })
    model({
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ prop(prop.err)
    })
  }
  suppressMessages(
    fit <- nlmixr(vModel, data = nlmixr2data::theo_sd, est = "pknca")
  )
  expect_false(isTRUE(all.equal(fit$ui$theta[["v"]], 3.45)))
})

test_that("pkncaTransformedNames", {
  allowed <- c("depot", "center", "ka", "cl", "vc")
  tn <- function(x) sort(pkncaTransformedNames(x, allowed = allowed))
  # plain uses
  expect_equal(tn(quote(-cl / vc * center)), character())
  expect_equal(tn(quote(ka * depot - cl / vc * center)), character())
  expect_equal(tn(quote(cp <- center / vc)), character())
  expect_equal(tn(quote(cp <- linCmt(ka, cl, vc))), character())
  expect_equal(
    tn(quote(d / dt(center) <- ka * depot - cl / vc * center)),
    character()
  )
  # transformations, scaling and shifts
  expect_equal(
    tn(quote(CL <- exp(cl + 0.75 * log(WT / 70)))),
    sort(c("cl", "WT"))
  )
  expect_true("cl" %in% tn(quote(base::exp(cl) / vc)))
  expect_equal(tn(quote(cl * WT / vc)), sort(c("cl", "WT", "vc")))
  expect_equal(tn(quote(-cl / 70 / vc * center)), sort(c("cl", "vc", "center")))
  expect_equal(
    tn(quote(-cl * (1 / 70) * f1 / vc * center)),
    sort(c("cl", "f1", "vc", "center"))
  )
  expect_equal(
    tn(quote(-cl / F1 / vc * center)),
    sort(c("cl", "F1", "vc", "center"))
  )
  expect_equal(tn(quote(ka * depot - k2 * center)), sort(c("k2", "center")))
  expect_equal(tn(quote(CL <- cl + dcl)), sort(c("cl", "dcl")))
  expect_true("cl" %in% tn(quote(cp <- center / vc + (cl) + bsl)))
  expect_true("cl" %in% tn(quote(cp <- center / vc + -cl)))
})

test_that("ini_transform with logit() parameters", {
  model <- function() {
    ini({
      tka <- 0
      tf <- 0.5
      eta.ka ~ 0.1
      eta.f ~ 0.1
      tcl <- 1
      tvc <- 3
      prop.err <- 0.5
    })
    model({
      ka <- logit(tka + eta.ka)
      fx <- logit(tf, -1, 2) * exp(eta.f)
      cl <- exp(tcl)
      vc <- exp(tvc)
      cp <- linCmt()
      cp ~ prop(prop.err)
    })
  }
  suppressMessages(ui <- rxode2::rxode(model))
  suppressMessages(newmod <- ini_transform(ui, ka = 1.5, fx = 0.3))
  expect_equal(newmod$theta[["tka"]], rxode2::expit(1.5))
  expect_equal(newmod$theta[["tf"]], rxode2::expit(0.3, -1, 2))
})

test_that("pkncaParamMap skips thetas also used elsewhere", {
  model <- function() {
    ini({
      tpop <- 1
      tvc <- 3
      prop.err <- 0.5
    })
    model({
      ka <- tpop
      cl <- tpop * WT
      vc <- tvc
      cp <- linCmt()
      cp ~ prop(prop.err)
    })
  }
  suppressMessages(ui <- rxode2::rxode(model))
  paramMap <- pkncaParamMap(ui)
  expect_false("tpop" %in% paramMap$theta)
  expect_true("tvc" %in% paramMap$theta)
})

test_that("pkncaParamMap keeps a theta that also has an alias", {
  model <- function() {
    ini({
      tka <- 0
      tcl <- 1
      tvc <- 3
      eta.ka ~ 0.1
      prop.err <- 0.5
    })
    model({
      ka <- exp(tka + eta.ka)
      myTka <- tka
      cl <- exp(tcl)
      vc <- exp(tvc)
      cp <- linCmt()
      cp ~ prop(prop.err)
    })
  }
  suppressMessages(ui <- rxode2::rxode(model))
  paramMap <- pkncaParamMap(ui)
  expect_equal(paramMap$param[paramMap$theta == "tka"], "ka")
})

test_that("pkncaAssignedNames", {
  expect_equal(
    pkncaAssignedNames(list(
      quote(a <- 1),
      quote(
        if (x) {
          b <- 2
        } else {
          a <- 3
        }
      ),
      quote(y ~ add(z))
    )),
    c("a", "b", "a")
  )
})

test_that("pkncaSimplifyZeroEta", {
  expect_equal(
    pkncaSimplifyZeroEta(quote(tka * exp(eta.ka)), "eta.ka"),
    quote(tka)
  )
  expect_equal(
    pkncaSimplifyZeroEta(quote(exp(eta.ka) * tka), "eta.ka"),
    quote(tka)
  )
  expect_equal(
    pkncaSimplifyZeroEta(quote(exp(tka + eta.ka)), "eta.ka"),
    quote(exp(tka))
  )
  expect_equal(
    pkncaSimplifyZeroEta(quote(exp(tka - eta.ka)), "eta.ka"),
    quote(exp(tka))
  )
  expect_equal(
    pkncaSimplifyZeroEta(quote(tka / exp(eta.ka)), "eta.ka"),
    quote(tka)
  )
  expect_equal(pkncaSimplifyZeroEta(quote(eta.ka + tka), "eta.ka"), quote(tka))
  expect_equal(pkncaSimplifyZeroEta(quote(tka * WT), "eta.ka"), quote(tka * WT))
  expect_equal(
    pkncaSimplifyZeroEta(quote(tka * exp(-eta.ka)), "eta.ka"),
    quote(tka)
  )
  expect_equal(
    pkncaSimplifyZeroEta(quote(expit(tf) * exp(eta.f)), "eta.f"),
    quote(expit(tf))
  )
  expect_equal(pkncaSimplifyZeroEta(3, "eta.ka"), 3)
  expect_equal(
    pkncaSimplifyZeroEta(quote(tka * exp(0.5 * eta.ka)), "eta.ka"),
    quote(tka)
  )
  expect_equal(
    pkncaSimplifyZeroEta(quote(tka * exp(eta.ka / 2)), "eta.ka"),
    quote(tka)
  )
})


test_that("est='pknca' with covariates and mixed IV/oral dosing (#102)", {
  # Data from the issue: IDs 11-13 IV (CMT 2), IDs 21-23 oral (CMT 1)
  # fmt: skip
  dat <- data.frame(
    ID = rep(c(11, 12, 13, 21, 22, 23), c(7, 7, 8, 7, 7, 7)),
    TIME = c(
      0, 0.05, 0.25, 0.5, 1, 3, 5, 0, 0.05, 0.25, 0.5, 1, 3, 5,
      0, 0.05, 0.25, 0.5, 1, 3, 5, 8, rep(c(0, 0.25, 0.5, 1, 3, 5, 8), 3)
    ),
    DV = c(
      NA, 2017.85, 1323.74, 792.5, 822.72, 36.27, 3.33,
      NA, 1702, 1290.75, 1095.95, 907.6, 125.44, 14.44,
      NA, 1933.04, 1242.43, 661.22, 193.52, 1.75, NA, NA,
      NA, 706.58, 1063.14, 2257.62, 941.33, 629.69, 100,
      NA, 1462.95, 2217.76, 2739.5, 705.3, 108.47, 8.75,
      NA, 211.66, 467.23, 174.24, 153.6, 27.07, 2.81
    ),
    AMT = rep(rep(c(1, 5), each = 3), c(7, 7, 8, 7, 7, 7)) *
      (c(1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
         rep(c(1, 0, 0, 0, 0, 0, 0), 3))),
    EVID = c(1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
             rep(c(1, 0, 0, 0, 0, 0, 0), 3)),
    CMT = c(rep(2, 22), rep(c(1, 2, 2, 2, 2, 2, 2), 3)),
    DOSE = rep(c(1, 5), c(22, 21)),
    ROUTE = rep(c(1, 2), c(22, 21))
  )
  modA <- function() {
    ini({
      tka <- 0.45
      tcl <- -7
      tv <- -8
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      prop.sd <- 0.7
    })
    model({
      if (ROUTE != 1) {
        ka <- exp(tka + eta.ka)
      }
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ prop(prop.sd)
    })
  }
  ctl <- pkncaControl(
    concu = "mg/L",
    doseu = "mg/kg",
    timeu = "hr",
    volumeu = "L/kg"
  )
  suppressMessages(suppressWarnings(
    ret <- nlmixr2est::nlmixr(
      object = modA,
      data = dat,
      est = "pknca",
      control = ctl
    )
  ))
  expect_s3_class(ret, "pkncaEst")
  ncaRes <- as.data.frame(ret$nca)
  ivId <- c(11, 12, 13)
  # tmax only from the extravascular doses
  tmaxRes <- ncaRes[ncaRes$PPTESTCD == "tmax", ]
  expect_setequal(tmaxRes$ID, c(21, 22, 23))
  # vc and cl only from the intravascular doses
  for (nm in c("cmax.dn", "cl.last")) {
    expect_setequal(ncaRes$ID[ncaRes$PPTESTCD == nm], ivId)
    expect_true(all(!is.na(ncaRes$PPORRES[ncaRes$PPTESTCD == nm])))
  }
  # C0 is back-extrapolated for the IV doses (higher than the first observed
  # concentration)
  cmaxIv <- ncaRes[ncaRes$PPTESTCD == "cmax" & ncaRes$ID %in% ivId, ]
  firstConc <- c(2017.85, 1702, 1933.04)
  expect_true(all(cmaxIv$PPORRES[order(cmaxIv$ID)] > firstConc))
  # tcl and tv (as "v") are updated
  expect_equal(
    ret$ui$iniDf$est[ret$ui$iniDf$name == "tcl"],
    log(median(ncaRes$PPORRES[ncaRes$PPTESTCD == "cl.last"]))
  )
  expect_equal(
    ret$ui$iniDf$est[ret$ui$iniDf$name == "tv"],
    log(1 / median(ncaRes$PPORRES[ncaRes$PPTESTCD == "cmax.dn"]))
  )
  # Covariates are kept in the NCA data
  expect_true("ROUTE" %in% names(as.data.frame(ret$nca$data$conc)))

  # The same C0 when the data are not sorted by time
  datShuffle <- dat[rev(seq_len(nrow(dat))), ]
  suppressMessages(suppressWarnings(
    retShuffle <- nlmixr2est::nlmixr(
      object = modA,
      data = datShuffle,
      est = "pknca",
      control = ctl
    )
  ))
  ncaShuffle <- as.data.frame(retShuffle$nca)
  cmaxShuffle <- ncaShuffle[
    ncaShuffle$PPTESTCD == "cmax" & ncaShuffle$ID %in% ivId,
  ]
  expect_equal(
    cmaxShuffle$PPORRES[order(cmaxShuffle$ID)],
    cmaxIv$PPORRES[order(cmaxIv$ID)]
  )

  # Extravascular AUClast starts from zero concentration at the time of dosing
  obs21 <- dat[dat$ID == 21 & dat$EVID == 0, ]
  expect_equal(
    ncaRes$PPORRES[ncaRes$PPTESTCD == "auclast" & ncaRes$ID == 21],
    PKNCA::pk.calc.auc.last(conc = c(0, obs21$DV), time = c(0, obs21$TIME))
  )

  # IV infusions are not back-extrapolated
  datInf <- dat[dat$ROUTE == 1, ]
  datInf$RATE <- ifelse(datInf$EVID == 1, 100, 0)
  suppressMessages(suppressWarnings(
    retInf <- nlmixr2est::nlmixr(
      object = modA,
      data = datInf,
      est = "pknca",
      control = ctl
    )
  ))
  ncaInf <- as.data.frame(retInf$nca)
  cmaxInf <- ncaInf[ncaInf$PPTESTCD == "cmax", ]
  expect_equal(cmaxInf$PPORRES[order(cmaxInf$ID)], firstConc)
  expect_true(all(!is.na(ncaInf$PPORRES[ncaInf$PPTESTCD == "cl.last"])))

  # Extravascular only works without a concentration at the time of dosing
  datOral <- dat[dat$ROUTE == 2, ]
  suppressMessages(suppressWarnings(
    retOral <- nlmixr2est::nlmixr(
      object = modA,
      data = datOral,
      est = "pknca",
      control = ctl
    )
  ))
  ncaOral <- as.data.frame(retOral$nca)
  expect_true(all(
    !is.na(ncaOral$PPORRES[
      ncaOral$PPTESTCD %in% c("tmax", "cmax.dn", "cl.last")
    ])
  ))
})

test_that("pkncaIntervals without extravascular-only doses (#102)", {
  dose <- data.frame(
    ID = c(1, 3, 3),
    TIME = 0,
    pkncaRoute = c("intravascular", "intravascular", "extravascular"),
    pkncaBolus = c(TRUE, TRUE, FALSE)
  )
  intervals <- data.frame(
    ID = c(1, 3),
    start = 0,
    end = Inf,
    cmax = TRUE,
    tmax = TRUE,
    auclast = TRUE
  )
  ret <- pkncaIntervals(
    intervals = intervals,
    dose = dose,
    groupCols = "ID",
    timeCol = "TIME"
  )
  expect_equal(ret$tmax, c(FALSE, FALSE))
  expect_equal(ret$cl.last, c(TRUE, FALSE))
  # No tmax with intravascular only
  ret <- pkncaIntervals(
    intervals = intervals,
    dose = dose[1:2, ],
    groupCols = "ID",
    timeCol = "TIME"
  )
  expect_equal(ret$tmax, c(FALSE, FALSE))
  expect_equal(ret$cl.last, c(TRUE, TRUE))
  # No ka estimate without tmax
  est <- ncaToEst(
    tmax = NULL,
    cmaxdn = c(1, 2, 3),
    cl = c(1, 2, 3),
    control = pkncaControl(),
    unitConversions = c(vss.last = 1, cl.last = 1)
  )
  expect_null(est$ka)
})

test_that("pkncaIntervals route handling (#102)", {
  dose <- data.frame(
    ID = c(1, 2, 3, 3),
    TIME = 0,
    pkncaRoute = c(
      "intravascular",
      "extravascular",
      "intravascular",
      "extravascular"
    ),
    pkncaBolus = c(TRUE, FALSE, TRUE, FALSE)
  )
  intervals <- data.frame(
    ID = c(1, 2, 3),
    start = 0,
    end = Inf,
    cmax = TRUE,
    tmax = TRUE,
    auclast = TRUE,
    half.life = TRUE
  )
  ret <- pkncaIntervals(
    intervals = intervals,
    dose = dose,
    groupCols = "ID",
    timeCol = "TIME"
  )
  # tmax only from the extravascular dose
  expect_equal(ret$tmax, c(FALSE, TRUE, FALSE))
  expect_equal(ret$half.life, c(FALSE, TRUE, FALSE))
  # vc and cl only from the intravascular dose
  expect_equal(ret$cmax.dn, c(TRUE, FALSE, FALSE))
  expect_equal(ret$cl.last, c(TRUE, FALSE, FALSE))
  # Imputation of the start for all but the bolus-only interval
  expect_equal(is.na(ret$impute), c(TRUE, FALSE, FALSE))

  # Intravascular intervals starting from the prior trough are not used for vc
  # when others are available
  doseMulti <- data.frame(
    ID = 1,
    TIME = c(0, 12),
    pkncaRoute = "intravascular",
    pkncaBolus = TRUE,
    pkncaNoC0 = c(FALSE, TRUE)
  )
  intervalsMulti <- data.frame(
    ID = 1,
    start = c(0, 12),
    end = c(12, 24),
    cmax = TRUE,
    auclast = TRUE
  )
  ret <- pkncaIntervals(
    intervals = intervalsMulti,
    dose = doseMulti,
    groupCols = "ID",
    timeCol = "TIME"
  )
  expect_equal(ret$cmax.dn, c(TRUE, FALSE))
  # The AUC from the trough is still used for cl
  expect_equal(ret$cl.last, c(TRUE, TRUE))
  # Without others calculating the parameter, they are used
  intervalsNoAuc <- intervalsMulti
  intervalsNoAuc$auclast <- c(FALSE, TRUE)
  ret <- pkncaIntervals(
    intervals = intervalsNoAuc,
    dose = doseMulti,
    groupCols = "ID",
    timeCol = "TIME"
  )
  expect_equal(ret$cmax.dn, c(TRUE, FALSE))
  expect_equal(ret$cl.last, c(FALSE, TRUE))
  # Without others, they are used (cmax.dn only from the first)
  doseMulti$pkncaNoC0 <- TRUE
  ret <- pkncaIntervals(
    intervals = intervalsMulti,
    dose = doseMulti,
    groupCols = "ID",
    timeCol = "TIME"
  )
  expect_equal(ret$cmax.dn, c(TRUE, FALSE))
  expect_equal(ret$cl.last, c(TRUE, TRUE))
})

test_that("pkncaAddIvC0 (#102)", {
  dose <- data.frame(
    ID = 1,
    TIME = c(0, 12),
    pkncaRoute = "intravascular",
    pkncaBolus = TRUE
  )
  obs <- data.frame(ID = 1, TIME = c(6, 12, 18), DV = c(4, 2, 5))
  ret <- pkncaAddIvC0(
    obs = obs,
    dose = dose,
    groupCols = "ID",
    timeCol = "TIME",
    dvCol = "DV"
  )
  # The trough at the next dose is used for back-extrapolation; the second
  # dose has a concentration at the time of dosing
  expect_equal(ret$TIME, c(0, 6, 12, 18))
  expect_equal(ret$DV, c(8, 4, 2, 5))
  expect_equal(attr(ret, "noC0"), pkncaKey(dose[2, ], c("ID", "TIME")))
  # Only one C0 for simultaneous doses
  dose2 <- rbind(dose[1, ], dose[1, ])
  ret2 <- pkncaAddIvC0(
    obs = obs,
    dose = dose2,
    groupCols = "ID",
    timeCol = "TIME",
    dvCol = "DV"
  )
  expect_equal(sum(ret2$TIME == 0), 1)
  # No C0 when an extravascular dose is at the same time
  dose3 <- rbind(dose[1, ], dose[1, ])
  dose3$pkncaRoute[2] <- "extravascular"
  dose3$pkncaBolus[2] <- FALSE
  ret3 <- pkncaAddIvC0(
    obs = obs,
    dose = dose3,
    groupCols = "ID",
    timeCol = "TIME",
    dvCol = "DV"
  )
  expect_equal(ret3, obs, ignore_attr = TRUE)
  # C0 for a loading bolus with an infusion at the same time
  dose4 <- dose3
  dose4$pkncaRoute[2] <- "intravascular"
  ret4 <- pkncaAddIvC0(
    obs = obs,
    dose = dose4,
    groupCols = "ID",
    timeCol = "TIME",
    dvCol = "DV"
  )
  expect_equal(ret4$TIME, c(0, 6, 12, 18))
  expect_equal(ret4$DV[1], 8)
  # A single concentration uses the first concentration, flagged as not
  # log-linearly back-extrapolated
  obs5 <- data.frame(ID = 1, TIME = c(11.9, 18), DV = c(2, 5))
  ret5 <- pkncaAddIvC0(
    obs = obs5,
    dose = dose,
    groupCols = "ID",
    timeCol = "TIME",
    dvCol = "DV"
  )
  expect_equal(ret5$TIME, c(0, 11.9, 12, 18))
  expect_equal(ret5$DV, c(2, 2, 5, 5))
  expect_equal(attr(ret5, "noC0"), pkncaKey(dose, c("ID", "TIME")))
  # A predose concentration at the first dose is replaced by C0 (missing or
  # not)
  for (dv0 in c(0, NA)) {
    obs4 <- data.frame(ID = 1, TIME = c(0, 6, 12, 18), DV = c(dv0, 4, 2, 5))
    ret4 <- pkncaAddIvC0(
      obs = obs4,
      dose = dose,
      groupCols = "ID",
      timeCol = "TIME",
      dvCol = "DV"
    )
    expect_equal(ret4$TIME, c(0, 6, 12, 18))
    expect_equal(ret4$DV, c(8, 4, 2, 5))
  }
})

test_that("pkncaObsStates (#102)", {
  mOde <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      conc <- center
      cp <- conc / v
      cp ~ add(add.sd)
    })
  }
  mLin <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  mLinTilde <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  suppressMessages({
    expect_equal(pkncaObsStates(rxode2::rxode2(mOde)), "center")
    expect_equal(pkncaObsStates(rxode2::rxode2(mLin)), "central")
    expect_equal(pkncaObsStates(rxode2::rxode2(mLinTilde)), "central")
  })
})

test_that("est='pknca' oral without a CMT column is extravascular (#102)", {
  modelGood <- function() {
    ini({
      tka <- 0.45
      lcl <- 1
      lvc <- 3.45
      prop.err <- 0.5
    })
    model({
      ka <- exp(tka)
      cl <- exp(lcl)
      vc <- exp(lvc)
      cp <- linCmt()
      cp ~ prop(prop.err)
    })
  }
  d <- nlmixr2data::theo_sd
  suppressMessages(
    retCmt <- nlmixr2est::nlmixr(object = modelGood, data = d, est = "pknca")
  )
  d$CMT <- NULL
  suppressMessages(
    retNoCmt <- nlmixr2est::nlmixr(object = modelGood, data = d, est = "pknca")
  )
  feCmt <- setNames(retCmt$ui$iniDf$est, retCmt$ui$iniDf$name)
  feNoCmt <- setNames(retNoCmt$ui$iniDf$est, retNoCmt$ui$iniDf$name)
  expect_false(feNoCmt[["tka"]] == 0.45)
  expect_equal(feNoCmt[c("tka", "lcl", "lvc")], feCmt[c("tka", "lcl", "lvc")])
})

test_that("pkncaCmtOrder with linCmt() and ODEs (#102)", {
  mLinEff <- function() {
    ini({
      tka <- 0
      lcl <- 1
      lv <- 3
      lke0 <- 0
      add.sd <- 0.1
    })
    model({
      ka <- exp(tka)
      cl <- exp(lcl)
      v <- exp(lv)
      ke0 <- exp(lke0)
      cp <- linCmt()
      d / dt(eff) <- ke0 * (cp - eff)
      cp ~ add(add.sd)
    })
  }
  ui <- suppressMessages(rxode2::rxode2(mLinEff))
  expect_equal(pkncaCmtOrder(ui), c("depot", "central", "eff"))
})

test_that("pkncaAutoIntervals (#102)", {
  dose <- data.frame(
    ID = c(1, 2, 2, 2),
    TIME = c(0, 0, 24, 48),
    pkncaRoute = "extravascular"
  )
  obs <- data.frame(
    ID = c(1, 1, 2, 2, 2, 2, 2, 2, 2),
    TIME = c(1, 2, 1, 2, 4, 23.9, 49, 50, 52),
    DV = 1
  )
  ret <- pkncaAutoIntervals(
    obs = obs,
    dose = dose,
    groupCols = "ID",
    timeCol = "TIME",
    dvCol = "DV"
  )
  # Single dose uses the PKNCA defaults
  expect_equal(ret$start[ret$ID == 1], c(0, 0))
  expect_equal(ret$end[ret$ID == 1], c(24, Inf))
  # Multiple doses use intervals with enough concentrations (not only the
  # trough at 23.9 for 24 to 48), ending at the last concentration before the
  # next dose, and the last dosing interval for the last dose
  expect_equal(ret$start[ret$ID == 2], c(0, 48))
  expect_equal(ret$end[ret$ID == 2], c(23.9, 72))
  # AUC (for cl) only from the interval with concentrations covering most of
  # the dosing interval (48 to 72 only has concentrations until 52)
  expect_equal(ret$auclast[ret$ID == 2], c(TRUE, FALSE))
  expect_null(ret$pkncaCoverage)
  # Coverage is decided separately for intravascular intervals
  doseIv <- rbind(
    dose,
    data.frame(ID = 3, TIME = c(0, 24), pkncaRoute = "intravascular")
  )
  obsIv <- rbind(obs, data.frame(ID = 3, TIME = c(1, 2, 4, 25, 26, 28), DV = 1))
  ret <- pkncaAutoIntervals(
    obs = obsIv,
    dose = doseIv,
    groupCols = "ID",
    timeCol = "TIME",
    dvCol = "DV"
  )
  expect_equal(ret$auclast[ret$ID == 2], c(TRUE, FALSE))
  expect_equal(ret$auclast[ret$ID == 3], c(TRUE, TRUE))
  # Peak and trough sampling keeps the intervals
  obsPt <- data.frame(ID = 2, TIME = c(2, 24, 26, 48, 50, 72), DV = 1)
  ret <- pkncaAutoIntervals(
    obs = obsPt,
    dose = dose[dose$ID == 2, ],
    groupCols = "ID",
    timeCol = "TIME",
    dvCol = "DV"
  )
  expect_equal(ret$start, c(0, 24, 48))
  expect_equal(ret$end, c(24, 48, 72))
  # With too few concentrations everywhere, intervals with any concentrations
  # are kept
  obsFew <- data.frame(ID = 2, TIME = c(26, 50), DV = 1)
  ret <- pkncaAutoIntervals(
    obs = obsFew,
    dose = dose[dose$ID == 2, ],
    groupCols = "ID",
    timeCol = "TIME",
    dvCol = "DV"
  )
  expect_equal(ret$start, c(24, 48))
})

test_that("pkncaIntervals first usable cmax.dn per route (#102)", {
  dose <- data.frame(
    ID = 1,
    TIME = c(0, 24, 48, 72),
    pkncaRoute = c(
      "intravascular",
      "intravascular",
      "extravascular",
      "extravascular"
    ),
    pkncaBolus = c(TRUE, TRUE, FALSE, FALSE),
    pkncaNoC0 = c(TRUE, FALSE, FALSE, FALSE)
  )
  intervals <- data.frame(
    ID = 1,
    start = c(0, 24, 48, 72),
    end = c(24, 48, 72, 96),
    cmax = TRUE,
    tmax = TRUE,
    auclast = TRUE
  )
  ret <- pkncaIntervals(
    intervals = intervals,
    dose = dose,
    groupCols = "ID",
    timeCol = "TIME"
  )
  # The first IV dose has no log-linear C0, so the second IV dose is used
  expect_equal(ret$cmax.dn, c(FALSE, TRUE, FALSE, FALSE))
  # Without IV doses, the first oral dose is used
  ret <- pkncaIntervals(
    intervals = intervals[3:4, ],
    dose = dose[3:4, ],
    groupCols = "ID",
    timeCol = "TIME"
  )
  expect_equal(ret$cmax.dn, c(TRUE, FALSE))
})

test_that("pkncaCollapseDose (#102)", {
  dose <- data.frame(
    ID = c(1, 1, 1, 2),
    TIME = c(0, 0, 12, 0),
    AMT = c(1, 2, 3, 4),
    pkncaRoute = c(
      "intravascular",
      "extravascular",
      "intravascular",
      "intravascular"
    )
  )
  ret <- pkncaCollapseDose(
    dose = dose,
    groupCols = "ID",
    timeCol = "TIME",
    amtCol = "AMT"
  )
  expect_equal(ret$TIME, c(0, 12, 0))
  expect_equal(ret$AMT, c(3, 3, 4))
  expect_equal(
    ret$pkncaRoute,
    c("extravascular", "intravascular", "intravascular")
  )
  expect_equal(
    pkncaCollapseDose(
      dose = dose[-1, ],
      groupCols = "ID",
      timeCol = "TIME",
      amtCol = "AMT"
    ),
    dose[-1, ]
  )
})

test_that("est='pknca' with simultaneous IV and oral doses (#102)", {
  mod <- function() {
    ini({
      tka <- 0.45
      lcl <- 1
      lvc <- 3.45
      prop.err <- 0.5
    })
    model({
      ka <- exp(tka)
      cl <- exp(lcl)
      vc <- exp(lvc)
      cp <- linCmt()
      cp ~ prop(prop.err)
    })
  }
  d <- nlmixr2data::theo_sd
  extra <- d[d$ID == 1 & d$EVID != 0, ]
  extra$CMT <- 2
  d <- rbind(d, extra)
  d <- d[order(d$ID, d$TIME, -d$EVID), ]
  suppressMessages(suppressWarnings(
    ret <- nlmixr2est::nlmixr(object = mod, data = d, est = "pknca")
  ))
  expect_s3_class(ret, "pkncaEst")
  ncaRes <- as.data.frame(ret$nca)
  # The combined doses for ID 1 are one PKNCA dose
  doseNca <- as.data.frame(ret$nca$data$dose$data)
  expect_equal(sum(doseNca$ID == 1), 1)
  # ID 1 (IV and oral at the same time) is not used for ka since other
  # subjects only have oral doses
  expect_false(1 %in% ncaRes$ID[ncaRes$PPTESTCD == "tmax"])
  expect_setequal(
    ncaRes$ID[ncaRes$PPTESTCD == "tmax"],
    setdiff(unique(d$ID), 1)
  )
})

test_that("est='pknca' multiple-dose oral without dose-time conc (#102)", {
  mod <- function() {
    ini({
      tka <- 0
      lcl <- 1
      lv <- 3
      add.sd <- 0.1
    })
    model({
      ka <- exp(tka)
      cl <- exp(lcl)
      v <- exp(lv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  d <- nlmixr2data::Oral_1CPT
  d <- d[d$SD == 0 & d$ID <= 40, ]
  suppressMessages(suppressWarnings(
    ret <- nlmixr2est::nlmixr(object = mod, data = d, est = "pknca")
  ))
  fe <- setNames(ret$ui$iniDf$est, ret$ui$iniDf$name)
  trueCl <- log(median(d$CL[!duplicated(d$ID)]))
  trueV <- log(median(d$V[!duplicated(d$ID)]))
  expect_lt(abs(fe[["lcl"]] - trueCl), 0.5)
  expect_lt(abs(fe[["lv"]] - trueV), 0.5)
})
