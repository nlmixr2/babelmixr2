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
