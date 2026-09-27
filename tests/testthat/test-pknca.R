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

test_that("est='pknca' with covariates and mixed IV/oral dosing (#102)", {
  dat <- structure(list(ID = c(11, 11, 11, 11, 11, 11, 11, 12, 12, 12,
                               12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 21, 21, 21, 21,
                               21, 21, 21, 22, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23),
                        TIME = c(0, 0.05, 0.25, 0.5, 1, 3, 5, 0, 0.05, 0.25, 0.5,
                                 1, 3, 5, 0, 0.05, 0.25, 0.5, 1, 3, 5, 8, 0, 0.25, 0.5, 1, 3,
                                 5, 8, 0, 0.25, 0.5, 1, 3, 5, 8, 0, 0.25, 0.5, 1, 3, 5, 8),
                        DV = c(NA,2017.85, 1323.74, 792.5, 822.72, 36.27, 3.33, NA, 1702, 1290.75,
                               1095.95, 907.6, 125.44, 14.44, NA, 1933.04, 1242.43, 661.22,
                               193.52, 1.75, NA, NA, NA, 706.58, 1063.14, 2257.62, 941.33, 629.69,
                               100, NA, 1462.95, 2217.76, 2739.5, 705.3, 108.47, 8.75, NA, 211.66,
                               467.23, 174.24, 153.6, 27.07, 2.81),
                        AMT = c(1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0,
                                0, 0, 5, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0),
                        EVID = c(1,0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                                 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
                        CMT = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 2,
                                2, 2, 2, 2, 2),
                        DOSE = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
                        ROUTE = c(1, 1, 1, 1, 1, 1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)),
                   row.names = c(NA, -43L), class = "data.frame")
  modA <- function() {
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
      if (ROUTE != 1) ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ prop(prop.sd)
    })
  }
  ctl <- pkncaControl(concu = "mg/L", doseu = "mg/kg", timeu = "hr", volumeu = "L/kg")
  suppressMessages(suppressWarnings(
    ret <- nlmixr2est::nlmixr(object = modA, data = dat, est = "pknca", control = ctl)
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

  # IV infusions are not back-extrapolated
  datInf <- dat[dat$ROUTE == 1, ]
  datInf$RATE <- ifelse(datInf$EVID == 1, 100, 0)
  suppressMessages(suppressWarnings(
    retInf <- nlmixr2est::nlmixr(object = modA, data = datInf, est = "pknca", control = ctl)
  ))
  ncaInf <- as.data.frame(retInf$nca)
  cmaxInf <- ncaInf[ncaInf$PPTESTCD == "cmax", ]
  expect_equal(cmaxInf$PPORRES[order(cmaxInf$ID)], firstConc)
  expect_true(all(!is.na(ncaInf$PPORRES[ncaInf$PPTESTCD == "cl.last"])))

  # Extravascular only works without a concentration at the time of dosing
  datOral <- dat[dat$ROUTE == 2, ]
  suppressMessages(suppressWarnings(
    retOral <- nlmixr2est::nlmixr(object = modA, data = datOral, est = "pknca", control = ctl)
  ))
  ncaOral <- as.data.frame(retOral$nca)
  expect_true(all(!is.na(ncaOral$PPORRES[ncaOral$PPTESTCD %in% c("tmax", "cmax.dn", "cl.last")])))
})
