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

  # The same C0 when the data are not sorted by time
  datShuffle <- dat[rev(seq_len(nrow(dat))), ]
  suppressMessages(suppressWarnings(
    retShuffle <- nlmixr2est::nlmixr(object = modA, data = datShuffle, est = "pknca", control = ctl)
  ))
  ncaShuffle <- as.data.frame(retShuffle$nca)
  cmaxShuffle <- ncaShuffle[ncaShuffle$PPTESTCD == "cmax" & ncaShuffle$ID %in% ivId, ]
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

test_that("pkncaIntervals without extravascular-only doses (#102)", {
  dose <- data.frame(
    ID = c(1, 3, 3),
    TIME = 0,
    pkncaRoute = c("intravascular", "intravascular", "extravascular"),
    pkncaBolus = c(TRUE, TRUE, FALSE)
  )
  intervals <- data.frame(ID = c(1, 3), start = 0, end = Inf, cmax = TRUE, tmax = TRUE, auclast = TRUE)
  ret <- pkncaIntervals(intervals = intervals, dose = dose, groupCols = "ID", timeCol = "TIME")
  expect_equal(ret$tmax, c(FALSE, FALSE))
  expect_equal(ret$cl.last, c(TRUE, FALSE))
  # No tmax with intravascular only
  ret <- pkncaIntervals(intervals = intervals, dose = dose[1:2, ], groupCols = "ID", timeCol = "TIME")
  expect_equal(ret$tmax, c(FALSE, FALSE))
  expect_equal(ret$cl.last, c(TRUE, TRUE))
  # No ka estimate without tmax
  est <- ncaToEst(
    tmax = NULL, cmaxdn = c(1, 2, 3), cl = c(1, 2, 3),
    control = pkncaControl(),
    unitConversions = c(vss.last = 1, cl.last = 1)
  )
  expect_null(est$ka)
  expect_equal(est$v, est$vc)
})

test_that("pkncaIntervals route handling (#102)", {
  dose <- data.frame(
    ID = c(1, 2, 3, 3),
    TIME = 0,
    pkncaRoute = c("intravascular", "extravascular", "intravascular", "extravascular"),
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
  ret <- pkncaIntervals(intervals = intervals, dose = dose, groupCols = "ID", timeCol = "TIME")
  # tmax only from the extravascular dose
  expect_equal(ret$tmax, c(FALSE, TRUE, FALSE))
  expect_equal(ret$half.life, c(FALSE, TRUE, FALSE))
  # vc and cl only from the intravascular dose
  expect_equal(ret$cmax.dn, c(TRUE, FALSE, FALSE))
  expect_equal(ret$cl.last, c(TRUE, FALSE, FALSE))
  # Imputation of the start for all but the bolus-only interval
  expect_equal(is.na(ret$impute), c(TRUE, FALSE, FALSE))

  # Intravascular intervals starting from the prior trough are not used for vc
  # and cl when others are available
  doseMulti <- data.frame(
    ID = 1,
    TIME = c(0, 12),
    pkncaRoute = "intravascular",
    pkncaBolus = TRUE,
    pkncaNoC0 = c(FALSE, TRUE)
  )
  intervalsMulti <- data.frame(ID = 1, start = c(0, 12), end = c(12, 24), cmax = TRUE, auclast = TRUE)
  ret <- pkncaIntervals(intervals = intervalsMulti, dose = doseMulti, groupCols = "ID", timeCol = "TIME")
  expect_equal(ret$cmax.dn, c(TRUE, FALSE))
  expect_equal(ret$cl.last, c(TRUE, FALSE))
  doseMulti$pkncaNoC0 <- TRUE
  ret <- pkncaIntervals(intervals = intervalsMulti, dose = doseMulti, groupCols = "ID", timeCol = "TIME")
  expect_equal(ret$cmax.dn, c(TRUE, TRUE))
})

test_that("pkncaAddIvC0 (#102)", {
  dose <- data.frame(
    ID = 1,
    TIME = c(0, 12),
    pkncaRoute = "intravascular",
    pkncaBolus = TRUE
  )
  obs <- data.frame(ID = 1, TIME = c(6, 12, 18), DV = c(4, 2, 5))
  ret <- pkncaAddIvC0(obs = obs, dose = dose, groupCols = "ID", timeCol = "TIME", dvCol = "DV")
  # The trough at the next dose is used for back-extrapolation; the second
  # dose has a concentration at the time of dosing
  expect_equal(ret$TIME, c(0, 6, 12, 18))
  expect_equal(ret$DV, c(8, 4, 2, 5))
  expect_equal(attr(ret, "noC0"), pkncaKey(dose[2, ], c("ID", "TIME")))
  # Only one C0 for simultaneous doses
  dose2 <- rbind(dose[1, ], dose[1, ])
  ret2 <- pkncaAddIvC0(obs = obs, dose = dose2, groupCols = "ID", timeCol = "TIME", dvCol = "DV")
  expect_equal(sum(ret2$TIME == 0), 1)
  # No C0 when an extravascular dose is at the same time
  dose3 <- rbind(dose[1, ], dose[1, ])
  dose3$pkncaRoute[2] <- "extravascular"
  dose3$pkncaBolus[2] <- FALSE
  ret3 <- pkncaAddIvC0(obs = obs, dose = dose3, groupCols = "ID", timeCol = "TIME", dvCol = "DV")
  expect_equal(ret3, obs, ignore_attr = TRUE)
  # A single concentration uses the first concentration, flagged as not
  # log-linearly back-extrapolated
  obs5 <- data.frame(ID = 1, TIME = c(11.9, 18), DV = c(2, 5))
  ret5 <- pkncaAddIvC0(obs = obs5, dose = dose, groupCols = "ID", timeCol = "TIME", dvCol = "DV")
  expect_equal(ret5$TIME, c(0, 11.9, 12, 18))
  expect_equal(ret5$DV, c(2, 2, 5, 5))
  expect_equal(attr(ret5, "noC0"), pkncaKey(dose, c("ID", "TIME")))
  # A predose concentration at the first dose is replaced by C0 (missing or
  # not)
  for (dv0 in c(0, NA)) {
    obs4 <- data.frame(ID = 1, TIME = c(0, 6, 12, 18), DV = c(dv0, 4, 2, 5))
    ret4 <- pkncaAddIvC0(obs = obs4, dose = dose, groupCols = "ID", timeCol = "TIME", dvCol = "DV")
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
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
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
  suppressMessages({
    expect_equal(pkncaObsStates(rxode2::rxode2(mOde)), "center")
    expect_equal(pkncaObsStates(rxode2::rxode2(mLin)), "central")
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
  suppressMessages(retCmt <- nlmixr2est::nlmixr(object = modelGood, data = d, est = "pknca"))
  d$CMT <- NULL
  suppressMessages(retNoCmt <- nlmixr2est::nlmixr(object = modelGood, data = d, est = "pknca"))
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
      d/dt(eff) <- ke0 * (cp - eff)
      cp ~ add(add.sd)
    })
  }
  ui <- suppressMessages(rxode2::rxode2(mLinEff))
  expect_equal(pkncaCmtOrder(ui), c("depot", "central", "eff"))
})

test_that("pkncaAutoIntervals (#102)", {
  dose <- data.frame(ID = c(1, 2, 2, 2), TIME = c(0, 0, 24, 48))
  obs <- data.frame(
    ID = c(1, 1, 2, 2, 2, 2, 2, 2, 2),
    TIME = c(1, 2, 1, 2, 4, 23.9, 49, 50, 52),
    DV = 1
  )
  ret <- pkncaAutoIntervals(obs = obs, dose = dose, groupCols = "ID", timeCol = "TIME", dvCol = "DV")
  # Single dose uses the PKNCA defaults
  expect_equal(ret$start[ret$ID == 1], c(0, 0))
  expect_equal(ret$end[ret$ID == 1], c(24, Inf))
  # Multiple doses use intervals with enough concentrations (not only the
  # trough at 23.9 for 24 to 48), ending at the last concentration before the
  # next dose, the last dosing interval for the last dose, and the half-life
  # after the last dose
  expect_equal(ret$start[ret$ID == 2], c(0, 48, 48))
  expect_equal(ret$end[ret$ID == 2], c(23.9, 72, Inf))
  expect_equal(ret$auclast[ret$ID == 2], c(TRUE, TRUE, FALSE))
  # Cmax only from the first dose
  expect_equal(ret$cmax[ret$ID == 2], c(TRUE, FALSE, FALSE))
  expect_equal(ret$half.life[ret$ID == 2], c(FALSE, FALSE, TRUE))
})

test_that("pkncaCollapseDose (#102)", {
  dose <- data.frame(
    ID = c(1, 1, 1, 2),
    TIME = c(0, 0, 12, 0),
    AMT = c(1, 2, 3, 4),
    pkncaRoute = c("intravascular", "extravascular", "intravascular", "intravascular")
  )
  ret <- pkncaCollapseDose(dose = dose, groupCols = "ID", timeCol = "TIME", amtCol = "AMT")
  expect_equal(ret$TIME, c(0, 12, 0))
  expect_equal(ret$AMT, c(3, 3, 4))
  expect_equal(ret$pkncaRoute, c("extravascular", "intravascular", "intravascular"))
  expect_equal(pkncaCollapseDose(dose = dose[-1, ], groupCols = "ID", timeCol = "TIME", amtCol = "AMT"), dose[-1, ])
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
})

test_that("est='pknca' multiple-dose extravascular without concentrations at dosing (#102)", {
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
