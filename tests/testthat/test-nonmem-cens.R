withr::with_tempdir({
  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  .theo <- nlmixr2data::theo_sd
  .cens <- .theo
  .cens$CENS <- ifelse(.cens$DV < 1 & .cens$EVID == 0, 1, 0)
  .cens$DV[.cens$CENS == 1] <- 1

  .export <- function(data, name, model = one.cmt, ...) {
    suppressMessages(suppressWarnings(
      nlmixr2est::nlmixr(
        model,
        data = data,
        est = "nonmem",
        control = nonmemControl(runCommand = NA, modelName = name, ...)
      )
    ))
    .dir <- paste0(name, "-nonmem")
    list(
      ctl = readLines(file.path(.dir, paste0(name, ".nmctl"))),
      data = utils::read.csv(file.path(.dir, paste0(name, ".csv")))
    )
  }

  test_that("CENS uses M3 censoring with LAPLACIAN in NONMEM (#92)", {
    .r <- .export(.cens, "m3")
    expect_true("$INPUT ID TIME EVID AMT DV CMT CENS RXROW" %in% .r$ctl)
    expect_true(any(grepl("F_FLAG = 1", .r$ctl, fixed = TRUE)))
    expect_true(any(grepl("Y = PHI(CENS*(DV-IPRED)/W)", .r$ctl, fixed = TRUE)))
    expect_false(any(grepl("LIMIT", .r$ctl)))
    expect_true(any(grepl("^\\$ESTIMATION METHOD=1 LAPLACIAN INTER ", .r$ctl)))
    expect_equal(sum(.r$data$CENS), 16)
  })

  test_that("a LIMIT with no finite values is dropped (#92)", {
    .d <- .cens
    .d$LIMIT <- NA
    .r <- .export(.d, "m3na")
    expect_true("$INPUT ID TIME EVID AMT DV CMT CENS RXROW" %in% .r$ctl)
    expect_false(any(names(.r$data) == "LIMIT"))
  })

  test_that("CENS with LIMIT uses M4 censoring in NONMEM (#92)", {
    .d <- .cens
    .d$LIMIT <- ifelse(.d$CENS == 1, 0, NA)
    .r <- .export(.d, "m4")
    expect_true("$INPUT ID TIME EVID AMT DV CMT CENS LIMIT RXROW" %in% .r$ctl)
    expect_true(any(grepl(
      "Y = (CUM1-CUM2)/CUM3",
      .r$ctl,
      fixed = TRUE
    )))
    expect_true(any(grepl("^\\$ESTIMATION METHOD=1 LAPLACIAN INTER ", .r$ctl)))
    # missing limits are NONMEM's infinity
    expect_equal(.r$data$LIMIT[.r$data$CENS == 1], rep(0, 16))
    expect_true(all(.r$data$LIMIT[.r$data$CENS == 0] == -1000000))
  })

  test_that("posthoc with censoring uses the Laplacian method (#92)", {
    .d <- .cens
    .d$LIMIT <- ifelse(.d$CENS == 1, 0, NA)
    .r <- .export(.d, "m4post", est = "posthoc")
    expect_true(any(grepl(
      "^\\$ESTIMATION METHOD=1 LAPLACIAN INTER MAXEVALS=0 POSTHOC ",
      .r$ctl
    )))
  })

  test_that("LIMIT without CENS uses M2 in NONMEM (#92)", {
    .d <- .theo
    .d$LIMIT <- ifelse(.d$EVID == 0, 0.5, NA)
    .r <- .export(.d, "m2")
    expect_true("$INPUT ID TIME EVID AMT DV CMT CENS LIMIT RXROW" %in% .r$ctl)
    expect_true(any(grepl(
      "Y = Y/PHI(ABS(LIMIT-IPRED)/W)",
      .r$ctl,
      fixed = TRUE
    )))
    expect_true(all(.r$data$CENS == 0))
    expect_true(all(.r$data$LIMIT[.r$data$EVID == 0] == 0.5))
  })

  test_that("CENS=0 with a finite LIMIT uses M2 next to M3 (#92)", {
    .d <- .cens
    .d$LIMIT <- ifelse(.d$CENS == 0 & .d$EVID == 0, 0.5, NA)
    .r <- .export(.d, "m2m3")
    expect_true("$INPUT ID TIME EVID AMT DV CMT CENS LIMIT RXROW" %in% .r$ctl)
    expect_true(any(grepl(
      "Y = Y/PHI(ABS(LIMIT-IPRED)/W)",
      .r$ctl,
      fixed = TRUE
    )))
    .obs <- .r$data$EVID == 0
    expect_true(all(.r$data$LIMIT[.obs & .r$data$CENS == 0] == 0.5))
    expect_true(all(abs(.r$data$LIMIT[.r$data$CENS == 1]) == 1000000))
  })

  test_that("right censoring keeps finite limits, missing are infinite (#92)", {
    .d <- .cens
    .d$CENS[.d$CENS == 1] <- -1
    .d$LIMIT <- NA
    .d$LIMIT[which(.d$CENS == -1)[1]] <- 10
    .r <- .export(.d, "right")
    .l <- .r$data$LIMIT[.r$data$CENS == -1]
    expect_equal(.l[1], 10)
    expect_true(all(abs(.l[-1]) == 1000000))
  })

  test_that("finite limits that are NONMEM's infinity are refused (#92)", {
    .d <- .cens
    .d$LIMIT <- ifelse(.d$CENS == 1, -2000000, NA)
    expect_error(
      suppressMessages(bblDatToNonmem(one.cmt, .d)),
      "between -1000000 and 1000000"
    )
  })

  test_that("no censored values means no censoring in NONMEM (#92)", {
    .d <- .theo
    .d$CENS <- 0
    .r <- .export(.d, "none")
    expect_true("$INPUT ID TIME EVID AMT DV CMT RXROW" %in% .r$ctl)
    expect_false(any(grepl("F_FLAG", .r$ctl)))
    expect_true(any(grepl("^\\$ESTIMATION METHOD=1 INTER ", .r$ctl)))
  })

  test_that("censored and limited observations adjust the objective (#92)", {
    .nFlag <- function(data) {
      .ui <- rxode2::rxUiDecompress(rxode2::rxode2(one.cmt))
      assign("control", nonmemControl(), envir = .ui)
      .env <- new.env(parent = emptyenv())
      .d <- suppressMessages(suppressWarnings(bblDatToNonmem(
        .ui,
        data,
        env = .env
      )))
      .nonmemFormatCensData(.d, .ui)
      rxode2::rxGetControl(.ui, ".nFlag", NA_integer_)
    }
    expect_equal(.nFlag(.cens), 16L)
    .d <- .cens
    .d$LIMIT <- ifelse(.d$EVID == 0, 0.5, NA)
    # every observation is M2, M3 or M4
    expect_equal(.nFlag(.d), sum(.d$EVID == 0))
    .d <- .theo
    .d$CENS <- 0
    expect_equal(.nFlag(.d), 0L)
  })

  test_that("ITS with censoring uses LAPLACE (#92)", {
    .r <- .export(.cens, "m3its", est = "its")
    expect_true(any(grepl(
      "^\\$ESTIMATION METHOD=ITS LAPLACIAN INTERACTION ",
      .r$ctl
    )))
  })

  test_that("censored likelihoods are floored in NONMEM (#92)", {
    .d <- .cens
    .d$LIMIT <- ifelse(.d$CENS == 1, 0, NA)
    .r <- .export(.d, "floor")
    expect_equal(
      sum(grepl("IF (Y .LT. 1.0E-30) Y = 1.0E-30", .r$ctl, fixed = TRUE)),
      2L
    )
    expect_true(any(grepl(
      "IF (CUM3 .LT. 1.0E-30) CUM3 = 1.0E-30",
      .r$ctl,
      fixed = TRUE
    )))
  })

  test_that("censoring works with multiple endpoints (#92)", {
    pk.turnover.emax3 <- function() {
      ini({
        tktr <- log(1)
        tka <- log(1)
        tcl <- log(0.1)
        tv <- log(10)
        eta.ktr ~ 1
        eta.ka ~ 1
        eta.cl ~ 2
        eta.v ~ 1
        prop.err <- 0.1
        pkadd.err <- 0.1
        temax <- logit(0.8)
        tec50 <- log(0.5)
        tkout <- log(0.05)
        te0 <- log(100)
        eta.emax ~ .5
        eta.ec50 ~ .5
        eta.kout ~ .5
        eta.e0 ~ .5
        pdadd.err <- 10
      })
      model({
        ktr <- exp(tktr + eta.ktr)
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        emax <- expit(temax + eta.emax)
        ec50 <- exp(tec50 + eta.ec50)
        kout <- exp(tkout + eta.kout)
        e0 <- exp(te0 + eta.e0)
        DCP <- center / v
        PD <- 1 - emax * DCP / (ec50 + DCP)
        effect(0) <- e0
        kin <- e0 * kout
        d / dt(depot) <- -ktr * depot
        d / dt(gut) <- ktr * depot - ka * gut
        d / dt(center) <- ka * gut - cl / v * center
        d / dt(effect) <- kin * PD - kout * effect
        cp <- center / v
        cp ~ prop(prop.err) + add(pkadd.err)
        effect ~ add(pdadd.err) | pca
      })
    }
    .d <- nlmixr2data::warfarin
    .d$CENS <- ifelse(.d$dvid == "cp" & .d$evid == 0 & .d$dv < 2, 1, 0)
    .d$dv[.d$CENS == 1] <- 2
    .d$LIMIT <- ifelse(.d$CENS == 1, 0, NA)
    .r <- .export(.d, "multi", model = pk.turnover.emax3)
    expect_true(
      "$INPUT ID TIME EVID AMT DV CMT DVID CENS LIMIT RXROW" %in% .r$ctl
    )
    expect_equal(sum(.r$data$CENS), sum(.d$CENS))
    expect_true(all(.r$data$DVID[.r$data$CENS == 1] == 1))
  })

  test_that("F_FLAG rows are dropped from the NONMEM PRED check (#92)", {
    .fit <- list(env = new.env(parent = emptyenv()))
    expect_equal(.nonmemFlagRows(.fit), integer(0))
    .fit$env$nonmemData <- data.frame(
      EVID = c(1, 0, 0, 0, 0),
      CENS = c(0, 1, 0, 0, -1),
      LIMIT = c(-1000000, -1000000, 0.5, -1000000, 1000000),
      nlmixrRowNums = 1:5
    )
    expect_equal(.nonmemFlagRows(.fit), c(2L, 3L, 5L))
    .fit$env$nonmemData$LIMIT <- NULL
    expect_equal(.nonmemFlagRows(.fit), c(2L, 5L))
  })

  test_that("missing CENS values are not censored in NONMEM (#92)", {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(one.cmt))
    .d <- data.frame(
      ID = 1,
      TIME = 0:3,
      EVID = c(1, 0, 0, 0),
      AMT = c(1, 0, 0, 0),
      DV = c(NA, 1, 2, 3),
      CMT = 1,
      CENS = c(NA, 1, NA, 0),
      nlmixrRowNums = 1:4
    )
    .r <- .nonmemFormatCensData(.d, .ui)
    expect_equal(.r$CENS, c(0, 1, 0, 0))
    expect_equal(rxode2::rxGetControl(.ui, ".nFlag", NA_integer_), 1L)
  })

  test_that("censoring with a transformed endpoint is refused (#92)", {
    .m <- rxode2::model(one.cmt, cp ~ lnorm(add.sd))
    expect_error(
      suppressMessages(
        nlmixr2est::nlmixr(
          .m,
          data = .cens,
          est = "nonmem",
          control = nonmemControl(runCommand = NA, modelName = "lnorm")
        )
      ),
      "transformed endpoints"
    )
  })
})
