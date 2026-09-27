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
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  .theo <- nlmixr2data::theo_sd
  .cens <- .theo
  .cens$CENS <- ifelse(.cens$DV < 1 & .cens$EVID == 0, 1, 0)
  .cens$DV[.cens$CENS == 1] <- 1

  .export <- function(data, name, model=one.cmt, ...) {
    suppressMessages(suppressWarnings(
      nlmixr2est::nlmixr(model, data=data, est="nonmem",
                         control=nonmemControl(runCommand=NA,
                                               modelName=name, ...))))
    .dir <- paste0(name, "-nonmem")
    list(ctl=readLines(file.path(.dir, paste0(name, ".nmctl"))),
         data=utils::read.csv(file.path(.dir, paste0(name, ".csv"))))
  }

  test_that("CENS uses M3 censoring with LAPLACIAN in NONMEM (#92)", {
    .r <- .export(.cens, "m3")
    expect_true("$INPUT ID TIME EVID AMT DV CMT CENS RXROW" %in% .r$ctl)
    expect_true(any(grepl("F_FLAG = 1", .r$ctl, fixed=TRUE)))
    expect_true(any(grepl("Y = PHI(CENS*(DV-IPRED)/W)", .r$ctl, fixed=TRUE)))
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
    expect_true(any(grepl("Y = (CUM1-CUM2)/PHI(-CENS*(LIMIT-IPRED)/W)",
                          .r$ctl, fixed=TRUE)))
    expect_true(any(grepl("^\\$ESTIMATION METHOD=1 LAPLACIAN INTER ", .r$ctl)))
    # missing limits are NONMEM's infinity
    expect_equal(.r$data$LIMIT[.r$data$CENS == 1], rep(0, 16))
    expect_true(all(.r$data$LIMIT[.r$data$CENS == 0] == -1000000))
  })

  test_that("posthoc with censoring uses the Laplacian method (#92)", {
    .d <- .cens
    .d$LIMIT <- ifelse(.d$CENS == 1, 0, NA)
    .r <- .export(.d, "m4post", est="posthoc")
    expect_true(any(grepl("^\\$ESTIMATION METHOD=1 LAPLACIAN INTER MAXEVALS=0 POSTHOC ",
                          .r$ctl)))
  })

  test_that("LIMIT without CENS uses M2 in NONMEM (#92)", {
    .d <- .theo
    .d$LIMIT <- ifelse(.d$EVID == 0, 0.5, NA)
    .r <- .export(.d, "m2")
    expect_true("$INPUT ID TIME EVID AMT DV CMT CENS LIMIT RXROW" %in% .r$ctl)
    expect_true(any(grepl("Y = Y/PHI(ABS(LIMIT-IPRED)/W)", .r$ctl, fixed=TRUE)))
    expect_true(all(.r$data$CENS == 0))
    expect_true(all(.r$data$LIMIT[.r$data$EVID == 0] == 0.5))
  })

  test_that("right censoring keeps finite limits and missing limits are infinite (#92)", {
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
    expect_error(suppressMessages(bblDatToNonmem(one.cmt, .d)),
                 "between -1000000 and 1000000")
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
      assign("control", nonmemControl(), envir=.ui)
      .env <- new.env(parent=emptyenv())
      .d <- suppressMessages(suppressWarnings(bblDatToNonmem(.ui, data, env=.env)))
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

  test_that("censoring with a transformed endpoint is refused in NONMEM (#92)", {
    .m <- rxode2::model(one.cmt, cp ~ lnorm(add.sd))
    expect_error(
      suppressMessages(
        nlmixr2est::nlmixr(.m, data=.cens, est="nonmem",
                           control=nonmemControl(runCommand=NA,
                                                 modelName="lnorm"))),
      "transformed endpoints"
    )
  })

})
