.hasLinCmtMicro <- function() {
  "linCmtMicro" %in% getNamespaceExports("rxode2")
}

.linCmtOne <- function() {
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
    alag(depot) <- 0.1
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

.linCmtTwoIv <- function() {
  ini({
    tcl <- log(4)
    tv <- log(40)
    tq <- log(8)
    tvp <- log(80)
    eta.cl ~ 0.1
    eta.v ~ 0.1
    eta.q ~ 0.1
    eta.vp ~ 0.1
    prop.sd <- 0.1
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    q <- exp(tq + eta.q)
    vp <- exp(tvp + eta.vp)
    k <- 2
    cp <- linCmt()
    cp ~ prop(prop.sd)
  })
}

.linCmtEffect <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    tke0 <- log(0.5)
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    eta.ke0 ~ 0.1
    add.sd <- 0.7
    add.e <- 0.3
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    ke0 <- exp(tke0 + eta.ke0)
    cp <- linCmt()
    d/dt(ce) <- ke0 * (cp - ce)
    cp ~ add(add.sd)
    ce ~ add(add.e)
  })
}

.linCmtExport <- function(model, data, est, control) {
  withr::with_tempdir({
    suppressMessages(nlmixr2est::nlmixr2(model, data, est, control))
    if (est == "nonmem") {
      readLines(file.path("x-nonmem", "x.nmctl"))
    } else {
      c(readLines("x-monolix.txt"), readLines("x-monolix.mlxtran"))
    }
  })
}

.theoLin <- function() {
  .d <- nlmixr2data::theo_sd
  .d$EVID <- ifelse(.d$EVID == 0, 0L, 1L)
  .d
}

.ivLin <- function() {
  .d <- .theoLin()
  .d$CMT <- 1
  .d
}

test_that("state dependence of the model lines", {
  .e <- list(quote(ka <- exp(tka)),
             quote(d/dt(central) <- -k * central),
             quote(cp <- central / v),
             quote(eff <- cp * 2),
             quote(tt <- t * 2),
             quote(f(depot) <- 0.5),
             quote(z <- ka + 1))
  .d <- .bblLinCmtStateDep(.e, c("depot", "central"))
  expect_equal(.d$dep, c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(.d$cmt, c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE))
  expect_true(all(c("cp", "eff", "tt") %in% .d$vars))
})

test_that("linCmt() is translated to NONMEM's closed-form ADVAN", {
  skip_on_cran()
  skip_if_not(.hasLinCmtMicro())
  .ctl <- nonmemControl(runCommand=NA, modelName="x")

  .nm <- .linCmtExport(.linCmtOne, .theoLin(), "nonmem", .ctl)
  expect_true("$SUBROUTINES ADVAN2 TRANS1" %in% .nm)
  expect_false(any(grepl("^\\$MODEL|^\\$DES|DADT", .nm)))
  expect_true(any(grepl("^  K=CL/V", .nm)))
  expect_true(any(grepl("^  ALAG1=0.1", .nm)))
  expect_true(any(grepl("A\\(2\\)/\\(V\\)", .nm)))

  .nm <- .linCmtExport(.linCmtTwoIv, .ivLin(), "nonmem", .ctl)
  expect_true("$SUBROUTINES ADVAN3 TRANS1" %in% .nm)
  expect_true(any(grepl("^  K12=", .nm)))
  expect_true(any(grepl("^  K21=", .nm)))
  # the model variable k is not NONMEM's K
  expect_true(any(grepl("^  RXR[0-9]+=2", .nm)))
  expect_true(any(grepl("A\\(1\\)/\\(V\\)", .nm)))

  # ODEs when asked for
  .nm <- .linCmtExport(.linCmtOne, .theoLin(), "nonmem",
                       nonmemControl(runCommand=NA, modelName="x", linCmt="ode"))
  expect_true(any(grepl("^\\$SUBROUTINES ADVAN13", .nm)))
  expect_true(any(grepl("DADT\\(2\\)", .nm)))
})

test_that("linCmt() with other ODEs is translated to ODEs", {
  skip_on_cran()
  .d <- .theoLin()
  .d$CMT <- ifelse(.d$EVID == 0, "cp", "depot")
  .d2 <- .d[.d$EVID == 0, ]
  .d2$CMT <- "ce"
  .d <- rbind(.d, .d2)
  .d <- .d[order(.d$ID, .d$TIME, -.d$EVID), ]
  .nm <- .linCmtExport(.linCmtEffect, .d, "nonmem",
                       nonmemControl(runCommand=NA, modelName="x"))
  expect_true(any(grepl("^\\$SUBROUTINES ADVAN13", .nm)))
  expect_true(any(grepl("DADT\\(3\\) = KE0", .nm)))
  .mlx <- .linCmtExport(.linCmtEffect, .d, "monolix",
                        monolixControl(runCommand=NA, modelName="x"))
  expect_true(any(grepl("ddt_ce", .mlx)))
  expect_false(any(grepl("pkmodel", .mlx)))
})

test_that("linCmt() is translated to Monolix's pkmodel()", {
  skip_on_cran()
  skip_if_not(.hasLinCmtMicro())
  .ctl <- monolixControl(runCommand=NA, modelName="x")
  .mlx <- .linCmtExport(.linCmtOne, .theoLin(), "monolix", .ctl)
  expect_true(any(grepl("rx_cc = pkmodel\\(V=rx_v, k=rx_k, ka=rx_ka, Tlag=rx_tlag\\)", .mlx)))
  expect_false(any(grepl("ddt_|^PK:", .mlx)))
  expect_true(any(grepl("^ *cp = rx_cc", .mlx)))

  .mlx <- .linCmtExport(.linCmtTwoIv, .ivLin(), "monolix", .ctl)
  expect_true(any(grepl("pkmodel\\(V=rx_v, k=rx_k, k12=rx_k12, k21=rx_k21\\)", .mlx)))

  .mlx <- .linCmtExport(.linCmtOne, .theoLin(), "monolix",
                        monolixControl(runCommand=NA, modelName="x", linCmt="ode"))
  expect_true(any(grepl("ddt_central", .mlx)))
  expect_false(any(grepl("pkmodel", .mlx)))
})

test_that("a direct linCmt() translation asks for nlmixr2() or linToOde()", {
  .ui <- rxode2::rxode2(.linCmtOne)
  expect_error(.ui$nonmemModel, "linToOde")
  expect_error(.ui$monolixModel, "linToOde")
})

test_that("NONMEM writes modeled rate/duration and keeps the dosing items", {
  skip_on_cran()
  .iv <- function() {
    ini({
      tcl <- log(4)
      tv <- log(40)
      tdur <- log(2)
      eta.cl ~ 0.1
      eta.v ~ 0.1
      prop.sd <- 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      dur(central) <- exp(tdur)
      d/dt(central) <- -cl / v * central
      cp <- central / v
      cp ~ prop(prop.sd)
    })
  }
  .d <- .ivLin()
  .d$RATE <- ifelse(.d$EVID == 1, -2, 0)
  .nm <- .linCmtExport(.iv, .d, "nonmem", nonmemControl(runCommand=NA, modelName="x"))
  expect_true(any(grepl("^  D1=", .nm)))
  expect_false(any(grepl("DUR1", .nm)))
  expect_true(any(grepl("^\\$INPUT.* RATE", .nm)))

  .iv <- function() {
    ini({
      tcl <- log(4)
      tv <- log(40)
      trate <- log(250)
      eta.cl ~ 0.1
      eta.v ~ 0.1
      prop.sd <- 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      rate(central) <- exp(trate)
      d/dt(central) <- -cl / v * central
      cp <- central / v
      cp ~ prop(prop.sd)
    })
  }
  .d$RATE <- ifelse(.d$EVID == 1, -1, 0)
  .nm <- .linCmtExport(.iv, .d, "nonmem", nonmemControl(runCommand=NA, modelName="x"))
  expect_true(any(grepl("^  R1=", .nm)))

  # an infusion given by rate in the data
  .d <- .ivLin()
  .d$RATE <- ifelse(.d$EVID == 1, 100, 0)
  .nm <- .linCmtExport(.linCmtTwoIv, .d, "nonmem", nonmemControl(runCommand=NA, modelName="x"))
  expect_true(any(grepl("^\\$INPUT.* RATE", .nm)))
  # steady state
  .d <- .theoLin()
  .d$SS <- ifelse(.d$EVID == 1, 1, 0)
  .d$II <- ifelse(.d$EVID == 1, 24, 0)
  .nm <- .linCmtExport(.linCmtOne, .d, "nonmem", nonmemControl(runCommand=NA, modelName="x"))
  expect_true(any(grepl("^\\$INPUT.* II .* SS", .nm)))
})

test_that("erf() (from probitInv()) translates", {
  .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.linCmtOne))
  expect_equal(.rxToNonmem(quote(erf(x)), .ui),
               "(2*PHI(1.414213562373095145475*(X))-1)")
  expect_equal(.rxToMonolix(quote(erf(x)), .ui),
               "(2*normcdf(1.414213562373095145475*(x))-1)")
})

test_that("Monolix macros and distributions separate their arguments", {
  .adm <- data.frame(adm=2L, cmt=1L, type=factor("empty"), f=NA, dur=NA, lag=NA, rate=NA)
  expect_equal(.monolixGetPkMacrosForI(1, .adm, "central"),
               "empty(adm=2, target=central)")
})

test_that("rxode2's normalized powers translate", {
  .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.linCmtOne))
  expect_equal(.rxToNonmem(quote(Rx_pow_di(x, 2)), .ui), "(X)**(2)")
  expect_equal(.rxToMonolix(quote(Rx_pow_di(x, 2)), .ui), "(x)^(2)")
})

test_that("Monolix properties only apply to dosed compartments", {
  .ui <- rxode2::rxUiDecompress(rxode2::rxode2(.linCmtOne))
  rxode2::rxAssignControlValue(.ui, ".adm",
                               data.frame(adm=1L, cmt=2L, type=factor("bolus"),
                                          f=NA_character_, dur=NA_character_,
                                          lag=NA_character_, rate=NA_character_))
  # depot (cmt 1) has no doses: nothing changes
  .monolixSetAdm(.ui, "depot", "0.5", "f")
  expect_equal(nrow(.monolixGetAdm(.ui)), 1L)
  expect_true(is.na(.monolixGetAdm(.ui)$f))
  .monolixSetAdm(.ui, "central", "0.5", "f")
  expect_equal(.monolixGetAdm(.ui)$f, "0.5")
})
