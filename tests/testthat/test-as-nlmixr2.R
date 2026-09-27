.nonmem2rx <- function(...) suppressWarnings(suppressMessages(nonmem2rx::nonmem2rx(...)))

.as.nonmem2rx <- function(...) suppressWarnings(suppressMessages(nonmem2rx::as.nonmem2rx(...)))

.as.nlmixr2 <- .as.nlmixr <- function(...) suppressWarnings(suppressMessages(as.nlmixr(...)))

test_that("nlmixr2 translation from nonmem2x", {
  skip_on_cran()

  mod <- .nonmem2rx(system.file("mods/cpt/runODE032.ctl", package="nonmem2rx"),
                   determineError=FALSE, lst=".res", save=FALSE)

  mod2 <-function() {
    ini({
      lcl <- 1.37034036528946
      lvc <- 4.19814911033061
      lq <- 1.38003493562413
      lvp <- 3.87657341967489
      RSV <- c(0, 0.196446108190896, 1)
      eta.cl ~ 0.101251418415006
      eta.v ~ 0.0993872449483344
      eta.q ~ 0.101302674763154
      eta.v2 ~ 0.0730497519364148
    })
    model({
      cmt(CENTRAL)
      cmt(PERI)
      cl <- exp(lcl + eta.cl)
      v <- exp(lvc + eta.v)
      q <- exp(lq + eta.q)
      v2 <- exp(lvp + eta.v2)
      v1 <- v
      scale1 <- v
      k21 <- q/v2
      k12 <- q/v
      d/dt(CENTRAL) <- k21 * PERI - k12 * CENTRAL - cl * CENTRAL/v1
      d/dt(PERI) <- -k21 * PERI + k12 * CENTRAL
      f <- CENTRAL/scale1
      f ~ prop(RSV)
    })
  }

  new <- .as.nonmem2rx(mod2, mod)

  expect_true(inherits(.as.nlmixr(new), "nlmixr2FitData"))

  mod <- .nonmem2rx(system.file("mods/cpt/runODE032.ctl", package="nonmem2rx"),
                   determineError=TRUE, lst=".res", save=FALSE)

  new <- .as.nonmem2rx(mod2, mod)

  fit <- .as.nlmixr(new)
  expect_true(inherits(fit, "nlmixr2FitData"))
  expect_true(any(names(fit$time) == "NONMEM"))

  # tableControl(cwres=TRUE) adds nlmixr2's FOCEi objective, but the
  # imported NONMEM objective stays in use (#94)
  fit <- .as.nlmixr(new, table = tableControl(cwres = TRUE))
  expect_true("CWRES" %in% names(fit))
  expect_setequal(row.names(fit$objDf), c("nonmem2rx", "FOCEi"))
  expect_equal(fit$ofvType, "nonmem2rx")
  expect_equal(fit$objective, fit$objDf["nonmem2rx", "OBJF"])
  expect_equal(AIC(fit), fit$objDf["nonmem2rx", "AIC"])

  # a different model right after the cwres=TRUE import must not start
  # from the etas of the FOCEi objective's nlmixr2() fit (#94)
  rx <- .nonmem2rx(system.file("mods/err/run006.lst", package="nonmem2rx"))
  fit <- .as.nlmixr(rx)
  expect_true(inherits(fit, "nlmixr2FitData"))
  expect_true(any(names(fit$time) == "NONMEM"))

})

test_that("as.nlmixr2 gives a clear error for untranslated eps/err (#95)", {
  skip_on_cran()

  mod <- .nonmem2rx(
    system.file("mods/cpt/runODE032.ctl", package = "nonmem2rx"),
    determineError = FALSE,
    lst = ".res",
    save = FALSE
  )

  expect_error(as.nlmixr2(mod), "'eps1'")

  mod2 <- function() {
    ini({
      lcl <- 1.37034036528946
      lvc <- 4.19814911033061
      lq <- 1.38003493562413
      lvp <- 3.87657341967489
      RSV <- c(0, 0.196446108190896, 1)
      eta.cl ~ 0.101251418415006
      eta.v ~ 0.0993872449483344
      eta.q ~ 0.101302674763154
      eta.v2 ~ 0.0730497519364148
    })
    model({
      cmt(CENTRAL)
      cmt(PERI)
      cl <- exp(lcl + eta.cl)
      v <- exp(lvc + eta.v)
      q <- exp(lq + eta.q)
      v2 <- exp(lvp + eta.v2)
      k21 <- q / v2
      k12 <- q / v
      d / dt(CENTRAL) <- k21 * PERI - k12 * CENTRAL - cl * CENTRAL / v
      d / dt(PERI) <- -k21 * PERI + k12 * CENTRAL
      f <- CENTRAL / v
      y <- f + f * eps1
      f ~ prop(RSV)
    })
  }

  new <- .as.nonmem2rx(mod2, mod)

  expect_error(
    as.nlmixr2(new),
    "untranslated NONMEM residual variable\\(s\\): 'eps1'"
  )

  expect_error(
    .nonmem2rxAssertNoEps(list(
      allCovs = c("WT", "err1", "eps2", "eps1x"),
      nonmemData = data.frame(WT = 1)
    )),
    "'err1', 'eps2'\\n"
  )
  expect_error(
    .nonmem2rxAssertNoEps(list(
      allCovs = c("WT", "eps1"),
      nonmemData = data.frame(WT = 1, eps1 = 0)
    )),
    NA
  )
})

.monolix2rx <- function(...) suppressWarnings(suppressMessages(monolix2rx::monolix2rx(...)))

test_that("nlmixr2 translation from monolix2rx", {
  skip_on_cran()

  pkgTheo <- system.file("theo/theophylline_project.mlxtran", package="monolix2rx")

  mod <- .monolix2rx(pkgTheo)

  fit <- .as.nlmixr2(mod)

  expect_true(inherits(fit, "nlmixr2FitData"))

  # the imported Monolix objective stays in use with cwres=TRUE (#94)
  fit <- .as.nlmixr2(mod, table = tableControl(cwres = TRUE))
  expect_true("CWRES" %in% names(fit))
  expect_true("FOCEi" %in% row.names(fit$objDf))
  expect_equal(nrow(fit$objDf), 2L)
  expect_false(fit$ofvType == "FOCEi")
  expect_equal(fit$objective, fit$objDf[fit$ofvType, "OBJF"])
  expect_equal(AIC(fit), fit$objDf[fit$ofvType, "AIC"])
  expect_equal(BIC(fit), fit$objDf[fit$ofvType, "BIC"])
  expect_equal(
    as.numeric(logLik(fit)),
    fit$objDf[fit$ofvType, "Log-likelihood"]
  )

  # a different model right after the Monolix cwres=TRUE import must not
  # start from the etas of the FOCEi objective's nlmixr2() fit (#94)
  rx <- .nonmem2rx(system.file("mods/err/run006.lst", package = "nonmem2rx"))
  fit <- .as.nlmixr2(rx)
  expect_true(inherits(fit, "nlmixr2FitData"))
})

test_that(".importEtaMat() only uses etas that match the model (#94)", {
  .b <- loadNamespace("babelmixr2")
  .ui <- rxode2::rxode2(function() {
    ini({
      tcl <- 1
      tv <- 2
      eta.v ~ 0.1
      eta.cl ~ 0.1
      add.sd <- 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  })
  .obf <- data.frame(
    ID = 1:2,
    eta.cl = c(0.1, 0.2),
    eta.v = c(0.3, 0.4),
    OBJI = NA_real_
  )
  # columns follow the model's eta order, not the etaObf order
  expect_equal(
    .b$.importEtaMat(.ui, .obf, 2L),
    matrix(c(0.3, 0.4, 0.1, 0.2), 2, 2)
  )
  # a subject dropped from the processed data (for example one without
  # observations), missing etas or NA etas give zero etas
  expect_equal(.b$.importEtaMat(.ui, .obf, 1L), matrix(0, 1, 2))
  expect_equal(
    .b$.importEtaMat(.ui, .obf[, c("ID", "eta.cl", "OBJI")], 2L),
    matrix(0, 2, 2)
  )
  .na <- .obf
  .na$eta.v[1] <- NA_real_
  expect_equal(.b$.importEtaMat(.ui, .na, 2L), matrix(0, 2, 2))
  # no etas, no matrix
  expect_null(.b$.importEtaMat(
    rxode2::rxode2(function() {
      ini({
        tcl <- 1
        add.sd <- 0.1
      })
      model({
        cp <- exp(tcl)
        cp ~ add(add.sd)
      })
    }),
    .obf[, c("ID", "OBJI")],
    2L
  ))
})
