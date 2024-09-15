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

  rx <- .nonmem2rx(system.file("mods/err/run006.lst", package="nonmem2rx"))
  fit <- .as.nlmixr(rx)
  expect_true(inherits(fit, "nlmixr2FitData"))
  expect_true(any(names(fit$time) == "NONMEM"))

})

.monolix2rx <- function(...) suppressWarnings(suppressMessages(monolix2rx::monolix2rx(...)))

test_that("nlmixr2 translation from monolix2rx", {
  skip_on_cran()

  pkgTheo <- system.file("theo/theophylline_project.mlxtran", package="monolix2rx")

  mod <- .monolix2rx(pkgTheo)

  fit <- .as.nlmixr2(mod)

  expect_true(inherits(fit, "nlmixr2FitData"))

})
