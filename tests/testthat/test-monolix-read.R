test_that("test monolix reading for 2019, 2020, and 2021", {

  pk.turnover.emax3 <- function() {
    ini({
      tktr <- log(1)
      tka <- log(1)
      tcl <- log(0.1)
      tv <- log(10)
      ##
      eta.ktr ~ 1
      eta.ka ~ 1
      eta.cl ~ 2
      eta.v ~ 1
      prop.err <- 0.1
      pkadd.err <- 0.1
      ##
      temax <- logit(0.8)
      tec50 <- log(0.5)
      tkout <- log(0.05)
      te0 <- log(100)
      ##
      eta.emax ~ .5
      eta.ec50  ~ .5
      eta.kout ~ .5
      eta.e0 ~ .5
      ##
      pdadd.err <- 10
    })
    model({
      ktr <- exp(tktr + eta.ktr)
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      emax = expit(temax+eta.emax)
      ec50 =  exp(tec50 + eta.ec50)
      kout = exp(tkout + eta.kout)
      e0 = exp(te0 + eta.e0)
      ##
      DCP = center/v
      PD=1-emax*DCP/(ec50+DCP)
      ##
      effect(0) = e0
      kin = e0*kout
      ##
      d/dt(depot) = -ktr * depot
      d/dt(gut) =  ktr * depot -ka * gut
      d/dt(center) =  ka * gut - cl / v * center
      d/dt(effect) = kin*PD -kout*effect
      ##
      cp = center / v
      cp ~ prop(prop.err) + add(pkadd.err)
      effect ~ add(pdadd.err) | pca
    })
  }

  if (file.exists("pk.turnover.emax3-2019.zip")) {
    .path <- normalizePath("pk.turnover.emax3-2019.zip")
    withr::with_tempdir({
      unzip(.path)
      f <- nlmixr2::nlmixr(pk.turnover.emax3, nlmixr2data::warfarin, "monolix")
      expect_true(inherits(f, "nlmixr2FitData"))
    })
  }

  if (file.exists("pk.turnover.emax3-2020.zip")) {
    .path <- normalizePath("pk.turnover.emax3-2020.zip")
    withr::with_tempdir({
      unzip(.path)
      f <- nlmixr2::nlmixr(pk.turnover.emax3, nlmixr2data::warfarin, "monolix")
      expect_true(inherits(f, "nlmixr2FitData"))
    })
  }

  if (file.exists("pk.turnover.emax3-2021.zip")) {
    .path <- normalizePath("pk.turnover.emax3-2021.zip")
    withr::with_tempdir({
      unzip(.path)
      f <- nlmixr2::nlmixr(pk.turnover.emax3, nlmixr2data::warfarin, "monolix")
      expect_true(inherits(f, "nlmixr2FitData"))
    })
  }

})


test_that("test more nlmixr2/monolix features", {

  pk.turnover.emax4 <- function() {
    ini({
      tktr <- log(1)
      tka <- log(1)
      tcl <- log(0.1)
      tv <- log(10)
      ##
      eta.ktr ~ 1
      eta.ka ~ 1
      eta.cl ~ 2
      eta.v ~ 1
      prop.err <- 0.1
      pkadd.err <- 0.1
      ##
      temax <- logit(0.8)
      tec50 <- log(0.5)
      tkout <- log(0.05)
      te0 <- log(100)
      ##
      eta.emax + eta.ec50 ~ c(.5,
                              0.05, .5)
      eta.kout ~ .5
      eta.e0 ~ .5
      ##
      pdadd.err <- 10
    })
    model({
      ktr <- exp(tktr + eta.ktr)
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      emax = expit(temax+eta.emax)
      ec50 =  exp(tec50 + eta.ec50)
      kout = exp(tkout + eta.kout)
      e0 = exp(te0 + eta.e0)
      ##
      DCP = center/v
      PD=1-emax*DCP/(ec50+DCP)
      ##
      effect(0) = e0
      kin = e0*kout
      ##
      d/dt(depot) = -ktr * depot
      d/dt(gut) =  ktr * depot -ka * gut
      d/dt(center) =  ka * gut - cl / v * center
      d/dt(effect) = kin*PD -kout*effect
      ##
      cp = center / v
      cp ~ prop(prop.err) + add(pkadd.err)
      effect ~ add(pdadd.err) | pca
    })
  }


  if (file.exists("pk.turnover.emax4-2021.zip")) {
    .path <- normalizePath("pk.turnover.emax4-2021.zip")
    withr::with_tempdir({
      unzip(.path)
      f <- nlmixr2::nlmixr(pk.turnover.emax4, nlmixr2data::warfarin, "monolix")
      expect_true(inherits(f, "nlmixr2FitData"))
    })
  }
})


test_that("test pheno", {

  pheno <- function() {
    ini({
      tcl <- log(0.008) # typical value of clearance
      tv <-  log(0.6)   # typical value of volume
      ## var(eta.cl)
      eta.cl + eta.v ~ c(1,
                         0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
      # interindividual variability on clearance and volume
      add.err <- 0.1    # residual variability
    })
    model({
      cl <- exp(tcl + eta.cl) # individual value of clearance
      v <- exp(tv + eta.v)    # individual value of volume
      ke <- cl / v            # elimination rate constant
      d/dt(A1) = - ke * A1    # model differential equation
      cp = A1 / v             # concentration in plasma
      cp ~ add(add.err)       # define error model
    })
  }


  if (file.exists("pheno-2021.zip")) {
    .path <- normalizePath("pheno-2021.zip")
    withr::with_tempdir({
      unzip(.path)
      f <- nlmixr2::nlmixr(pheno, nlmixr2data::pheno_sd, "monolix")
      expect_true(inherits(f, "nlmixr2FitData"))
    })
  }
})
