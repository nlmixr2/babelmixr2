test_that("warfarin NONMEM reading", {

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

  if (file.exists("pk.turnover.emax3.zip")) {
    .path <- normalizePath("pk.turnover.emax3.zip")
    withr::with_tempdir({
      unzip(.path)
      # This has rounding errors
      expect_error(nlmixr(pk.turnover.emax3, nlmixr2data::warfarin, "nonmem",
                                   nonmemControl(readRounding=FALSE)))

      # Can still load the model to get information (possibly pipe) and create a new model
      f <- nlmixr(pk.turnover.emax3, nlmixr2data::warfarin, "nonmem",
                  nonmemControl(readRounding=TRUE))

      expect_true(inherits(f, "nlmixr2FitData"))

      # Will still error if you try to read this with readRounding=FALSE
      expect_error(nlmixr(pk.turnover.emax3, nlmixr2data::warfarin, "nonmem",
                          nonmemControl(readRounding=FALSE)))

      # Note this shouldn't have a covariance step so you can add it (at least a nlmixr2 covariance step)
      getVarCov(f)

      # nlmixr2 is more generous in what constitutes a covariance
      # step, in this case it is |r|,|s| which should be regarded with
      # caution but can give some clues on why this is not working in
      # NONMEM.

      # Here you can see the shrinkage is high for temax tktr and tka,
      # so they could be dropped with a model in nonmem that is more
      # likely to converge in NONMEM, starting from the model

      # In addition to dropping the problematic parameters, this will
      # restart the fit at the final initial estimates

      f %>% model(ktr <- exp(tktr)) %>%
        model(ka <- exp(tka)) %>%
        model(kout <- exp(tkout + eta.kout)) %>%
        nlmixr(data=nlmixr2data::warfarin, est="nonmem", control=nonmemControl(readRounding=FALSE))

    })
  }
})

test_that("pheno NONMEM reading", {

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

  if (file.exists("pheno.zip")) {
    .path <- normalizePath("pheno.zip")
    withr::with_tempdir({
      unzip(.path)
      f <- nlmixr2::nlmixr(pheno, nlmixr2data::pheno_sd, "nonmem")
      expect_true(inherits(f, "nlmixr2FitData"))
    })
  }

})

test_that("wbc NONMEM reading", {

  wbc <- function() {
    ini({
      ## Note that the UI can take expressions
      ## Also note that these initial estimates should be provided on the log-scale
      log_CIRC0 <- log(7.21)
      log_MTT <- log(124)
      log_SLOPU <- log(28.9)
      log_GAMMA <- log(0.239)
      ## Initial estimates should be high for SAEM ETAs
      eta.CIRC0  ~ .1
      eta.MTT  ~ .03
      eta.SLOPU ~ .2
      ##  Also true for additive error (also ignored in SAEM)
      prop.err <- 10
    })
    model({
      CIRC0 =  exp(log_CIRC0 + eta.CIRC0)
      MTT =  exp(log_MTT + eta.MTT)
      SLOPU =  exp(log_SLOPU + eta.SLOPU)
      GAMMA = exp(log_GAMMA)

      # PK parameters from input dataset
      CL = CLI;
      V1 = V1I;
      V2 = V2I;
      Q = 204;

      CONC = A_centr/V1;

      # PD parameters
      NN = 3;
      KTR = (NN + 1)/MTT;
      EDRUG = 1 - SLOPU * CONC;
      FDBK = (CIRC0 / A_circ)^GAMMA;

      CIRC = A_circ;

      A_prol(0) = CIRC0;
      A_tr1(0) = CIRC0;
      A_tr2(0) = CIRC0;
      A_tr3(0) = CIRC0;
      A_circ(0) = CIRC0;

      d/dt(A_centr) = A_periph * Q/V2 - A_centr * (CL/V1 + Q/V1);
      d/dt(A_periph) = A_centr * Q/V1 - A_periph * Q/V2;
      d/dt(A_prol) = KTR * A_prol * EDRUG * FDBK - KTR * A_prol;
      d/dt(A_tr1) = KTR * A_prol - KTR * A_tr1;
      d/dt(A_tr2) = KTR * A_tr1 - KTR * A_tr2;
      d/dt(A_tr3) = KTR * A_tr2 - KTR * A_tr3;
      d/dt(A_circ) = KTR * A_tr3 - KTR * A_circ;

      CIRC ~ prop(prop.err)
    })
  }

  .path <- normalizePath("wbc.zip")
  withr::with_tempdir({
      unzip(.path)
      f <- nlmixr2(wbc, nlmixr2data::wbcSim, "nonmem")
  })


})
