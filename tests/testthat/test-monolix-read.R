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

test_that("test monolix reading for 2019, 2020, and 2021", {
  
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


  skip_if_not(file.exists("pk.turnover.emax4-2021.zip"))
  .path <- normalizePath("pk.turnover.emax4-2021.zip")
  withr::with_tempdir({
    unzip(.path)
    f <- nlmixr2::nlmixr(pk.turnover.emax4, nlmixr2data::warfarin, "monolix")
    expect_true(inherits(f, "nlmixr2FitData"))
  })
  
})

test_that("test Monolix pheno", {

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

  skip_if_not(file.exists("pheno-2021.zip"))
  .path <- normalizePath("pheno-2021.zip")
  withr::with_tempdir({
    unzip(.path)
    f <- nlmixr2::nlmixr(pheno, nlmixr2data::pheno_sd, "monolix")
    expect_true(inherits(f, "nlmixr2FitData"))
  })
  
})

test_that("pbpk mavoglurant", {

  pbpk <- function(){
    ini({
      ##theta=exp(c(1.1, .3, 2, 7.6, .003, .3))
      lKbBR = 1.1
      lKbMU = 0.3
      lKbAD = 2

      lCLint = 7.6
      lKbBO = 0.03
      lKbRB = 0.3
      eta.LClint ~ 4
      lnorm.sd <- 1
    })
    model({
      KbBR = exp(lKbBR)
      KbMU = exp(lKbMU)
      KbAD = exp(lKbAD)
      CLint= exp(lCLint + eta.LClint)
      KbBO = exp(lKbBO)
      KbRB = exp(lKbRB)

      ## Regional blood flows
      CO  = (187.00*WT^0.81)*60/1000;         # Cardiac output (L/h) from White et al (1968)
      QHT = 4.0 *CO/100;
      QBR = 12.0*CO/100;
      QMU = 17.0*CO/100;
      QAD = 5.0 *CO/100;
      QSK = 5.0 *CO/100;
      QSP = 3.0 *CO/100;
      QPA = 1.0 *CO/100;
      QLI = 25.5*CO/100;
      QST = 1.0 *CO/100;
      QGU = 14.0*CO/100;
      QHA = QLI - (QSP + QPA + QST + QGU); # Hepatic artery blood flow
      QBO = 5.0 *CO/100;
      QKI = 19.0*CO/100;
      QRB = CO - (QHT + QBR + QMU + QAD + QSK + QLI + QBO + QKI);
      QLU = QHT + QBR + QMU + QAD + QSK + QLI + QBO + QKI + QRB;

      ## Organs' volumes = organs' weights / organs' density
      VLU = (0.76 *WT/100)/1.051;
      VHT = (0.47 *WT/100)/1.030;
      VBR = (2.00 *WT/100)/1.036;
      VMU = (40.00*WT/100)/1.041;
      VAD = (21.42*WT/100)/0.916;
      VSK = (3.71 *WT/100)/1.116;
      VSP = (0.26 *WT/100)/1.054;
      VPA = (0.14 *WT/100)/1.045;
      VLI = (2.57 *WT/100)/1.040;
      VST = (0.21 *WT/100)/1.050;
      VGU = (1.44 *WT/100)/1.043;
      VBO = (14.29*WT/100)/1.990;
      VKI = (0.44 *WT/100)/1.050;
      VAB = (2.81 *WT/100)/1.040;
      VVB = (5.62 *WT/100)/1.040;
      VRB = (3.86 *WT/100)/1.040;

      ## Fixed parameters
      BP = 0.61;      # Blood:plasma partition coefficient
      fup = 0.028;    # Fraction unbound in plasma
      fub = fup/BP;   # Fraction unbound in blood

      KbLU = exp(0.8334);
      KbHT = exp(1.1205);
      KbSK = exp(-.5238);
      KbSP = exp(0.3224);
      KbPA = exp(0.3224);
      KbLI = exp(1.7604);
      KbST = exp(0.3224);
      KbGU = exp(1.2026);
      KbKI = exp(1.3171);


      ##-----------------------------------------
      S15 = VVB*BP/1000;
      C15 = Venous_Blood/S15

      ##-----------------------------------------
      d/dt(Lungs) = QLU*(Venous_Blood/VVB - Lungs/KbLU/VLU);
      d/dt(Heart) = QHT*(Arterial_Blood/VAB - Heart/KbHT/VHT);
      d/dt(Brain) = QBR*(Arterial_Blood/VAB - Brain/KbBR/VBR);
      d/dt(Muscles) = QMU*(Arterial_Blood/VAB - Muscles/KbMU/VMU);
      d/dt(Adipose) = QAD*(Arterial_Blood/VAB - Adipose/KbAD/VAD);
      d/dt(Skin) = QSK*(Arterial_Blood/VAB - Skin/KbSK/VSK);
      d/dt(Spleen) = QSP*(Arterial_Blood/VAB - Spleen/KbSP/VSP);
      d/dt(Pancreas) = QPA*(Arterial_Blood/VAB - Pancreas/KbPA/VPA);
      d/dt(Liver) = QHA*Arterial_Blood/VAB + QSP*Spleen/KbSP/VSP + QPA*Pancreas/KbPA/VPA + QST*Stomach/KbST/VST + QGU*Gut/KbGU/VGU - CLint*fub*Liver/KbLI/VLI - QLI*Liver/KbLI/VLI;
      d/dt(Stomach) = QST*(Arterial_Blood/VAB - Stomach/KbST/VST);
      d/dt(Gut) = QGU*(Arterial_Blood/VAB - Gut/KbGU/VGU);
      d/dt(Bones) = QBO*(Arterial_Blood/VAB - Bones/KbBO/VBO);
      d/dt(Kidneys) = QKI*(Arterial_Blood/VAB - Kidneys/KbKI/VKI);
      d/dt(Arterial_Blood) = QLU*(Lungs/KbLU/VLU - Arterial_Blood/VAB);
      d/dt(Venous_Blood) = QHT*Heart/KbHT/VHT + QBR*Brain/KbBR/VBR + QMU*Muscles/KbMU/VMU + QAD*Adipose/KbAD/VAD + QSK*Skin/KbSK/VSK + QLI*Liver/KbLI/VLI + QBO*Bones/KbBO/VBO + QKI*Kidneys/KbKI/VKI + QRB*Rest_of_Body/KbRB/VRB - QLU*Venous_Blood/VVB;
      d/dt(Rest_of_Body) = QRB*(Arterial_Blood/VAB - Rest_of_Body/KbRB/VRB);

      C15 ~ lnorm(lnorm.sd)
    })
  }

  expect_error(bblDatToMonolix(pbpk, nlmixr2data::mavoglurant), NA)

  skip_if_not(file.exists("pbpk-2021.zip"))
  .path <- normalizePath("pbpk-2021.zip")
  withr::with_tempdir({
    unzip(.path)
    f <- nlmixr2::nlmixr(pbpk, nlmixr2data::mavoglurant, "monolix")
    expect_true(inherits(f, "nlmixr2FitData"))
  })
  
})

test_that("nimo test", {
  
  nimo <- function() {
    ini({
      ## Note that the UI can take expressions
      ## Also note that these initial estimates should be provided on the log-scale
      tcl <- log(0.001)
      tv1 <- log(1.45)
      tQ <- log(0.004)
      tv2 <- log(44)
      tkss <- log(12)
      tkint <- log(0.3)
      tksyn <- log(1)
      tkdeg <- log(7)
      ## Initial estimates should be high for SAEM ETAs
      eta.cl  ~ 2
      eta.v1  ~ 2
      eta.kss ~ 2
      ##  Also true for additive error (also ignored in SAEM)
      lnorm.sd <- 10
    })
    model({
      cl <- exp(tcl + eta.cl)
      v1 <- exp(tv1 + eta.v1)
      Q  <- exp(tQ)
      v2 <- exp(tv2)
      kss <- exp(tkss + eta.kss)
      kint <- exp(tkint)
      ksyn <- exp(tksyn)
      kdeg <- exp(tkdeg)

      k <- cl/v1
      k12 <- Q/v1
      k21 <- Q/v2

      eff(0) <- ksyn/kdeg ##initializing compartment

      ## Concentration is calculated
      conc = 0.5*(central/v1-eff-kss)+0.5*sqrt((central/v1-eff-kss)**2+4*kss*central/v1)

      d/dt(central)  = -(k+k12)*conc*v1+k21*peripheral-kint*eff*conc*v1/(kss+conc)
      d/dt(peripheral) = k12*conc*v1-k21*peripheral  ##Free Drug second compartment amount
      d/dt(eff) = ksyn - kdeg*eff - (kint-kdeg)*conc*eff/(kss+conc)

      conc ~ lnorm(lnorm.sd)

    })
  }

  tmp <- nlmixr2data::nimoData

  tmp$DV <-exp(tmp$DV)

  skip_if_not(file.exists("nimo-2021.zip"))
  .path <- normalizePath("nimo-2021.zip")
  withr::with_tempdir({
    unzip(.path)
    f <- suppressWarnings(nlmixr2::nlmixr2(nimo, tmp, "monolix"))
    expect_true(inherits(f, "nlmixr2FitData"))
  })
  
})

# WBC Model tests ####

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

test_that("Monolix wbc", {
  skip_if_not(file.exists("wbc-2021.zip"))
  skip_if_not(file.exists("wbc-charts-2021.zip"))
  .path <- normalizePath("wbc-2021.zip")
  .pathCharts <- normalizePath("wbc-charts-2021.zip")
  withr::with_tempdir({
    unzip(.path)
    
    f <- suppressWarnings(nlmixr2::nlmixr2(wbc, nlmixr2data::wbcSim, "monolix"))
    expect_true(inherits(f, "nlmixr2FitData"))
    expect_equal(f$env$parHist, NULL)
    
    unzip(.pathCharts)
    
    f <- suppressWarnings(nlmixr2::nlmixr2(wbc, nlmixr2data::wbcSim, "monolix"))
    expect_true(inherits(f, "nlmixr2FitData"))
    expect_true(inherits(f$env$parHist, "data.frame"))
    
    f <- suppressWarnings(nlmixr2::nlmixr2(wbc, nlmixr2data::wbcSim, "monolix"))
    expect_true(inherits(f, "nlmixr2FitData"))
    expect_true(inherits(f$env$parHist, "data.frame"))
    
  })
  
  withr::with_tempdir({
    unzip(.path)
    unzip(.pathCharts)
    f <- suppressWarnings(nlmixr2::nlmixr2(wbc, nlmixr2data::wbcSim, "monolix"))
    expect_true(inherits(f, "nlmixr2FitData"))
    expect_true(inherits(f$env$parHist, "data.frame"))
  })
})

test_that("Monolix wbc test 2", {
  skip_if_not(file.exists("x-2021.zip"))
  .path <- normalizePath("x-2021.zip")
  withr::with_tempdir({
    unzip(.path)
    # with piping, the original model name is scrubbed, and is changed to x
    p1 <- wbc %>%
      ini(prop.err=15) %>%
      nlmixr2(., nlmixr2data::wbcSim, "monolix")
    
    expect_true(inherits(p1, "nlmixr2FitData"))
    
    p2 <-wbc %>%
      ini(prop.err=7) %>%
      nlmixr2(., nlmixr2data::wbcSim, "monolix")
    
    expect_true(inherits(p2, "nlmixr2FitData"))
  })
})
