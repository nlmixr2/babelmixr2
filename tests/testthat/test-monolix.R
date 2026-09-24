.b <- loadNamespace("babelmixr2")

if (FALSE) {
  test_that("pure mu refrence parsing", {

    # this was moved to rxode2

    expect_equal(.b$getPureMuRef(quote(cl <- tcl),
                                 muRefCurEval=data.frame(parameter="tcl", curEval="",
                                                         low=NA_character_, hi=NA_character_)),
                 c(tcl="cl"))

    expect_equal(.b$.getPureMuRef(quote(cl <- tcl),
                                  muRefCurEval=data.frame(parameter="tcl", curEval="exp",
                                                          low=NA_character_, hi=NA_character_)), NULL)

    expect_equal(.b$.getPureMuRef(quote(cl <- exp(tcl)),
                                  muRefCurEval=data.frame(parameter="tcl", curEval="exp",
                                                          low=NA_character_, hi=NA_character_)),
                 c(tcl="cl"))

    expect_equal(.b$.getPureMuRef(quote(cl <- exp(tcl)),
                                  muRefCurEval=data.frame(parameter="tcl", curEval="",
                                                          low=NA_character_, hi=NA_character_)),
                 NULL)


    expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl, 0, 1)),
                                  muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                          low=NA_character_, hi=NA_character_)),
                 c(tcl="cl"))


    expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl, 0, 2)),
                                  muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                          low=NA_character_, hi=NA_character_)),
                 NULL)


    expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl, 0, 2)),
                                  muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                          low=0, hi=2)),
                 c(tcl="cl"))


    expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl, 0.5, 1)),
                                  muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                          low=NA_character_, hi=NA_character_)),
                 NULL)

    expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl, 0.5, 1)),
                                  muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                          low=0.5, hi=NA_character_)),
                 c(tcl="cl"))


    expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl, 0.5)),
                                  muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                          low=0.5, hi=NA_character_)),
                 c(tcl="cl"))


    expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl, 0.5)),
                                  muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                          low=0.4, hi=NA_character_)),
                 NULL)

    expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl)),
                                  muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                          low=0, hi=1)),
                 c(tcl="cl"))

    expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl)),
                                  muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                          low=NA_real_, hi=1)),
                 c(tcl="cl"))

    expect_equal(.b$.getPureMuRef(quote(cl <- expit(tcl)),
                                  muRefCurEval=data.frame(parameter="tcl", curEval="expit",
                                                          low=0, hi=NA_real_)),
                 c(tcl="cl"))

    expect_equal(.b$.getPureMuRef(quote(cl(0) <- tcl),
                                  muRefCurEval=data.frame(parameter="tcl", curEval="",
                                                          low=NA_real_, hi=NA_real_)),
                 NULL)
  })

}


test_that("individual distribution switch", {
  expect_equal(.b$.mlxTranCurEvalToDistribution("exp"),
               "distribution=logNormal")
  expect_equal(.b$.mlxTranCurEvalToDistribution("expit"),
               "distribution=logitNormal")
  expect_equal(.b$.mlxTranCurEvalToDistribution("probitInv"),
               "distribution=probitNormal")
  expect_equal(.b$.mlxTranCurEvalToDistribution(""),
               "distribution=normal")
  expect_error(.b$.mlxTranCurEvalToDistribution("log"))
})

test_that("can determine if parameter is population only", {
  .df <- data.frame(theta = c("tktr", "tka", "tcl", "tv", "tkout", "te0", "tdepot"),
                    eta = c("eta.ktr", "eta.ka", "eta.cl", "eta.v", "eta.kout", "eta.e0", "eta.depot"),
                    level="id")
  expect_true(.b$.mlxTranIsPopOnly("temax", .df))
  expect_false(.b$.mlxTranIsPopOnly("tka", .df))
})

test_that("get variability component", {

  .df <- data.frame(theta = c("tktr", "tka", "tcl", "tv", "tkout", "te0", "tdepot"),
                    eta = c("eta.ktr", "eta.ka", "eta.cl", "eta.v", "eta.kout", "eta.e0", "eta.depot"),
                    level="id")
  expect_equal(.b$.mlxTranGetVaraibility("emax", "temax", .df),
               "no-variability")

  expect_equal(.b$.mlxTranGetVaraibility("ka", "tka", .df),
               "sd=omega_ka")

})

test_that("test datafile use", {

  one.cmt <- function() {
    ini({
      tka <- 0.45 ; label("Ka")
      tcl <- log(c(0, 2.7, 100)) ; label("Log Cl")
      tv <- 3.45; label("log V")
      cl.wt <- 0
      v.wt <- 0
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + wt * cl.wt)
      v <- exp(tv + eta.v) + wt2 ^ 2 * v.wt
      linCmt() ~ add(add.sd)
    })
  }

  ui <- rxode2::rxode2(one.cmt)

  expect_equal(.b$.monolixMapDataUse("ID", ui), "ID = {use=identifier}")
  expect_equal(.b$.monolixMapDataUse("TIME", ui), "TIME = {use=time}")
  expect_equal(.b$.monolixMapDataUse("EVID", ui), "EVID = {use=eventidentifier}")
  expect_equal(.b$.monolixMapDataUse("AMT", ui), "AMT = {use=amount}")
  expect_equal(.b$.monolixMapDataUse("II", ui), "II = {use=interdoseinterval}")
  expect_equal(.b$.monolixMapDataUse("DV", ui), "DV = {use=observation, name=rx_prd_rxLinCmt, type=continuous}")
  expect_equal(.b$.monolixMapDataUse("CENS", ui), "CENS = {use=censored}")
  expect_equal(.b$.monolixMapDataUse("LIMIT", ui), "LIMIT = {use=limit}")
  expect_equal(.b$.monolixMapDataUse("YTYPE", ui), "YTYPE = {use=observationtype}")
  expect_equal(.b$.monolixMapDataUse("ADM", ui), "ADM = {use=administration}")
  expect_equal(.b$.monolixMapDataUse("SS", ui), "SS = {use=steadystate}")
  expect_equal(.b$.monolixMapDataUse("wt2", ui), "wt2 = {use=regressor}")
  expect_equal(.b$.monolixMapDataUse("wt", ui), "wt = {use=covariate, type=continuous}")
  expect_equal(.b$.monolixMapDataUse("nlmixrRowNums", ui), "")

  one.cmt <- function() {
    ini({
      tka <- 0.45 ; label("Ka")
      tcl <- log(c(0, 2.7, 100)) ; label("Log Cl")
      tv <- 3.45; label("log V")
      cl.wt <- 0
      v.wt <- 0
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + WT * cl.wt)
      v <- exp(tv + eta.v)+ WT ^ 2 * v.wt
      linCmt() ~ add(add.sd)
    })
  }

  ui <- rxode2::rxode2(one.cmt)

  expect_equal(.monolixMapDataUse("wt2", ui), "")
  # This is only true with the new rxode2;
  expect_equal(.monolixMapDataUse("WT", ui), "WT = {use=regressor}")

})

test_that("monolix dsl", {

  one.cmt <- function() {
    ini({
      tka <- 0.45 ; label("Ka")
      tcl <- log(c(0, 2.7, 100)) ; label("Log Cl")
      tv <- 3.45; label("log V")
      cl.wt <- 0
      v.wt <- 0
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + WT * cl.wt)
      v <- exp(tv + eta.v)+ WT ^ 2 * v.wt
      linCmt() ~ add(add.sd)
    })
  }

  ui <- rxode2::rxode2(one.cmt)

  .rxToM <- function(x) {
    rxToMonolix(x, ui)
  }
  .ee <- function(x, y) {
    .x <- gsub(" +", " ", x)
    .x <- gsub("^ +", "", .x)
    .x <- gsub(" +$", "", .x)
    .x <- gsub(" *\n *", "\n", .x)

    .y <- gsub(" +", " ", y)
    .y <- gsub("^ +", "", .y)
    .y <- gsub(" +$", "", .y)
    .y <- gsub(" *\n *", "\n", .y)
    expect_equal(.x, .y)
  }

  .ee(.rxToM("sqrt(a)"), "sqrt(a)")
  .ee(.rxToM("max(a,b)"), "max(a,b)")
  expect_error(.rxToM("max(a,b,c)"))
  expect_error(.rxToM("max(a)"))
  .ee(.rxToM("sum(a,b,c,d)"), "((a)+(b)+(c)+(d))")
  .ee(.rxToM("prod(a,b,c,d)"), "((a)*(b)*(c)*(d))")
  .ee(.rxToM("a<-1+b"), "a = 1+b")
  .ee(.rxToM("a~1+b"), "a = 1+b")
  .ee(.rxToM("a=1+b"), "a = 1+b")
  .ee(.rxToM("expit(a)"), "1/(1+exp(-(a)))")
  .ee(.rxToM("expit(a,b)"), "(1.0-(b))*(1/(1+exp(-(a))))+(b)")
  .ee(.rxToM("expit(a,b,c)"), "((c)-(b))*(1/(1+exp(-(a))))+(b)")
  .ee(.rxToM("logit(a)"), "-log(1/(a)-1)")
  .ee(.rxToM("logit(a,b)"), "-log(1/(((a)-(b))/(1.0-(b)))-1)")
  .ee(.rxToM("logit(a,b,c)"), "-log(1/(((a)-(b))/((c)-(b)))-1)")
  .ee(.rxToM("probitInv(a)"), "normcdf(a)")
  .ee(.rxToM("probitInv(a,b)"), "(1.0-(b))*(normcdf(a))+(b)")
  .ee(.rxToM("probitInv(a,b,c)"), "((c)-(b))*(normcdf(a))+(b)")
  .ee(.rxToM("probit(a)"), "probit(a)")
  .ee(.rxToM("probit(a,b)"), "probit(((a)-(b))/(1.0-(b)))")
  .ee(.rxToM("probit(a,b,c)"), "probit(((a)-(b))/((c)-(b)))")
  .ee(.rxToM("d/dt(depot)=-depot*kel"), "ddt_depot = - depot*kel")
  .ee(.rxToM("depot(0)=50"), "depot_0 = 50")
  .ee(.rxToM("f(depot)=3"), ";f defined in PK section")
  .ee(.rxToM("a**b"), "a^b")
  .ee(.rxToM("if (a<=b){c=1} else if (a==4) {c=2} else {c=4}"), "if a<=b\n  c = 1\nelseif a==4\n  c = 2\nelse \n  c = 4\nend\n")
  .ee(.rxToM("if (a<=b){c=1} else if (a==4) {c=2} else if (a==30) {c=4} else {c=100}"), "if a<=b\n  c = 1\nelseif a==4\n  c = 2\nelseif a==30\n  c = 4\nelse \n  c = 100\nend\n")
  .ee(.rxToM("if (a<=b){c=1} else if (a==4) {c=2}"), "if a<=b\n  c = 1\nelseif a==4\n  c = 2\nend\n")
  .ee(.rxToM("if (a<=b){c=1}"), "if a<=b\n  c = 1\nend\n")
  .ee(.rxToM("time"), "t")
  expect_error(.rxToM("NA"))
  expect_error(.rxToM("newind"))
  .ee(.rxToM("log1pmx(a)"), "(log(1+a)-(a))")

  .ee(.rxToM("4.3"), "4.3")
  .ee(.rxToM("add.sd"), "add__sd")

})

test_that("monolix model creation without running", {
  withr::with_tempdir({

    one.cmt <- function() {
      ini({
        tka <- 0.45 ; label("Ka")
        tcl <- log(c(0, 2.7, 100)) ; label("Log Cl")
        tv <- 3.45; label("log V")
        cl.wt <- 0
        v.wt <- 0
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + WT * cl.wt)
        v <- exp(tv + eta.v)+ WT ^ 2 * v.wt
        d/dt(depot) <- -depot*ka
        d/dt(central) <- depot*ka - cl*central/v
        cp <-central/v
        cp ~ add(add.sd)
      })
    }

    files <- c("monolixTest-monolix.csv", "monolixTest-monolix.mlxtran",
               "monolixTest-monolix.md5", "monolixTest-monolix.txt")

    nlmixr2(one.cmt, nlmixr2data::theo_sd, "monolix",
            monolixControl(runCommand=NA, modelName="monolixTest"))

    lapply(files, function(f) { expect_true(file.exists(f)) })

    nlmixr2(one.cmt, nlmixr2data::theo_sd, "monolix",
            monolixControl(runCommand=NA, modelName="monolixTest"))
    lapply(files, function(f) {
      expect_true(file.exists(f))
      unlink(f)
    })

  })
})


test_that("monolix treatment of +var()", {

  one.cmt <- function() {
    ini({
      tka <- 0.45 ; label("Ka")
      tcl <- log(c(0, 2.7, 100)) ; label("Log Cl")
      tv <- 3.45; label("log V")
      cl.wt <- 0
      v.wt <- 0
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + WT * cl.wt)
      v <- exp(tv + eta.v)+ WT ^ 2 * v.wt
      d/dt(depot) <- -depot*ka
      d/dt(central) <- depot*ka - cl*central/v
      cp <-central/v
      cp ~ add(add.sd)
    })
  }

  f <- one.cmt()

  par <- strsplit(f$mlxtranParameter, "\n")[[1]]
  par <- par[grepl("add__sd", par)]

  expect_equal(par, "add__sd={value=0.7, method=MLE}")

  f2 <- try(f %>% model({cp~add(add.sd) + var()}), silent=TRUE)
  if (inherits(f2, "try-error")) {
    skip("rxode2 doesn't support + var()")
  } else {
    par <- strsplit(f2$mlxtranParameter, "\n")[[1]]
    par <- par[grepl("add__sd", par)]
    expect_equal(par, paste0("add__sd={value=", sqrt(.7), ", method=MLE}"))
  }

})

test_that("dotted mu-referenced parameters use one Monolix name (#220)", {

  f <- function() {
    ini({
      tka <- 0.45; tcl <- log(2.7); tv <- 3.45; cl.wt <- 0
      eta.ka + eta.cl ~ c(0.6, 0.01, 0.3)
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka.x <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + WT * cl.wt)
      v <- exp(tv + eta.v + WT * cl.wt)
      d/dt(depot) <- -depot*ka.x
      d/dt(central) <- depot*ka.x - cl*central/v
      cp <- central/v
      cp ~ add(add.sd)
    })
  }

  u <- rxode2::rxode2(f)

  expect_equal(.b$.monolixMuRef(u),
               c(tka="ka__x", tcl="rx__tcl", tv="rx__tv", cl.wt="rx__cl__wt"))

  mod <- u$monolixModel
  expect_true(grepl("input={ka__x,rx__tcl,rx__tv,rx__cl__wt,WT}", mod, fixed=TRUE))
  expect_true(grepl("cl = exp(rx__tcl+WT*rx__cl__wt)", mod, fixed=TRUE))
  expect_true(grepl("ddt_depot = - depot*ka__x", mod, fixed=TRUE))

  ind <- u$mlxtranModelIndividual
  expect_false(grepl(".", ind, fixed=TRUE))
  expect_true(grepl("ka__x = {distribution=logNormal, typical=ka__x_pop, sd=omega_ka__x}",
                    ind, fixed=TRUE))
  expect_true(grepl("rx__cl__wt = {distribution=normal, typical=rx__cl__wt_pop, no-variability}",
                    ind, fixed=TRUE))
  expect_true(grepl("r(rx__tcl, ka__x)=corr_rx__tcl_ka__x", ind, fixed=TRUE))

  par <- strsplit(u$mlxtranParameter, "\n")[[1]]
  expect_true(all(c("ka__x_pop", "rx__cl__wt_pop", "omega_ka__x",
                    "corr_rx__tcl_ka__x") %in% sub("=.*$", "", par)))

  # read-back: the readers look up the same names in Monolix's output
  pop <- data.frame(parameter=c("ka__x_pop", "rx__tcl_pop", "rx__tv_pop",
                                "rx__cl__wt_pop", "add__sd",
                                "omega_ka__x", "omega_rx__tcl", "omega_rx__tv",
                                "corr_rx__tcl_ka__x"),
                    value=c(1.5, 1, 3.5, 0.25, 0.7, 0.8, 0.5, 0.3, 0.1))
  ind <- data.frame(id=1:2,
                    eta_ka__x_SAEM=c(0.1, -0.1),
                    eta_rx__tcl_SAEM=c(0.2, -0.2),
                    eta_rx__tv_SAEM=c(0.3, -0.3))
  local_mocked_bindings(
    rxUiGet.monolixPopulationParameters=function(x, ...) pop,
    rxUiGet.monolixIndividualParameters=function(x, ...) ind,
    rxUiGet.monolixIndividualLL=function(x, ...) NULL,
    .package="babelmixr2")

  theta <- u$monolixFullTheta
  expect_equal(theta[["tka"]], log(1.5))
  expect_equal(theta[["cl.wt"]], 0.25)
  expect_false(anyNA(theta))

  omega <- u$monolixOmega
  expect_equal(omega["eta.ka", "eta.ka"], 0.8^2)
  expect_equal(omega["eta.ka", "eta.cl"], 0.1 * 0.8 * 0.5)

  etaObf <- u$monolixEtaObf
  expect_equal(etaObf$eta.ka, c(0.1, -0.1))
  expect_equal(etaObf$eta.v, c(0.3, -0.3))

  # a mu-referenced covariate coefficient on a dotted parameter
  expect_equal(.b$.monolixGetPopParValue(
    "ka.age", u$muRefCurEval, c(tka="ka__x"),
    data.frame(theta="tka", covariate="AGE", covariateParameter="ka.age"),
    data.frame(parameter="beta_ka__x_AGE", value=0.02)), 0.02)
})
