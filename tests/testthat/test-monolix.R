test_that("model to input information", {

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

  uif <- nlmixr(pk.turnover.emax3)

  expect_equal(babelmixr:::monolixMapData(theo_sd, uif),
               list(headerType = c(ID = "id", TIME = "ignore", DV = "observation",
                                   AMT = "amount", EVID = "evid", CMT = "obsid", WT = "ignore"),
                    regressors = "input={v, emax, ec50, e0, kout, ktr, ka, cl}"))


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
      cl0 <- exp(tcl + eta.cl)
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
      cl0 <- cl * (WT / 70) ^ 0.75
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

  uif <- nlmixr(pk.turnover.emax3)


  expect_equal(babelmixr:::monolixMapData(theo_sd, uif),
               list(headerType = c(ID = "id", TIME = "ignore", DV = "observation",
                                   AMT = "amount", EVID = "evid", CMT = "obsid",
                                   WT = "regressor"),
                    regressors = "input={v, emax, ec50, e0, kout, ktr, ka, cl}\nWT = {use=regressor}"))

})



test_that("monolix dsl", {
  expect_equal(rxToMonolix("sqrt(a)"), "sqrt(a)")
  expect_equal(rxToMonolix("max(a,b)"), "max(a,b)")
  expect_error(rxToMonolix("max(a,b,c)"))
  expect_error(rxToMonolix("max(a)"))
  expect_equal(rxToMonolix("sum(a,b,c,d)"), "((a)+(b)+(c)+(d))")
  expect_equal(rxToMonolix("prod(a,b,c,d)"), "((a)*(b)*(c)*(d))")
  expect_equal(rxToMonolix("a<-1+b"), "a = 1+b")
  expect_equal(rxToMonolix("a~1+b"), "a = 1+b")
  expect_equal(rxToMonolix("a=1+b"), "a = 1+b")
  expect_equal(rxToMonolix("expit(a)"), "1/(1+exp(-(a)))")
  expect_equal(rxToMonolix("expit(a,b)"), "(1.0-(b))*(1/(1+exp(-(a))))+(b)")
  expect_equal(rxToMonolix("expit(a,b,c)"), "((c)-(b))*(1/(1+exp(-(a))))+(b)")
  expect_equal(rxToMonolix("logit(a)"), "-log(1/(a)-1)")
  expect_equal(rxToMonolix("logit(a,b)"), "-log(1/(((a)-(b))/(1.0-(b)))-1)")
  expect_equal(rxToMonolix("logit(a,b,c)"), "-log(1/(((a)-(b))/((c)-(b)))-1)")
  expect_equal(rxToMonolix("probitInv(a)"), "normcdf(a)")
  expect_equal(rxToMonolix("probitInv(a,b)"), "(1.0-(b))*(normcdf(a))+(b)")
  expect_equal(rxToMonolix("probitInv(a,b,c)"), "((c)-(b))*(normcdf(a))+(b)")
  expect_equal(rxToMonolix("probit(a)"), "probit(a)")
  expect_equal(rxToMonolix("probit(a,b)"), "probit(((a)-(b))/(1.0-(b)))")
  expect_equal(rxToMonolix("probit(a,b,c)"), "probit(((a)-(b))/((c)-(b)))")
  expect_equal(rxToMonolix("d/dt(depot)=-depot*kel"), "ddt_depot = - depot*kel")
  expect_equal(rxToMonolix("depot(0)=50"), "depot_0 = 50")
  expect_equal(rxToMonolix("f(depot)=3"), ";f defined in PK section")
  expect_equal(rxToMonolix("a**b"), "a^b")
  expect_equal(rxToMonolix("if (a<=b){c=1} else if (a==4) {c=2} else {c=4}"), "if a<=b\n  c = 1\nelseif a==4\n  c = 2\nelse \n  c = 4\nend\n")
  expect_equal(rxToMonolix("if (a<=b){c=1} else if (a==4) {c=2} else if (a==30) {c=4} else {c=100}"), "if a<=b\n  c = 1\nelseif a==4\n  c = 2\nelseif a==30\n  c = 4\nelse \n  c = 100\nend\n")
  expect_equal(rxToMonolix("if (a<=b){c=1} else if (a==4) {c=2}"), "if a<=b\n  c = 1\nelseif a==4\n  c = 2\nend\n")
  expect_equal(rxToMonolix("if (a<=b){c=1}"), "if a<=b\n  c = 1\nend\n")
  expect_equal(rxToMonolix("time"), "t")
  expect_error(rxToMonolix("NA"))
  expect_error(rxToMonolix("newind"))
})
