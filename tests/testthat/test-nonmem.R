test_that("NONMEM dsl", {

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
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl/v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }

  ui <- rxode2::rxode2(one.cmt)

  .rxToN <- function(x) {
    rxToNonmem(x, ui)
  }

  expect_equal(.rxToN("sqrt(a)"), "DSQRT(RX001)")
  expect_equal(.rxToN("max(a,b)"), "MAX(RX001,RX002)")
  expect_equal(.rxToN("max(c,b,a)"), "MAX(RX003,RX002,RX001)")
  expect_equal(.rxToN("sum(a,b,c,d)"), "((RX001)+(RX002)+(RX003)+(RX004))")
  expect_equal(.rxToN("prod(a,b,c,d)"), "((RX001)*(RX002)*(RX003)*(RX004))")
  expect_equal(.rxToN("a<-1+b"), "RX001 = 1+RX002")
  expect_equal(.rxToN("a~1+b"), "RX001 = 1+RX002")
  expect_equal(.rxToN("a=1+b"), "RX001 = 1+RX002")
  expect_equal(.rxToN("expit(a)"), "1/(1+DEXP(-(RX001)))")
  expect_equal(.rxToN("expit(a,b)"), "(1.0-(RX002))*(1/(1+DEXP(-(RX001))))+(RX002)")
  expect_equal(.rxToN("expit(a,b,c)"), "((RX003)-(RX002))*(1/(1+DEXP(-(RX001))))+(RX002)")
  expect_equal(.rxToN("logit(a)"), "-DLOG(1/(RX001)-1)")
  expect_equal(.rxToN("logit(a,b)"), "-DLOG(1/(((RX001)-(RX002))/(1.0-(RX002)))-1)")
  expect_equal(.rxToN("logit(a,b,c)"), "-DLOG(1/(((RX001)-(RX002))/((RX003)-(RX002)))-1)")
  expect_equal(.rxToN("probitInv(a)"), "PHI(RX001)")
  expect_equal(.rxToN("probitInv(a,b)"), "(1.0-(RX002))*(PHI(RX001))+(RX002)")
  expect_equal(.rxToN("probitInv(a,b,c)"), "((RX003)-(RX002))*(PHI(RX001))+(RX002)")
  expect_error(.rxToN("probit(a)"))
  expect_error(.rxToN("probit(a,b)"))
  expect_error(.rxToN("probit(a,b,c)"))
  expect_equal(.rxToN("d/dt(depot)=-depot*kel"), "  DADT(1)=- A(1)*KEL")
  expect_equal(.rxToN("f(depot)=3"), ";f defined in PK section")
  expect_equal(.rxToN("a**b"), "RX001**RX002")
  expect_equal(.rxToN("a^b"), "RX001**RX002")
  expect_equal(.rxToN("if (a<=b){c=1} else if (a==4) {c=2} else {c=4}"),
               "  IF (RX001.LE.RX002) THEN\n    RX003=1\n  ELSE IF (RX001.EQ.4) THEN\n    RX003=2\n  ELSE\n    RX003=4\n  END IF\n")
  expect_equal(.rxToN("if (a<=b){c=1} else if (a==4) {c=2} else if (a==30) {c=4} else {c=100}"),
               "  IF (RX001.LE.RX002) THEN\n    RX003=1\n  ELSE IF (RX001.EQ.4) THEN\n    RX003=2\n  ELSE IF (RX001.EQ.30) THEN\n    RX003=4\n  ELSE\n    RX003=100\n  END IF\n")
  expect_equal(.rxToN("if (a<=b){c=1} else if (a==4) {c=2}"),
               "  IF (RX001.LE.RX002) THEN\n    RX003=1\n  ELSE IF (RX001.EQ.4) THEN\n    RX003=2\n  END IF\n")
  expect_equal(.rxToN("if (a<=b){c=1}"), "  IF (RX001.LE.RX002) THEN\n    RX003=1\n  END IF\n")
  expect_equal(.rxToN("time"), "TIME")
  expect_error(.rxToN("NA"))
  expect_error(.rxToN("newind"))
  expect_equal(.rxToN("log1pmx(a)"), "(DLOG(1+RX001)-(RX001))")

  expect_equal(.rxToN("4.3"), "4.3")
  expect_equal(.rxToN("add.sd2"), "ADD_SD2")
  expect_equal(.rxToN("add.sd"), "THETA(6)")

  expect_equal(.rxToN("v.wt"), "THETA(5)")
  expect_equal(.rxToN("eta.cl"), "ETA(2)")

})


pk.turnover.emax3 <- function() {
  ini({
    tktr <- log(1)
    tka <- log(1)
    tcl <- log(0.1)
    tv <- fix(log(10))
    ##
    eta.ktr ~ 1
    eta.ka ~ 1
    eta.cl + eta.v ~ c(2,
                       0.01, 1)
    prop.err <- 0.1
    pkadd.err <- 0.1
    ##
    temax <- logit(0.8)
    tec50 <- log(0.5)
    tkout <- log(0.05)
    te0 <- log(100)
    cl.wt <- c(-10, 0.1, 10)
    cl.sex <- c(-Inf, 0.1, 10)
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
    cl <- exp(tcl + eta.cl + WT * cl.wt + SEXF * cl.sex)
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

w <- pk.turnover.emax3()


