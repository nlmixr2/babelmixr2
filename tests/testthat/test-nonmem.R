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
  expect_equal(.rxToN("depot(0)=50"), "A_0(1) = 50")
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
  expect_equal(.rxToN("add.sd"), "add_sd")

})
