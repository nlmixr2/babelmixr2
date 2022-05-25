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
  expect_equal(.rxToN("a<-1+b"), "  RX001=1+RX002")
  expect_equal(.rxToN("a~1+b"), "  RX001=1+RX002")
  expect_equal(.rxToN("a=1+b"), "  RX001=1+RX002")
  expect_equal(.rxToN("expit(a)"), "1/(1+DEXP(-(RX001)))")
  expect_equal(.rxToN("expit(a,b)"), "(1.0-(RX002))*(1/(1+DEXP(-(RX001))))+(RX002)")
  expect_equal(.rxToN("expit(a,b,c)"), "((RX003)-(RX002))*(1/(1+DEXP(-(RX001))))+(RX002)")
  expect_equal(.rxToN("expit(a,1,2)"), "((2.0)-(1.0))*(1/(1+DEXP(-(RX001))))+(1.0)")
  expect_equal(.rxToN("expit(a,0)"), "(1.0-(0.0))*(1/(1+DEXP(-(RX001))))+(0.0)")
  expect_equal(.rxToN("logit(a)"), "-DLOG(1/(RX001)-1)")
  expect_equal(.rxToN("logit(a,b)"), "-DLOG(1/(((RX001)-(RX002))/(1.0-(RX002)))-1)")
  expect_equal(.rxToN("logit(a,0)"), "-DLOG(1/(((RX001)-(0.0))/(1.0-(0.0)))-1)")
  expect_equal(.rxToN("logit(a,b,c)"), "-DLOG(1/(((RX001)-(RX002))/((RX003)-(RX002)))-1)")
  expect_equal(.rxToN("logit(a,0,1)"), "-DLOG(1/(((RX001)-(0.0))/((1.0)-(0.0)))-1)")
  expect_equal(.rxToN("probitInv(a)"), "PHI(RX001)")
  expect_equal(.rxToN("probitInv(a,b)"), "(1.0-(RX002))*(PHI(RX001))+(RX002)")
  expect_equal(.rxToN("probitInv(a,b,c)"), "((RX003)-(RX002))*(PHI(RX001))+(RX002)")
  expect_error(.rxToN("probit(a)"))
  expect_error(.rxToN("probit(a,b)"))
  expect_error(.rxToN("probit(a,b,c)"))
  expect_equal(.rxToN("d/dt(depot)=-depot*kel"), "  DADT(1)=- A(1)*KEL")
  expect_equal(.rxToN("f(depot)=3"), "  ;f defined in $PK block")
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

test_that("tbs tests", {

  pheno <- function() {
    ini({
      tcl <- log(0.008) # typical value of clearance
      tv <-  log(0.6)   # typical value of volume
      ## var(eta.cl)
      eta.cl + eta.v ~ c(1,
                         0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
      # interindividual variability on clearance and volume
      add.err <- 0.1    # residual variability
      lambda <- c(-2, 0, 2)
    })
    model({
      cl <- exp(tcl + eta.cl) # individual value of clearance
      v <- exp(tv + eta.v)    # individual value of volume
      ke <- cl / v            # elimination rate constant
      d/dt(A1) = - ke * A1    # model differential equation
      cp = A1 / v             # concentration in plasma
      cp ~ add(add.err) + boxCox(lambda)# define error model
    })
  }

  p <- pheno()

  expect_equal(p$nonmemCcontra,
               "      subroutine ccontr (icall,c1,c2,c3,ier1,ier2)\n      USE ROCM_REAL,   ONLY: theta=>THETAC,y=>DV_ITM2\n      USE NM_INTERFACE,ONLY: CELS\n!      parameter (lth=40,lvr=30,no=50)\n!      common /rocm0/ theta (lth)\n!      common /rocm4/ y\n!      double precision c1,c2,c3,theta,y,w,one,two\n      double precision c1,c2,c3,w,one,two\n      dimension c2(:),c3(:,:)\n      data one,two/1.,2./\n      if (icall.le.1) return\n      w=y(1)\n      if(theta(4).eq.0) y(1)=log(y(1))\n      if(theta(4).ne.0) y(1)=(y(1)**theta(4)-one)/theta(4)\n      call cels (c1,c2,c3,ier1,ier2)\n      y(1)=w\n      c1=c1-two*(theta(4)-one)*log(y(1))\n      return\n      end")

  expect_equal(p$nonmemErr,
               "$ERROR\n  IPRED = RX_PRED\n  W     = DSQRT(W)\n  IF (W .EQ. 0.0) W = 1.0\n  IF (THETA(4) .EQ. 0.0 .AND. IPRED .NE. 0.0) THEN\n     IPRED = DLOG(IPRED)\n  ELSE IF (THETA(4) .EQ. 0.0 .AND. IPRED .EQ. 0.0) THEN\n     IPRED = -1/THETA(4)\n  ELSE IF (THETA(4) .NE. 0.0 .AND. IPRED .NE. 0.0) THEN\n     IPRED = (IPRED**THETA(4) - 1.0)/THETA(4)\n  ELSE IF (THETA(4) .NE. 0.0 .AND. IPRED .EQ. 0.0) THEN\n     IPRED = -1000000000\n  END IF\n  Y     = IPRED + EPS(1)*W")


   pheno <- function() {
    ini({
      tcl <- log(0.008) # typical value of clearance
      tv <-  log(0.6)   # typical value of volume
      ## var(eta.cl)
      eta.cl + eta.v ~ c(1,
                         0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
      # interindividual variability on clearance and volume
      add.err <- 0.1    # residual variability
      lambda <- c(-2, 0, 2)
    })
    model({
      cl <- exp(tcl + eta.cl) # individual value of clearance
      v <- exp(tv + eta.v)    # individual value of volume
      ke <- cl / v            # elimination rate constant
      d/dt(A1) = - ke * A1    # model differential equation
      cp = A1 / v             # concentration in plasma
      cp ~ add(add.err) + yeoJohnson(lambda)# define error model
    })
  }

  p <- pheno()

  expect_equal(p$nonmemCcontra, "      subroutine ccontr (icall,c1,c2,c3,ier1,ier2)\n      USE ROCM_REAL,   ONLY: theta=>THETAC,y=>DV_ITM2\n      USE NM_INTERFACE,ONLY: CELS\n!      parameter (lth=40,lvr=30,no=50)\n!      common /rocm0/ theta (lth)\n!      common /rocm4/ y\n!      double precision c1,c2,c3,theta,y,w,one,two\n      double precision c1,c2,c3,w,one,two\n      dimension c2(:),c3(:,:)\n      data one,two/1.,2./\n      if (icall.le.1) return\n      w=y(1)\n      if (theta(4) .eq. 1.0) then\n         y(1) = y(1)\n      else if (y(1) .gt. 0.0) then\n         if (theta(4) .eq. 0.0) then\n            y(1) = log(y(1) + one)\n         else\n            y(1) = ((y(1)+one)**theta(4)-one)/theta(4)\n         end if\n      else\n         if (theta(4) .eq. 2.0) then\n            y(1) = -log(one - y(1))\n         else\n            y(1) = (1.0 - (1.0- y(1))**(2.0-theta(4)))/(2.0 - theta(4))\n         end if\n      end if\n      call cels (c1,c2,c3,ier1,ier2)\n      y(1)=w\n      if (y(1) .ge. 0) then\n         c1=c1-two*(theta(4)-one)*log(one+y(1))\n      else\n         c1=c1-two*(one-theta(4))*log(one-y(1))\n      end if\n      return\n      end")

  expect_equal(p$nonmemErr,
               "$ERROR\n  IPRED = RX_PRED\n  W     = DSQRT(W)\n  IF (W .EQ. 0.0) W = 1.0\n  IF (IPRED .GE. 0.0) THEN\n     IF (THETA(4) .EQ. 0.0) THEN\n        IPRED = DLOG(IPRED + 1.0)\n     ELSE IF (THETA(4) .EQ. 1) THEN\n        IPRED = IPRED\n     ELSE\n        IPRED = ((IPRED+1.0)**THETA(4) - 1.0)/THETA(4)\n     END IF \n  ELSE\n     IF (THETA(4) .EQ. 2.0) THEN\n        IPRED = -DLOG(1.0 - IPRED)\n     ELSE IF  (THETA(4) .EQ. 1.0) THEN\n        IPRED = IPRED\n     ELSE\n        IPRED = (1.0 - (1.0 - IPRED)**(2.0 - THETA(4)))/(2.0 - THETA(4))\n     END IF\n  END IF\n  Y = IPRED + EPS(1)*W")

  pheno <- function() {
    ini({
      tcl <- log(0.008) # typical value of clearance
      tv <-  log(0.6)   # typical value of volume
      ## var(eta.cl)
      eta.cl + eta.v ~ c(1,
                         0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
      # interindividual variability on clearance and volume
      lnorm.err <- 0.1    # residual variability
    })
    model({
      cl <- exp(tcl + eta.cl) # individual value of clearance
      v <- exp(tv + eta.v)    # individual value of volume
      ke <- cl / v            # elimination rate constant
      d/dt(A1) = - ke * A1    # model differential equation
      cp = A1 / v             # concentration in plasma
      cp ~ lnorm(lnorm.err)# define error model
    })
  }

  p <- pheno()

  expect_equal(p$nonmemCcontra,
               "      subroutine ccontr (icall,c1,c2,c3,ier1,ier2)\n      USE ROCM_REAL,   ONLY: theta=>THETAC,y=>DV_ITM2\n      USE NM_INTERFACE,ONLY: CELS\n!      parameter (lth=40,lvr=30,no=50)\n!      common /rocm0/ theta (lth)\n!      common /rocm4/ y\n!      double precision c1,c2,c3,theta,y,w,one,two\n      double precision c1,c2,c3,w,one,two\n      dimension c2(:),c3(:,:)\n      data one,two/1.,2./\n      if (icall.le.1) return\n      w=y(1)\n      y(1)=log(y(1))\n      call cels (c1,c2,c3,ier1,ier2)\n      y(1)=w\n      c1=c1+two*log(y(1))\n      return\n      end")

  expect_equal(p$nonmemErr,
               "$ERROR\n  IPRED = RX_PRED\n  W     = DSQRT(W)\n  IF (W .EQ. 0.0) W = 1.0\n  IF (IPRED .GE. 0.0) THEN\n     IF (THETA(4) .EQ. 0.0) THEN\n        IPRED = DLOG(IPRED + 1.0)\n     ELSE IF (THETA(4) .EQ. 1) THEN\n        IPRED = IPRED\n     ELSE\n        IPRED = ((IPRED+1.0)**THETA(4) - 1.0)/THETA(4)\n     END IF \n  ELSE\n     IF (THETA(4) .EQ. 2.0) THEN\n        IPRED = -DLOG(1.0 - IPRED)\n     ELSE IF  (THETA(4) .EQ. 1.0) THEN\n        IPRED = IPRED\n     ELSE\n        IPRED = (1.0 - (1.0 - IPRED)**(2.0 - THETA(4)))/(2.0 - THETA(4))\n     END IF\n  END IF\n  Y = IPRED + EPS(1)*W")


  pheno <- function() {
    ini({
      tcl <- log(0.008) # typical value of clearance
      tv <-  log(0.6)   # typical value of volume
      ## var(eta.cl)
      eta.cl + eta.v ~ c(1,
                         0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
      # interindividual variability on clearance and volume
      lnorm.err <- 0.1    # residual variability
    })
    model({
      cl <- exp(tcl + eta.cl) # individual value of clearance
      v <- exp(tv + eta.v)    # individual value of volume
      ke <- cl / v            # elimination rate constant
      d/dt(A1) = - ke * A1    # model differential equation
      cp = A1 / v             # concentration in plasma
      cp ~ logitNorm(lnorm.err, -0.1, 70)# define error model
    })
  }

  p <- pheno()

  expect_equal(p$nonmemCcontra,
               "      subroutine ccontr (icall,c1,c2,c3,ier1,ier2)\n      USE ROCM_REAL,   ONLY: theta=>THETAC,y=>DV_ITM2\n      USE NM_INTERFACE,ONLY: CELS\n!      parameter (lth=40,lvr=30,no=50)\n!      common /rocm0/ theta (lth)\n!      common /rocm4/ y\n!      double precision c1,c2,c3,theta,y,w,one,two\n      double precision c1,c2,c3,w,one,two,xl,hl,hl2\n      dimension c2(:),c3(:,:)\n      data one,two/1.,2./\n      if (icall.le.1) return\n      w=y(1)\n      xl = (y(1)-(-0.1))/((70.0)-(-0.1))\n      y(1) = -log(one/xl-one)\n      call cels (c1,c2,c3,ier1,ier2)\n      y(1)=w\n      xl = (y(1)-(-0.1))\n      hl = ((70.0) - (-0.1))\n      hl2 = hl-xl\n      c1=c1-two*(log(hl)-log(xl)-log(hl2))\n      return\n      end")

  expect_error(p$nonmemError,
  "$ERROR\n  IPRED = RX_PRED\n  W     = DSQRT(W)\n  XL  = (IPRED - (-0.1))/((70.0) - (-0.1))\n  IPRED = -DLOG(1.0/XL - 1.0)\n  Y     = IPRED + EPS(1)*W")

  pheno <- function() {
    ini({
      tcl <- log(0.008) # typical value of clearance
      tv <-  log(0.6)   # typical value of volume
      ## var(eta.cl)
      eta.cl + eta.v ~ c(1,
                         0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
      # interindividual variability on clearance and volume
      lnorm.err <- 0.1    # residual variability
      lambda <- c(-2, 0, 2)
    })
    model({
      cl <- exp(tcl + eta.cl) # individual value of clearance
      v <- exp(tv + eta.v)    # individual value of volume
      ke <- cl / v            # elimination rate constant
      d/dt(A1) = - ke * A1    # model differential equation
      cp = A1 / v             # concentration in plasma
      cp ~ logitNorm(lnorm.err, -0.1, 70) + yeoJohnson(lambda)# define error model
    })
  }

  p <- pheno()

  expect_equal(p$nonmemErr,
               "$ERROR\n  IPRED = RX_PRED\n  W     = DSQRT(W)\n  XL  = (IPRED - (-0.1))/((70.0) - (-0.1))\n  XL  = -DLOG(1.0/XL - 1.0)\n  IF (THETA(4) .EQ. 1.0) THEN\n     IPRED = XL\n  ELSE IF (XL .GE. 0.0) THEN\n     IF (THETA(4) .EQ. 0.0) THEN\n        IPRED = DLOG(1.0 + XL)\n     ELSE\n        IPRED = ((XL + 1.0)**THETA(4) - 1.0)/THETA(4)\n     END IF\n  ELSE\n     IF (THETA(4) .EQ. 2.0) THEN\n        IPRED = -LOG(1.0 - XL)\n     ELSE\n        HL = 2.0 - THETA(4)\n        IPRED = (1.0 (1.0 - XL)**HL)/HL\n     END IF\n  END IF\n  Y     = IPRED + EPS(1)*W")

  expect_equal(p$nonmemCcontra,
               "      subroutine ccontr (icall,c1,c2,c3,ier1,ier2)\n      USE ROCM_REAL,   ONLY: theta=>THETAC,y=>DV_ITM2\n      USE NM_INTERFACE,ONLY: CELS\n!      parameter (lth=40,lvr=30,no=50)\n!      common /rocm0/ theta (lth)\n!      common /rocm4/ y\n!      double precision c1,c2,c3,theta,y,w,one,two\n      double precision c1,c2,c3,w,one,two,xl,hl,hl2,pd,pdd1,pdd2,p\n      dimension c2(:),c3(:,:)\n      data one,two/1.,2./\n      if (icall.le.1) return\n      w=y(1)\n      xl = (y(1)-(-0.1))/((70.0)-(-0.1))\n      xl = -log(one/xl-one)\n      if (theta(4) .eq. 1.0) then\n         y(1) = xl\n      else if (xl .ge. 0) then\n         if (theta(4) .eq. 0) then\n            y(1) = log(one+xl)\n         else\n            y(1) = ((xl + one)**theta(4) - one)/theta(4)\n         end if\n      else\n         if (theta(4) .eq. 2.0) then\n            y(1) = -log(one-xl)\n         else\n            hl = two - theta(4)\n            y(1) = (one - (one - xl)**hl)/hl\n         end if\n      end if\n      call cels (c1,c2,c3,ier1,ier2)\n      y(1)=w\n      p = (y(1)-(-0.1))/((70.0)-(-0.1))\n      pd = -log(one/p-one)\n      if (theta(4) .eq. one) then\n         pdd1 = 1.0\n      else if (pd .ge. 0) then\n         if (theta(4) .eq. 0.0) then\n            pdd1 = one/(pd + one)\n         else\n            pdd1 = (pd + one)**(theta(4)-1.0)\n         end if\n      else \n         if (theta(4) .eq. 2.0) then\n            pdd1 =  -one/(one - pd);\n         else\n            pdd1 = (one - pd)**(one-theta(4))\n         end if\n      end if\n      xl = (y(1)-(-0.1))\n      hl = ((70.0) - (-0.1));\n      pdd2 = hl/(xl*(hl-xl));\n      c1=c1-two*(log(pdd1)+log(pdd2))\n      return\n      end")

})

