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

  expect_equal(.rxToN("sqrt(a)"), "DSQRT(RXDZ001)")
  expect_equal(.rxToN("max(a,b)"), "MAX(RXR1,RXR2)")
  expect_equal(.rxToN("max(c,b,a)"), "MAX(RXR3,RXR2,RXR1)")
  expect_equal(.rxToN("sum(a,b,c,d)"), "((RXR1)+(RXR2)+(RXR3)+(RXR4))")
  expect_equal(.rxToN("prod(a,b,c,d)"), "((RXR1)*(RXR2)*(RXR3)*(RXR4))")
  expect_equal(.rxToN("a<-1+b"), "              RXR1=1+RXR2")
  expect_equal(.rxToN("a~1+b"), "              RXR1=1+RXR2")
  expect_equal(.rxToN("a=1+b"), "  RXR1=1+RXR2")
  expect_equal(.rxToN("expit(a)"), "1/(1+DEXP(-(RXR1)))")
  expect_equal(.rxToN("expit(a,b)"), "(1.0-(RXR2))*(1/(1+DEXP(-(RXR1))))+(RXR2)")
  expect_equal(.rxToN("expit(a,b,c)"), "((RXR3)-(RXR2))*(1/(1+DEXP(-(RXR1))))+(RXR2)")
  expect_equal(.rxToN("expit(a,1,2)"), "((2.0)-(1.0))*(1/(1+DEXP(-(RXR1))))+(1.0)")
  expect_equal(.rxToN("expit(a,0)"), "(1.0-(0.0))*(1/(1+DEXP(-(RXR1))))+(0.0)")
  expect_equal(.rxToN("logit(a)"), "-DLOG(1/(RXR1)-1)")
  expect_equal(.rxToN("logit(a,b)"), "-DLOG(1/(((RXR1)-(RXR2))/(1.0-(RXR2)))-1)")
  expect_equal(.rxToN("logit(a,0)"), "-DLOG(1/(((RXR1)-(0.0))/(1.0-(0.0)))-1)")
  expect_equal(.rxToN("logit(a,b,c)"), "-DLOG(1/(((RXR1)-(RXR2))/((RXR3)-(RXR2)))-1)")
  expect_equal(.rxToN("logit(a,0,1)"), "-DLOG(1/(((RXR1)-(0.0))/((1.0)-(0.0)))-1)")
  expect_equal(.rxToN("probitInv(a)"), "PHI(RXR1)")
  expect_equal(.rxToN("probitInv(a,b)"), "(1.0-(RXR2))*(PHI(RXR1))+(RXR2)")
  expect_equal(.rxToN("probitInv(a,b,c)"), "((RXR3)-(RXR2))*(PHI(RXR1))+(RXR2)")
  expect_error(.rxToN("probit(a)"))
  expect_error(.rxToN("probit(a,b)"))
  expect_error(.rxToN("probit(a,b,c)"))
  expect_equal(.rxToN("d/dt(depot)=-depot*kel"), "  DADT(1)=- A(1)*KEL")
  expect_equal(.rxToN("f(depot)=3"), "  ;f defined in $PK block")
  expect_equal(.rxToN("a**b"), "RXDZ001**RXR2")
  expect_equal(.rxToN("a^b"), "RXDZ001**RXR2")
  expect_error(
    .rxToN("if (a<=b){c=1} else if (a==4) {c=2} else {c=4}"),
    # Prior result:
    # "  IF (RXR1.LE.RXR2) THEN\n    RXR3=1\n  ELSE IF (RXR1.EQ.4) THEN\n    RXR3=2\n  ELSE\n    RXR3=4\n  END IF\n"
    regexp="babelmixr2 will not allow `else if` or `else` statements in NONMEM models"
  )
  expect_error(
    .rxToN("if (a<=b){c=1} else if (a==4) {c=2} else if (a==30) {c=4} else {c=100}"),
    # Prior result:
    #"  IF (RXR1.LE.RXR2) THEN\n    RXR3=1\n  ELSE IF (RXR1.EQ.4) THEN\n    RXR3=2\n  ELSE IF (RXR1.EQ.30) THEN\n    RXR3=4\n  ELSE\n    RXR3=1R\n  END IF\n"
    regexp="babelmixr2 will not allow `else if` or `else` statements in NONMEM models"
  )
  expect_error(
    .rxToN("if (a<=b){c=1} else if (a==4) {c=2}"),
    # Prior result:
    # "  IF (RXR1.LE.RXR2) THEN\n    RXR3=1\n  ELSE IF (RXR1.EQ.4) THEN\n    RXR3=2\n  END IF\n"
    regexp="babelmixr2 will not allow `else if` or `else` statements in NONMEM models"
  )
  expect_equal(
    .rxToN("if (a<=b){c=1}"),
    paste(c(
    "        IF (RXR1.LE.RXR2) THEN",
    "          RXR3=1",
    "        END IF",
    ""
    ), collapse="\n")
  )
  expect_equal(.rxToN("time"), "TIME")
  expect_error(
    .rxToN("NA"),
    regexp="'NA' cannot be translated to NONMEM"
  )
  expect_error(
    .rxToN("newind"),
    regexp="'newind' cannot be translated to NONMEM"
  )
  expect_equal(.rxToN("log1pmx(a)"), "(DLOG(1+RXR1)-(RXR1))")

  expect_equal(.rxToN("4.3"), "4.3")
  expect_equal(.rxToN("add.sd2"), "ADD_SD2")
  expect_equal(.rxToN("add.sd"), "THETA(6)")

  expect_equal(.rxToN("v.wt"), "THETA(5)")
  expect_equal(.rxToN("eta.cl"), "ETA(2)")
})

# pk.turnover.emax3 <- function() {
#   ini({
#     tktr <- log(1)
#     tka <- log(1)
#     tcl <- log(0.1)
#     tv <- fix(log(10))
#     ##
#     eta.ktr ~ 1
#     eta.ka ~ 1
#     eta.cl + eta.v ~ c(2,
#                        0.01, 1)
#     prop.err <- 0.1
#     pkadd.err <- 0.1
#     ##
#     temax <- logit(0.8)
#     tec50 <- log(0.5)
#     tkout <- log(0.05)
#     te0 <- log(100)
#     cl.wt <- c(-10, 0.1, 10)
#     cl.sex <- c(-Inf, 0.1, 10)
#     ##
#     eta.emax ~ .5
#     eta.ec50  ~ .5
#     eta.kout ~ .5
#     eta.e0 ~ .5
#     ##
#     pdadd.err <- 10
#   })
#   model({
#     ktr <- exp(tktr + eta.ktr)
#     ka <- exp(tka + eta.ka)
#     cl <- exp(tcl + eta.cl + WT * cl.wt + SEXF * cl.sex)
#     v <- exp(tv + eta.v)
#     emax = expit(temax+eta.emax)
#     ec50 =  exp(tec50 + eta.ec50)
#     kout = exp(tkout + eta.kout)
#     e0 = exp(te0 + eta.e0)
#     ##
#     DCP = center/v
#     PD=1-emax*DCP/(ec50+DCP)
#     ##
#     effect(0) = e0
#     kin = e0*kout
#     ##
#     d/dt(depot) = -ktr * depot
#     d/dt(gut) =  ktr * depot -ka * gut
#     d/dt(center) =  ka * gut - cl / v * center
#     d/dt(effect) = kin*PD -kout*effect
#     ##
#     cp = center / v
#     cp ~ prop(prop.err) + add(pkadd.err)
#     effect ~ add(pdadd.err) | pca
#   })
# }
# 
# w <- pk.turnover.emax3()

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

  expect_equal(
    p$nonmemCcontra,
    "      subroutine ccontr (icall,c1,c2,c3,ier1,ier2)\n      USE ROCM_REAL,   ONLY: theta=>THETAC,y=>DV_ITM2\n      USE NM_INTERFACE,ONLY: CELS\n!      parameter (lth=40,lvr=30,no=50)\n!      common /rocm0/ theta (lth)\n!      common /rocm4/ y\n!      double precision c1,c2,c3,theta,y,w,one,two\n      double precision c1,c2,c3,w,one,two\n      dimension c2(:),c3(:,:)\n      data one,two/1.,2./\n      if (icall.le.1) return\n      w=y(1)\n      if(theta(4).eq.0) y(1)=log(y(1))\n      if(theta(4).ne.0) y(1)=(y(1)**theta(4)-one)/theta(4)\n      call cels (c1,c2,c3,ier1,ier2)\n      y(1)=w\n      c1=c1-two*(theta(4)-one)*log(y(1))\n      return\n      end"
  )

  # Confirmed accurate relative to
  # Dosne, AG., Bergstrand, M. & Karlsson, M.O. A strategy for residual error
  # modeling incorporating scedasticity of variance and distribution shape. J
  # Pharmacokinet Pharmacodyn 43, 137–151 (2016).
  # https://doi.org/10.1007/s10928-015-9460-y
  expect_equal(
    p$nonmemErrF,
    paste(c(
      "",
      "  ; Write out expressions for ipred and w",
      "  RX_IP1 = RX_PF1",
      "  IF (THETA(4) .EQ. 0.0 .AND. RX_IP1 .NE. 0.0) THEN",
      "     RX_IP1 = DLOG(RX_IP1)",
      "  ELSE IF (THETA(4) .EQ. 0.0 .AND. RX_IP1 .EQ. 0.0) THEN",
      "     RX_IP1 = -1/THETA(4)",
      "  ELSE IF (THETA(4) .NE. 0.0 .AND. RX_IP1 .NE. 0.0) THEN",
      "     RX_IP1 = (RX_IP1**THETA(4) - 1.0)/THETA(4)",
      "  ELSE IF (THETA(4) .NE. 0.0 .AND. RX_IP1 .EQ. 0.0) THEN",
      "     RX_IP1 = -1000000000",
      "  END IF",
      "  RX_P1 = RX_IP1",
      "  W1=DSQRT((THETA(3))**2)",
      "  IF (W1 .EQ. 0.0) W1 = 1",
      "  IPRED = RX_IP1",
      "  W     = W1",
      "  Y     = IPRED + W*EPS(1)",
      ""
    ), collapse="\n")
  )

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

  expect_equal(
    p$nonmemCcontra,
    "      subroutine ccontr (icall,c1,c2,c3,ier1,ier2)\n      USE ROCM_REAL,   ONLY: theta=>THETAC,y=>DV_ITM2\n      USE NM_INTERFACE,ONLY: CELS\n!      parameter (lth=40,lvr=30,no=50)\n!      common /rocm0/ theta (lth)\n!      common /rocm4/ y\n!      double precision c1,c2,c3,theta,y,w,one,two\n      double precision c1,c2,c3,w,one,two\n      dimension c2(:),c3(:,:)\n      data one,two/1.,2./\n      if (icall.le.1) return\n      w=y(1)\n      if (theta(4) .eq. 1.0) then\n         y(1) = y(1)\n      else if (y(1) .gt. 0.0) then\n         if (theta(4) .eq. 0.0) then\n            y(1) = log(y(1) + one)\n         else\n            y(1) = ((y(1)+one)**theta(4)-one)/theta(4)\n         end if\n      else\n         if (theta(4) .eq. 2.0) then\n            y(1) = -log(one - y(1))\n         else\n            y(1) = (1.0 - (1.0- y(1))**(2.0-theta(4)))/(2.0 - theta(4))\n         end if\n      end if\n      call cels (c1,c2,c3,ier1,ier2)\n      y(1)=w\n      if (y(1) .ge. 0) then\n         c1=c1-two*(theta(4)-one)*log(one+y(1))\n      else\n         c1=c1-two*(one-theta(4))*log(one-y(1))\n      end if\n      return\n      end"
  )

  # Confirmed as an accurate Yeo-Johnson transformation by cross-checking with
  # https://en.wikipedia.org/wiki/Power_transform#Yeo%E2%80%93Johnson_transformation
  expect_equal(
    p$nonmemErrF,
    paste(c(
      "",
      "  ; Write out expressions for ipred and w",
      "  RX_IP1 = RX_PF1", 
      "  IF (RX_IP1 .GE. 0.0) THEN",
      "     IF (THETA(4) .EQ. 0.0) THEN",
      "        RX_IP1 = DLOG(RX_IP1 + 1.0)",
      "     ELSE IF (THETA(4) .EQ. 1.0) THEN",
      "        RX_IP1 = RX_IP1",
      "     ELSE",
      "        RX_IP1 = ((RX_IP1+1.0)**THETA(4) - 1.0)/THETA(4)",
      "     END IF ",
      "  ELSE",
      "     IF (THETA(4) .EQ. 2.0) THEN",
      "        RX_IP1 = -DLOG(1.0 - RX_IP1)",
      "     ELSE IF  (THETA(4) .EQ. 1.0) THEN",
      "        RX_IP1 = RX_IP1",
      "     ELSE",
      "        RX_IP1 = (1.0 - (1.0 - RX_IP1)**(2.0 - THETA(4)))/(2.0 - THETA(4))",
      "     END IF",
      "  END IF",
      "  RX_P1 = RX_IP1",
      "  W1=DSQRT((THETA(3))**2)",
      "  IF (W1 .EQ. 0.0) W1 = 1",
      "  IPRED = RX_IP1",
      "  W     = W1",
      "  Y     = IPRED + W*EPS(1)",
      ""
    ), collapse="\n")
  )

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

  expect_equal(
    p$nonmemErrF,
    paste(c(
      "",
      "  ; Write out expressions for ipred and w",
      "  RX_IP1 = RX_PF1",
      "  IF (RX_IP1 .EQ. 0.0) THEN",
      "     RX_IP1 = -1000000000",
      "  ELSE",
      "     RX_IP1 = DLOG(RX_IP1)",
      "  END IF",
      "  RX_P1 = RX_IP1",
      "  W1=DSQRT((THETA(3))**2)",
      "  IF (W1 .EQ. 0.0) W1 = 1",
      "  IPRED = RX_IP1",
      "  W     = W1",
      "  Y     = IPRED + W*EPS(1)",
      ""
    ), collapse="\n")
  )

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

  expect_equal(
    p$nonmemErrF,
    paste(c(
      "",
      "  ; Write out expressions for ipred and w",
      "  RX_IP1 = RX_PF1",
      "  XL  = (RX_IP1 - (-0.1))/((70.0) - (-0.1))",
      "  RX_IP1 = -DLOG(1.0/XL - 1.0)",
      "  RX_P1 = RX_IP1",
      "  W1=DSQRT((THETA(3))**2)",
      "  IF (W1 .EQ. 0.0) W1 = 1",
      "  IPRED = RX_IP1",
      "  W     = W1",
      "  Y     = IPRED + W*EPS(1)",
      ""
    ), collapse="\n")
  )

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

  expect_equal(p$nonmemCcontra,
               "      subroutine ccontr (icall,c1,c2,c3,ier1,ier2)\n      USE ROCM_REAL,   ONLY: theta=>THETAC,y=>DV_ITM2\n      USE NM_INTERFACE,ONLY: CELS\n!      parameter (lth=40,lvr=30,no=50)\n!      common /rocm0/ theta (lth)\n!      common /rocm4/ y\n!      double precision c1,c2,c3,theta,y,w,one,two\n      double precision c1,c2,c3,w,one,two,xl,hl,hl2,pd,pdd1,pdd2,p\n      dimension c2(:),c3(:,:)\n      data one,two/1.,2./\n      if (icall.le.1) return\n      w=y(1)\n      xl = (y(1)-(-0.1))/((70.0)-(-0.1))\n      xl = -log(one/xl-one)\n      if (theta(4) .eq. 1.0) then\n         y(1) = xl\n      else if (xl .ge. 0) then\n         if (theta(4) .eq. 0) then\n            y(1) = log(one+xl)\n         else\n            y(1) = ((xl + one)**theta(4) - one)/theta(4)\n         end if\n      else\n         if (theta(4) .eq. 2.0) then\n            y(1) = -log(one-xl)\n         else\n            hl = two - theta(4)\n            y(1) = (one - (one - xl)**hl)/hl\n         end if\n      end if\n      call cels (c1,c2,c3,ier1,ier2)\n      y(1)=w\n      p = (y(1)-(-0.1))/((70.0)-(-0.1))\n      pd = -log(one/p-one)\n      if (theta(4) .eq. one) then\n         pdd1 = 1.0\n      else if (pd .ge. 0) then\n         if (theta(4) .eq. 0.0) then\n            pdd1 = one/(pd + one)\n         else\n            pdd1 = (pd + one)**(theta(4)-1.0)\n         end if\n      else \n         if (theta(4) .eq. 2.0) then\n            pdd1 =  -one/(one - pd);\n         else\n            pdd1 = (one - pd)**(one-theta(4))\n         end if\n      end if\n      xl = (y(1)-(-0.1))\n      hl = ((70.0) - (-0.1));\n      pdd2 = hl/(xl*(hl-xl));\n      c1=c1-two*(log(pdd1)+log(pdd2))\n      return\n      end")

  expect_equal(
    p$nonmemErrF,
    paste(c(
      "",
      "  ; Write out expressions for ipred and w",
      "  RX_IP1 = RX_PF1",
      "  XL  = (RX_IP1 - (-0.1))/((70.0) - (-0.1))",
      "  XL  = -DLOG(1.0/XL - 1.0)",
      "  IF (THETA(4) .EQ. 1.0) THEN",
      "     RX_IP1 = XL",
      "  ELSE IF (XL .GE. 0.0) THEN",
      "     IF (THETA(4) .EQ. 0.0) THEN",
      "        RX_IP1 = DLOG(1.0 + XL)",
      "     ELSE",
      "        RX_IP1 = ((XL + 1.0)**THETA(4) - 1.0)/THETA(4)",
      "     END IF",
      "  ELSE",
      "     IF (THETA(4) .EQ. 2.0) THEN",
      "        RX_IP1 = -DLOG(1.0 - XL)",
      "     ELSE",
      "        HL = 2.0 - THETA(4)",
      "        RX_IP1 = (1.0 - (1.0 - XL)**HL)/HL",
      "     END IF",
      "  END IF",
      "  RX_P1 = RX_IP1",
      "  W1=DSQRT((THETA(3))**2)",
      "  IF (W1 .EQ. 0.0) W1 = 1",
      "  IPRED = RX_IP1",
      "  W     = W1",
      "  Y     = IPRED + W*EPS(1)",
      ""
    ), collapse="\n")
  )
})

test_that("NONMEM WBC model", {
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
  
  w <- wbc()
  
  expect_equal(
    w$nonmemErrF,
    paste(c(
      "",
      "  ; Write out expressions for ipred and w",
      "  RX_IP1 = RX_PF1",
      "  RX_P1 = RX_IP1",
      "  W1=DSQRT((RX_PF1*THETA(5))**2)",
      "  IF (W1 .EQ. 0.0) W1 = 1",
      "  IPRED = RX_IP1",
      "  W     = W1",
      "  Y     = IPRED + W*EPS(1)",
      ""
    ), collapse="\n")
  )
})

test_that("NONMEM multiple endpoint methods", {
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
      cpadd.sd <- 0.1
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
      lambda1 <- fix(0.5)
      lambda2 <- fix(-0.5)
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
      cp ~ add(cpadd.sd) + boxCox(lambda1)
      effect ~ add(pdadd.err) + yeoJohnson(lambda2)| pca
    })
  }

  w <- pk.turnover.emax3()
  
  expect_equal(
    w$nonmemErrF,
    paste(c(
      "",
      "  ; Write out expressions for ipred and w",
      "  RX_IP1 = RX_PF1",
      "  IF (THETA(10) .EQ. 0.0 .AND. RX_IP1 .NE. 0.0) THEN",
      "     RX_IP1 = DLOG(RX_IP1)",
      "  ELSE IF (THETA(10) .EQ. 0.0 .AND. RX_IP1 .EQ. 0.0) THEN",
      "     RX_IP1 = -1/THETA(10)",
      "  ELSE IF (THETA(10) .NE. 0.0 .AND. RX_IP1 .NE. 0.0) THEN",
      "     RX_IP1 = (RX_IP1**THETA(10) - 1.0)/THETA(10)",
      "  ELSE IF (THETA(10) .NE. 0.0 .AND. RX_IP1 .EQ. 0.0) THEN",
      "     RX_IP1 = -1000000000",
      "  END IF",
      "  RX_P1 = RX_IP1",
      "  W1=DSQRT((THETA(5))**2)",
      "  IF (W1 .EQ. 0.0) W1 = 1",
      "  RX_IP2 = RX_PF2",
      "  IF (RX_IP2 .GE. 0.0) THEN",
      "     IF (THETA(11) .EQ. 0.0) THEN",
      "        RX_IP2 = DLOG(RX_IP2 + 1.0)",
      "     ELSE IF (THETA(11) .EQ. 1.0) THEN",
      "        RX_IP2 = RX_IP2",
      "     ELSE",
      "        RX_IP2 = ((RX_IP2+1.0)**THETA(11) - 1.0)/THETA(11)",
      "     END IF ",
      "  ELSE",
      "     IF (THETA(11) .EQ. 2.0) THEN",
      "        RX_IP2 = -DLOG(1.0 - RX_IP2)",
      "     ELSE IF  (THETA(11) .EQ. 1.0) THEN",
      "        RX_IP2 = RX_IP2",
      "     ELSE",
      "        RX_IP2 = (1.0 - (1.0 - RX_IP2)**(2.0 - THETA(11)))/(2.0 - THETA(11))",
      "     END IF",
      "  END IF",
      "  RX_P2 = RX_IP2",
      "  W2=DSQRT((THETA(12))**2)",
      "  IF (W2 .EQ. 0.0) W2 = 1",
      "  IPRED = RX_IP1",
      "  W     = W1",
      "  IF (DVID .EQ. 2) THEN",
      "    IPRED = RX_IP2",
      "    W     = W2",
      "  END IF",
      "  Y     = IPRED + W*EPS(1)",
      ""
    ), collapse="\n")
  )
})
