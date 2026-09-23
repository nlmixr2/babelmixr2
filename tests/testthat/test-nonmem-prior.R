test_that("est='nonmem' declares NWPRI prior support (#205)", {
  expect_equal(attr(nlmixr2Est.nonmem, "nlmixr2Priors"), "nwpri")
})

withr::with_tempdir({

  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1.0
      tv <- 3.45
      add.sd <- 0.7
      eta.cl + eta.v ~ c(0.3,
                         0.01, 0.1)
      eta.ka ~ 0.6
      tka ~ 0.01
      prior(tcl) ~ dnorm(1.1, 0.2)
      prior(eta.cl, eta.v) ~ invWishart(200)
      prior(eta.ka) ~ invWishart(4)
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl/v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }

  test_that("ini() priors become $PRIOR NWPRI records (#205)", {
    ui <- suppressMessages(one.cmt())
    expect_equal(
      ui$nonmemModel,
      paste(
        c(
          "$PROBLEM  translated from babelmixr2",
          "; comments show mu referenced model in ui$getSplitMuModel",
          "",
          "$DATA one.cmt.csv IGNORE=@",
          "",
          "$INPUT ID TIME EVID AMT DV CMT RXROW",
          "",
          "$SUBROUTINES ADVAN13 TOL=6 ATOL=12 SSTOL=6 SSATOL=12",
          "",
          "$PRIOR NWPRI",
          "",
          "$MODEL NCOMPARTMENTS=2",
          "     COMP(DEPOT, DEFDOSE) ; depot",
          "     COMP(CENTRAL) ; central",
          "",
          "$PK",
          "  MU_3=THETA(1)",
          "  MU_1=THETA(2)",
          "  MU_2=THETA(3)",
          "  KA=DEXP(MU_3+ETA(3)) ; ka <- exp(tka)",
          "  CL=DEXP(MU_1+ETA(1)) ; cl <- exp(tcl)",
          "  V=DEXP(MU_2+ETA(2)) ; v <- exp(tv)",
          "",
          "$DES",
          "  DADT(1) = - KA*A(1) ; d/dt(depot) = -ka * depot",
          "  DADT(2) = KA*A(1)-CL/V*A(2) ; d/dt(central) = ka * depot - cl/v * central",
          "  CP=A(2)/V ; cp = central/v",
          "",
          "$ERROR",
          "  ;Redefine LHS in $DES by prefixing with on RXE_ for $ERROR",
          "  RXE_CP=A(2)/V ; cp = central/v",
          "  RX_PF1=RXE_CP ; rx_pf1 ~ cp",
          "  ; Write out expressions for ipred and w",
          "  RX_IP1 = RX_PF1",
          "  RX_P1 = RX_IP1",
          "  W1=DSQRT((THETA(4))**2) ; W1 ~ sqrt((add.sd)^2)",
          "  IF (W1 .EQ. 0.0) W1 = 1",
          "  IPRED = RX_IP1",
          "  W     = W1",
          "  Y     = IPRED + W*EPS(1)",
          "",
          "$THETA (0.45   ) ; 1 - tka   ",
          "       (1      ) ; 2 - tcl   ",
          "       (3.45   ) ; 3 - tv    ",
          "       (0,  0.7) ; 4 - add.sd",
          "",
          "$OMEGA BLOCK(2) ; eta.cl eta.v",
          "   0.3",
          "   0.01 0.1",
          "$OMEGA 0.6 ; eta.ka",
          "",
          "$SIGMA 1 FIX",
          "",
          "$THETAP (0.45 FIX) ; 1 - tka",
          "        (1.1 FIX ) ; 2 - tcl",
          "",
          "$THETAPV BLOCK(2) ; tka tcl",
          "   0.01",
          "   0 0.04  FIX",
          "",
          "$OMEGAP BLOCK(2) ; eta.cl eta.v",
          "   0.3",
          "   0.01 0.1  FIX",
          "$OMEGAPD (200 FIX) ; eta.cl eta.v",
          "$OMEGAP 0.6 FIX ; eta.ka",
          "$OMEGAPD (4 FIX) ; eta.ka",
          "",
          "$ESTIMATION METHOD=1 INTER MAXEVALS=10000 SIGDIG=3 SIGL=12 PRINT=1 NOABORT",
          "",
          "$COVARIANCE",
          "",
          "$TABLE ID ETAS(1:3) OBJI FIRSTONLY ONEHEADER NOPRINT",
          "     FORMAT=s1PE17.9 NOAPPEND FILE=one.cmt.eta",
          "",
          "$TABLE ID TIME IPRED PRED RXROW ONEHEADER NOPRINT",
          "    FORMAT=s1PE17.9 NOAPPEND FILE=one.cmt.pred",
          ""
        ),
        collapse="\n"
      )
    )
  })

  test_that("nlmixr(est='nonmem') exports a model with priors (#205)", {
    expect_message(
      nlmixr2est::nlmixr(one.cmt, data=nlmixr2data::Oral_1CPT, est="nonmem",
                         control=nonmemControl(runCommand=NA)),
      regexp="only exported NONMEM"
    )
    .ctl <- readLines(file.path("one.cmt-nonmem", "one.cmt.nmctl"))
    expect_true("$PRIOR NWPRI" %in% .ctl)
    expect_true(any(grepl("^\\$THETAP ", .ctl)))
    expect_true(any(grepl("^\\$OMEGAPD ", .ctl)))
  })

  thetaOnly <- function() {
    ini({
      tcl <- 1.0
      tv <- 3.45
      tka <- 0.45
      add.sd <- 0.7
      eta.cl ~ 0.3
      tcl + tv ~ c(0.1,
                   0.01, 0.2)
      prior(tka) ~ stdNormal()
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl/v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }

  test_that("joint and standard normal theta priors, no omega prior (#205)", {
    ui <- suppressMessages(thetaOnly())
    expect_equal(
      ui$nonmemPriorRecords,
      paste(
        c(
          "$THETAP (1 FIX   ) ; 1 - tcl",
          "        (3.45 FIX) ; 2 - tv ",
          "        (0 FIX   ) ; 3 - tka",
          "",
          "$THETAPV BLOCK(3) ; tcl tv tka",
          "   0.1",
          "   0.01 0.2",
          "   0 0 1  FIX",
          "",
          ""
        ),
        collapse="\n"
      )
    )
    expect_false(grepl("OMEGAP", ui$nonmemModel))
    # no omega prior: every eta is still in the model
    expect_true(grepl("ETAS(1:1)", ui$nonmemModel, fixed=TRUE))
  })

  test_that("models without priors are unchanged (#205)", {
    noPrior <- function() {
      ini({
        tcl <- 1.0
        eta.cl ~ 0.3
        add.sd <- 0.7
      })
      model({
        cl <- exp(tcl + eta.cl)
        d/dt(central) <- - cl * central
        cp <- central
        cp ~ add(add.sd)
      })
    }
    ui <- suppressMessages(noPrior())
    expect_null(ui$nonmemPriorSpec)
    expect_equal(ui$nonmemPrior, "")
    expect_equal(ui$nonmemPriorRecords, "")
    expect_false(grepl("PRIOR|THETAP|OMEGAP", ui$nonmemModel))
    expect_true(grepl("ETAS(1:LAST)", ui$nonmemModel, fixed=TRUE))
  })

  test_that("$OMEGA and $OMEGAP blocks are written row by row (#205)", {
    block3 <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.7
        eta.ka + eta.cl + eta.v ~ c(1,
                                    0.1, 2,
                                    0.2, 0.3, 3)
        prior(eta.ka, eta.cl, eta.v) ~ invWishart(10)
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl/v * central
        cp <- central / v
        cp ~ add(add.sd)
      })
    }
    ui <- suppressMessages(block3())
    expect_equal(ui$nonmemOmega,
                 paste0("$OMEGA BLOCK(3) ; eta.ka eta.cl eta.v\n",
                        "   1\n",
                        "   0.1 2\n",
                        "   0.2 0.3 3\n"))
    expect_true(grepl(paste0("$OMEGAP BLOCK(3) ; eta.ka eta.cl eta.v\n",
                             "   1\n",
                             "   0.1 2\n",
                             "   0.2 0.3 3  FIX\n"),
                      ui$nonmemPriorRecords, fixed=TRUE))
    # an omega prior alone: no theta prior records
    expect_false(grepl("THETAP", ui$nonmemModel, fixed=TRUE))
    expect_true(grepl("$OMEGAPD (10 FIX) ; eta.ka eta.cl eta.v", ui$nonmemModel, fixed=TRUE))
  })

  base <- function() {
    ini({
      tcl <- 1.0
      tv <- 3.45
      tka <- 0.45
      eta.cl + eta.v ~ c(0.3,
                         0.01, 0.1)
      eta.ka ~ 0.6
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl/v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }

  test_that("priors NWPRI cannot express are refused (#205)", {
    .ui <- function(pri) suppressMessages(eval(bquote(ini(base, .(substitute(pri))))))
    # NWPRI puts its priors on the first THETAs
    expect_error(.ui(prior(tv) ~ dnorm(3, 1))$nonmemModel,
                 "move the population parameter\\(s\\) with a prior \\('tv'\\) before the one\\(s\\) without \\('tcl'\\)")
    expect_error(.ui(prior(add.sd) ~ dnorm(0.7, 0.1))$nonmemModel,
                 "first THETAs")
    # ... and on the first omega blocks
    expect_error(.ui(prior(eta.ka) ~ invWishart(3))$nonmemModel,
                 "move the eta\\(s\\) with a prior \\('eta.ka'\\) before the one\\(s\\) without \\('eta.cl', 'eta.v'\\)")
    # a second prior on a block that already has one is not dropped
    expect_error(.ui({
      prior(eta.cl) ~ invWishart(10)
      prior(eta.v) ~ invWishart(20)
    })$nonmemModel,
    "already covered by the prior 'invWishart\\(10\\)'")
    # a prior whose parameters cannot be read back
    expect_error(.ui(prior(tcl) ~ dnorm(undefinedVariable, 0.2))$nonmemModel,
                 "its parameters could not be read back")
    # no Cauchy analogue
    expect_error(.ui(prior(tcl) ~ dcauchy(0, 1))$nonmemModel,
                 "'dcauchy\\(\\)' is not a normal prior")
    # a normal prior on an omega element is TNPRI, not NWPRI
    expect_error(.nonmemPriorSpec(.ui(prior(eta.cl) ~ dnorm(0.3, 0.1))),
                 "only 'invWishart\\(nu\\)' can be given to an omega block")
  })

  test_that("a joint normal prior reaching an omega element is refused (#205)", {
    # rxode2's omega-normal assertion does not see a joint block, so
    # this is the only thing refusing it
    joint <- function() {
      ini({
        tcl <- 1
        tv <- 3
        add.sd <- 0.7
        eta.cl ~ 0.1
        tcl + om.eta.cl ~ c(0.1,
                            0.01, 0.2)
      })
      model({
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d/dt(central) <- - cl/v * central
        cp <- central/v
        cp ~ add(add.sd)
      })
    }
    expect_error(suppressMessages(joint())$nonmemModel,
                 "a joint normal prior can only be on population parameters")
  })

  test_that("nlmixr(est='nonmem') refuses priors before writing files (#205)", {
    withr::with_tempdir({
      .fit <- function(pri) {
        .m <- suppressMessages(eval(bquote(ini(base, .(substitute(pri))))))
        nlmixr2est::nlmixr(.m, data=nlmixr2data::Oral_1CPT,
                           est="nonmem", control=nonmemControl(runCommand=NA))
      }
      expect_error(.fit(prior(tcl) ~ dcauchy(0, 1)), "dcauchy")
      expect_error(.fit(prior(eta.cl) ~ dnorm(0.3, 0.1)), "normal prior on the omega")
      expect_error(.fit(prior(tv) ~ dnorm(3, 1)), "first THETAs")
      expect_length(list.files(recursive=TRUE), 0L)
    })
  })

  test_that("NWPRI's own THETAs/OMEGAs are dropped when reading NONMEM output (#205)", {
    pheno <- function() {
      ini({
        tcl <- log(0.008)
        tv <-  log(0.6)
        add.err <- 0.1
        eta.cl + eta.v ~ c(1,
                           0.01, 1)
        prior(tcl) ~ dnorm(-4.8, 1)
        prior(eta.cl, eta.v) ~ invWishart(10)
      })
      model({
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        ke <- cl / v
        d/dt(A1) = - ke * A1
        cp = A1 / v
        cp ~ add(add.err)
      })
    }
    ui <- rxode2::rxUiDecompress(suppressMessages(pheno()))
    rxode2::rxAssignControlValue(ui, "modelName", "pheno")
    # NM-TRAN appends the prior values to the model's THETAs and OMEGAs:
    # here THETAP and OMEGAPD as THETA(4)/THETA(5), THETAPV and OMEGAP
    # as a 3x3 OMEGA after the model's 2x2 block
    .theta <- c(-4.9, -0.5, 2.8, -4.8, 10)
    .omega <- lotri::lotri(eta1 + eta2 ~ c(0.4, 0.05, 0.2),
                           eta3 ~ 1,
                           eta4 + eta5 ~ c(1, 0.01, 1))
    rxode2::rxAssignControlValue(ui, ".lstInfo",
                                 list(theta=setNames(.theta, paste0("theta", 1:5)),
                                      omega=.omega, objf=100))
    expect_equal(ui$nonmemFullTheta, c(tcl=-4.9, tv=-0.5, add.err=2.8))
    expect_equal(ui$nonmemOutputOmega,
                 matrix(c(0.4, 0.05, 0.05, 0.2), 2, 2,
                        dimnames=list(c("eta.cl", "eta.v"), c("eta.cl", "eta.v"))))
    expect_equal(ui$nonmemObjfType, "nonmem focei nwpri")

    .path <- ui$nonmemExportPath
    dir.create(.path)
    .om <- paste0("OMEGA(", unlist(lapply(1:5, function(i) paste0(i, ",", seq_len(i)))), ")")
    .cols <- c(paste0("THETA", 1:5), "SIGMA(1,1)", .om)
    .n <- length(.cols)
    .fmt <- function(x) formatC(x, format="E", digits=5, width=13)
    # .cov: an identifiable value in every cell
    .cov <- outer(seq_len(.n), seq_len(.n), function(i, j) 1 / (i + j))
    writeLines(c("TABLE NO.     1: First Order Conditional Estimation with Interaction: Problem=1",
                 paste(c(" NAME        ", sprintf("%-13s", .cols)), collapse=""),
                 vapply(seq_len(.n), function(i) {
                   paste0(" ", sprintf("%-12s", .cols[i]), paste(.fmt(.cov[i, ]), collapse=""))
                 }, character(1))),
               file.path(.path, "pheno.cov"))
    .c <- ui$nonmemCovariance
    expect_equal(dimnames(.c), list(c("tcl", "tv", "add.err"), c("tcl", "tv", "add.err")))
    expect_equal(unname(.c), .cov[1:3, 1:3], tolerance=1e-5)

    # .ext: iteration history keeps only the model parameters
    .val <- c(.theta, 1,
              unlist(lapply(1:5, function(i) .omega[i, seq_len(i)])))
    writeLines(c("TABLE NO.     1: First Order Conditional Estimation with Interaction: Goal Function=MINIMUM VALUE OF OBJECTIVE FUNCTION: Problem=1",
                 paste(c(" ITERATION   ", sprintf("%-13s", c(.cols, "OBJ"))), collapse=""),
                 paste0(sprintf("%13d", 1L), paste(.fmt(.val), collapse=""), "    100.5"),
                 paste0(sprintf("%13d", -1000000000L), paste(.fmt(.val), collapse=""), "    100.5")),
               file.path(.path, "pheno.ext"))
    .h <- ui$nonmemParHistory
    expect_equal(names(.h),
                 c("iter", "tcl", "tv", "add.err", "eta.cl", "(eta.cl,eta.v)", "eta.v",
                   "objf", "type"))
    expect_equal(.h$iter, 1)
    expect_equal(unlist(.h[1, c("tcl", "tv", "add.err", "eta.cl", "(eta.cl,eta.v)", "eta.v")],
                        use.names=FALSE),
                 c(-4.9, -0.5, 2.8, 0.4, 0.05, 0.2))
  })
})
