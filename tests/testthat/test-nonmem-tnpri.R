test_that("nonmemTnpri() checks and describes the prior (#206)", {
  expect_error(nonmemTnpri(1), "has to be a dataset")
  expect_error(nonmemTnpri("a.msf", plev=1), "fraction below 1")
  expect_error(nonmemTnpri("a.msf", plev=-0.1))
  expect_error(nonmemTnpri("a.msf", mode=3))
  expect_error(nonmemTnpri("a.msf", display=NA))
  .t <- nonmemTnpri("run/a.msf", mode=2, display=TRUE)
  expect_s3_class(.t, "nonmemTnpri")
  expect_equal(.t$type, "msf")
  # PLEV=0 is what NONMEM asks for when the prior is used for estimation
  expect_equal(.nonmemTnpriOptions(.t), "(PROBLEM 2) PLEV=0 MODE=2 DISPLAY")
  expect_equal(.nonmemTnpriOptions(nonmemTnpri("a.msf")), "(PROBLEM 2) PLEV=0")
  expect_equal(.nonmemTnpriOptions(nonmemTnpri("a.msf", plev=NULL)), "(PROBLEM 2)")
  expect_equal(.nonmemTnpriOptions(nonmemTnpri("a.msf", plev=0.999)),
               "(PROBLEM 2) PLEV=0.999")
  expect_output(print(.t), "run/a.msf")
  # already a nonmemTnpri
  expect_identical(nonmemTnpri(.t), .t)
  expect_equal(nonmemTnpri(data.frame(ID=1))$type, "data")
  # nonmemControl() makes one from what it is given
  expect_equal(nonmemControl(tnpri="a.msf")$tnpri, nonmemTnpri("a.msf"))
  expect_null(nonmemControl()$tnpri)
  expect_false(nonmemControl()$msfo)
  # and keeps it through the nlmixr2 control validation
  .ctl <- getValidNlmixrCtl.nonmem(list(nonmemControl(tnpri="a.msf", msfo=TRUE)))
  expect_equal(.ctl$tnpri, nonmemTnpri("a.msf"))
  expect_true(.ctl$msfo)
})

withr::with_tempdir({

  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1.0
      tv <- 3.45
      add.sd <- 0.7
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
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

  .ui <- function(control) {
    ui <- rxode2::rxUiDecompress(suppressMessages(rxode2::rxode2(one.cmt)))
    assign("control", control, envir=ui)
    ui
  }

  test_that("msfo=TRUE writes the model specification file (#206)", {
    ui <- .ui(nonmemControl(msfo=TRUE))
    expect_true(any(grepl("^\\$ESTIMATION METHOD=1 .* NOABORT MSFO=one.cmt.msf$",
                          strsplit(ui$nonmemModel, "\n")[[1]])))
    ui <- .ui(nonmemControl(msfo=TRUE, est="imp"))
    expect_true(any(grepl("MAPITER=1 NOABORT MSFO=one.cmt.msf$",
                          strsplit(ui$nonmemModel, "\n")[[1]])))
    ui <- .ui(nonmemControl())
    expect_false(grepl("MSFO", ui$nonmemModel))
  })

  test_that("a TNPRI prior gives a two problem control stream (#206)", {
    ui <- .ui(nonmemControl(tnpri="prior/run1.msf"))
    expect_equal(
      ui$nonmemModel,
      paste(
        c(
          "$PROBLEM one.cmt TNPRI prior from run1.msf",
          "; the prior is the estimates and covariance of an earlier fit of this model;",
          "; the model code in this problem is used by every problem",
          "",
          "$DATA one.cmt.csv IGNORE=@",
          "",
          "$INPUT ID TIME EVID AMT DV CMT RXROW",
          "",
          "$SUBROUTINES ADVAN13 TOL=6 ATOL=12 SSTOL=6 SSATOL=12",
          "",
          "$PRIOR TNPRI (PROBLEM 2) PLEV=0",
          "",
          "$MSFI run1.msf ONLYREAD",
          "",
          "$MODEL NCOMPARTMENTS=2",
          "     COMP(DEPOT, DEFDOSE) ; depot",
          "     COMP(CENTRAL) ; central",
          "",
          "$PK",
          "  MU_1=THETA(1)",
          "  MU_2=THETA(2)",
          "  MU_3=THETA(3)",
          "  KA=DEXP(MU_1+ETA(1)) ; ka <- exp(tka)",
          "  CL=DEXP(MU_2+ETA(2)) ; cl <- exp(tcl)",
          "  V=DEXP(MU_3+ETA(3)) ; v <- exp(tv)",
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
          "$PROBLEM one.cmt translated from babelmixr2",
          "; comments show mu referenced model in ui$getSplitMuModel",
          "",
          "$DATA one.cmt.csv IGNORE=@ REWIND",
          "",
          "$INPUT ID TIME EVID AMT DV CMT RXROW",
          "",
          "$THETA (0.45   ) ; 1 - tka   ",
          "       (1      ) ; 2 - tcl   ",
          "       (3.45   ) ; 3 - tv    ",
          "       (0,  0.7) ; 4 - add.sd",
          "",
          "$OMEGA 0.6 ; eta.ka",
          "$OMEGA 0.3 ; eta.cl",
          "$OMEGA 0.1 ; eta.v",
          "",
          "$SIGMA 1 FIX",
          "",
          "$ESTIMATION METHOD=1 INTER MAXEVALS=100000 SIGDIG=3 SIGL=12 PRINT=1 NOABORT",
          "",
          "$COVARIANCE",
          "",
          "$TABLE ID ETAS(1:LAST) OBJI FIRSTONLY ONEHEADER NOPRINT",
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

  test_that("the objective function type says it includes the TNPRI prior (#206)", {
    expect_equal(.ui(nonmemControl(tnpri="a.msf"))$nonmemObjfType, "nonmem focei tnpri")
    expect_equal(.ui(nonmemControl())$nonmemObjfType, "nonmem focei")
  })

  test_that("TNPRI is refused with the NONMEM 7 estimation methods (#206)", {
    for (.est in c("imp", "its")) {
      expect_error(.ui(nonmemControl(est=.est, tnpri="a.msf"))$nonmemModel,
                   paste0("method '", .est, "'"))
    }
    withr::with_tempdir({
      expect_error(
        suppressMessages(
          nlmixr2est::nlmixr(one.cmt, nlmixr2data::Oral_1CPT, est="nonmem",
                             control=nonmemControl(runCommand=NA, est="imp",
                                                   tnpri=nlmixr2data::Oral_1CPT))),
        "method 'imp'")
      # refused before the prior run was written
      expect_false(dir.exists("one.cmt_prior-nonmem"))
    })
    expect_true(grepl("POSTHOC",
                      .ui(nonmemControl(est="posthoc", tnpri="a.msf"))$nonmemModel))
  })

  test_that("ini() priors cannot be combined with a TNPRI prior (#206)", {
    one.cmt.prior <- function() {
      ini({
        tka <- 0.45
        tcl <- 1.0
        tv <- 3.45
        add.sd <- 0.7
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        prior(tka) ~ dnorm(0.45, 1)
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
    ui <- rxode2::rxUiDecompress(suppressMessages(rxode2::rxode2(one.cmt.prior)))
    assign("control", nonmemControl(tnpri="a.msf"), envir=ui)
    expect_error(ui$nonmemModel, "only one \\$PRIOR")
    expect_error(
      suppressMessages(
        nlmixr2est::nlmixr(one.cmt.prior, nlmixr2data::Oral_1CPT, est="nonmem",
                           control=nonmemControl(runCommand=NA, tnpri="a.msf"))),
      "only one \\$PRIOR")
  })

  .d <- nlmixr2data::Oral_1CPT
  .dA <- .d[.d$ID <= 60, ]
  .dB <- .d[.d$ID > 60, ]

  test_that("a dataset prior is exported as its own NONMEM run (#206)", {
    withr::with_tempdir({
      expect_message(
        nlmixr2est::nlmixr(one.cmt, .dB, est="nonmem",
                           control=nonmemControl(runCommand=NA, cov="r",
                                                 tnpri=nonmemTnpri(.dA, mode=1))),
        regexp="run the TNPRI prior 'one.cmt_prior-nonmem/one.cmt_prior.nmctl'")
      .prior <- readLines(file.path("one.cmt_prior-nonmem", "one.cmt_prior.nmctl"))
      # the prior run writes the MSF and has the covariance step TNPRI needs
      expect_true(any(grepl("NOABORT MSFO=one.cmt_prior.msf$", .prior)))
      expect_true("$COVARIANCE MATRIX=R" %in% .prior)
      expect_false(any(grepl("\\$PRIOR|\\$MSFI", .prior)))
      expect_equal(length(unique(read.csv(file.path("one.cmt_prior-nonmem", "one.cmt_prior.csv"))$ID)),
                   length(unique(.dA$ID)))
      .ctl <- readLines(file.path("one.cmt-nonmem", "one.cmt.nmctl"))
      expect_equal(sum(grepl("^\\$PROBLEM", .ctl)), 2L)
      expect_true("$PRIOR TNPRI (PROBLEM 2) PLEV=0 MODE=1" %in% .ctl)
      expect_true("$MSFI one.cmt_prior.msf ONLYREAD" %in% .ctl)
      expect_true("$COVARIANCE MATRIX=R" %in% .ctl)
      expect_false(any(grepl("MSFO", .ctl)))
      expect_equal(length(unique(read.csv(file.path("one.cmt-nonmem", "one.cmt.csv"))$ID)),
                   length(unique(.dB$ID)))
    })
  })

  test_that("a fit exported before its prior was run is found again afterwards (#206)", {
    withr::with_tempdir({
      .fit <- function() {
        suppressMessages(
          nlmixr2est::nlmixr(one.cmt, .dB, est="nonmem",
                             control=nonmemControl(runCommand=NA, tnpri=.dA)))
      }
      .fit()
      # NONMEM run by hand writes the prior's MSF...
      writeLines("msf", file.path("one.cmt_prior-nonmem", "one.cmt_prior.msf"))
      .fit()
      # ...and the fit is still the one that was exported, not a new one
      expect_false(dir.exists("one.cmt-001-nonmem"))
      expect_false(dir.exists("one.cmt_prior-001-nonmem"))
      expect_equal(readLines(file.path("one.cmt-nonmem", "one.cmt_prior.msf")), "msf")
      # while a different prior is a different fit
      suppressMessages(
        nlmixr2est::nlmixr(one.cmt, .dB, est="nonmem",
                           control=nonmemControl(runCommand=NA, tnpri=.dA[.dA$ID <= 30, ])))
      expect_true(dir.exists("one.cmt-001-nonmem"))
    })
  })

  test_that("an MSF prior is copied next to the control stream, and changes the fit (#206)", {
    withr::with_tempdir({
      dir.create("prior")
      writeLines("first", file.path("prior", "run1.msf"))
      .fit <- function() {
        suppressMessages(
          nlmixr2est::nlmixr(one.cmt, .dB, est="nonmem",
                             control=nonmemControl(runCommand=NA,
                                                   tnpri="prior/run1.msf")))
      }
      .fit()
      expect_equal(readLines(file.path("one.cmt-nonmem", "run1.msf")), "first")
      expect_true("$MSFI run1.msf ONLYREAD" %in%
                    readLines(file.path("one.cmt-nonmem", "one.cmt.nmctl")))
      # a different prior is a different fit
      writeLines("second", file.path("prior", "run1.msf"))
      .fit()
      expect_equal(readLines(file.path("one.cmt-001-nonmem", "run1.msf")), "second")
      # a missing MSF is only a warning when NONMEM is not run...
      expect_warning(
        suppressMessages(
          nlmixr2est::nlmixr(one.cmt, .dB, est="nonmem",
                             control=nonmemControl(runCommand=NA,
                                                   modelName="missing",
                                                   tnpri="prior/none.msf"))),
        "does not exist")
      # ...and an error when it is
      expect_error(
        suppressWarnings(suppressMessages(
          nlmixr2est::nlmixr(one.cmt, .dB, est="nonmem",
                             control=nonmemControl(runCommand=function(ctl, directory, ui) NULL,
                                                   modelName="missing2",
                                                   tnpri="prior/none.msf")))),
        "does not exist")
    })
  })

  test_that("a failed prior run stops the fit (#206)", {
    withr::with_tempdir({
      expect_error(
        suppressMessages(
          nlmixr2est::nlmixr(one.cmt, .dB, est="nonmem",
                             control=nonmemControl(runCommand=function(ctl, directory, ui) NULL,
                                                   tnpri=.dA))),
        "creates the TNPRI prior was not successful")
      # the fit itself was never written
      expect_false(dir.exists("one.cmt-nonmem"))
    })
  })

  test_that("a linCmt() fit prior is the same model when NONMEM uses ODEs (#206)", {
    lin.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- 1.0
        tv <- 3.45
        add.sd <- 0.7
        eta.cl ~ 0.3
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        cp <- linCmt()
        cp ~ add(add.sd)
      })
    }
    withr::with_tempdir({
      .f <- suppressMessages(suppressWarnings(
        nlmixr2est::nlmixr(lin.cmt, .dA, est="posthoc")))
      suppressMessages(
        nlmixr2est::nlmixr(lin.cmt, .dB, est="nonmem",
                           control=nonmemControl(runCommand=NA, linCmt="ode", tnpri=.f)))
      expect_true(any(grepl("^\\$DES",
                            readLines(file.path("lin.cmt_prior-nonmem", "lin.cmt_prior.nmctl")))))
      expect_true("$MSFI lin.cmt_prior.msf ONLYREAD" %in%
                    readLines(file.path("lin.cmt-nonmem", "lin.cmt.nmctl")))
    })
  })

  test_that("an nlmixr2 fit prior has to be a fit of the same model (#206)", {
    withr::with_tempdir({
      .f <- suppressMessages(suppressWarnings(
        nlmixr2est::nlmixr(one.cmt, .dA, est="posthoc")))
      suppressMessages(
        nlmixr2est::nlmixr(one.cmt, .dB, est="nonmem",
                           control=nonmemControl(runCommand=NA, tnpri=.f)))
      expect_equal(length(unique(read.csv(file.path("one.cmt_prior-nonmem", "one.cmt_prior.csv"))$ID)),
                   length(unique(.dA$ID)))
      .other <- suppressMessages(rxode2::model(rxode2::rxode2(one.cmt), v <- exp(tv + eta.v + 0.1)))
      expect_error(
        suppressMessages(
          nlmixr2est::nlmixr(.other, .dB, est="nonmem",
                             control=nonmemControl(runCommand=NA, tnpri=.f))),
        "has to be a fit of the same model")
    })
  })
})
