test_that("nonmemControl(est=) emits the matching $ESTIMATION method (#211)", {
  f <- function() {
    ini({
      tcl <- 1
      eta.cl ~ 0.1
      add.sd <- 0.7
    })
    model({
      cl <- exp(tcl + eta.cl)
      cp <- cl * t
      cp ~ add(add.sd)
    })
  }
  .ui <- function(est, ...) {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(f))
    assign("control", nonmemControl(est=est, ...), envir=.ui)
    list(.ui)
  }
  .check <- function(est, method, ofvType) {
    .x <- .ui(est)
    expect_match(rxUiGet.nonmemEst(.x),
                 paste0("^\\$ESTIMATION METHOD=", method, " "))
    expect_equal(rxUiGet.nonmemObjfType(.x), ofvType)
  }
  .check("its", "ITS", "nonmem its")
  .check("imp", "IMP", "nonmem imp")
  .check("focei", "1", "nonmem focei")
  .check("posthoc", "0", "nonmem focei")
  .its <- "$ESTIMATION METHOD=ITS INTERACTION PRINT=1 NITER=100 NOABORT\n"
  expect_equal(rxUiGet.nonmemEst(.ui("its")), .its)
  expect_equal(rxUiGet.nonmemEst(.ui("its", niter=50, print=5, noabort=FALSE)),
               "$ESTIMATION METHOD=ITS INTERACTION PRINT=5 NITER=50\n")
  # the importance sampling options do not apply to ITS
  expect_equal(rxUiGet.nonmemEst(.ui("its", seed=1, isample=5000, iaccept=0.5,
                                     iscaleMin=0.2, iscaleMax=5, df=2,
                                     mapiter=3)),
               .its)
})

test_that("IMP/ITS/posthoc NONMEM runs are read as successful", {
  f <- function() {
    ini({
      tcl <- 1
      eta.cl ~ 0.1
      add.sd <- 0.7
    })
    model({
      cl <- exp(tcl + eta.cl)
      cp <- cl * t
      cp ~ add(add.sd)
    })
  }
  .ok <- function(est, term) {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(f))
    assign("control", nonmemControl(est = est), envir = .ui)
    local_mocked_bindings(rxUiGet.nonmemTermMessage = function(x, ...) term)
    rxUiGet.nonmemSuccessful(list(.ui))
  }
  .notTested <- "\n OPTIMIZATION WAS NOT TESTED FOR CONVERGENCE\n\n"
  .success <- "\n0MINIMIZATION SUCCESSFUL\n NO. OF FUNCTION EVALUATIONS USED: 123\n"
  expect_true(.ok("focei", .success))
  expect_false(.ok(
    "focei",
    "\n0MINIMIZATION TERMINATED\n DUE TO ROUNDING ERRORS (ERROR=134)\n"
  ))
  # IMP/ITS without a convergence test
  expect_true(.ok("imp", .notTested))
  expect_true(.ok("its", .notTested))
  expect_true(.ok("its", "\n OPTIMIZATION WAS COMPLETED\n"))
  expect_false(.ok("imp", "\n OPTIMIZATION WAS NOT COMPLETED\n"))
  # MAXEVALS=0 writes no #TERM: block
  expect_true(.ok("posthoc", NA_character_))
  expect_false(.ok("focei", NA_character_))
  expect_false(.ok("focei", NULL))
})

test_that("a posthoc run has no NONMEM parameter history", {
  f <- function() {
    ini({
      tcl <- 1
      eta.cl ~ 0.1
      add.sd <- 0.7
    })
    model({
      cl <- exp(tcl + eta.cl)
      cp <- cl * t
      cp ~ add(add.sd)
    })
  }
  .ui <- rxode2::rxUiDecompress(rxode2::rxode2(f))
  assign("control", nonmemControl(est = "posthoc"), envir = .ui)
  withr::local_dir(withr::local_tempdir())
  # only the final estimates and their summaries (negative iterations)
  writeLines(
    c(
      "TABLE NO.     1: First Order (Evaluation): Goal Function=MINIMUM VALUE OF OBJECTIVE FUNCTION: Problem=1 Subproblem=0 Superproblem1=0 Iteration1=0 Superproblem2=0 Iteration2=0",
      " ITERATION    THETA1       THETA2       SIGMA(1,1)   OMEGA(1,1)   OBJ",
      "  -1000000000  1.00000E+00  7.00000E-01  1.00000E+00  1.00000E-01    141.25637",
      "  -1000000001  1.00000E-01  1.00000E-01  1.00000E+10  1.00000E-02    0.00000"
    ),
    "posthoc.ext"
  )
  local_mocked_bindings(
    rxUiGet.nonmemExportPath = function(x, ...) ".",
    rxUiGet.nonmemExt = function(x, ...) "posthoc.ext"
  )
  expect_null(rxUiGet.nonmemParHistory(list(.ui)))
})
