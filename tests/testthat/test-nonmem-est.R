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
  # the importance sampling options do not apply to ITS
  expect_equal(rxUiGet.nonmemEst(.ui("its", seed=1, isample=5000, iaccept=0.5,
                                     iscaleMin=0.2, iscaleMax=5, df=2,
                                     mapiter=3)),
               .its)
})
