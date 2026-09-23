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
  .ui <- function(est) {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(f))
    assign("control", nonmemControl(est=est), envir=.ui)
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
  # ITS must not carry the importance sampling options
  expect_false(grepl("SEED|ISAMPLE|IACCEPT|ISCALE_MIN|ISCALE_MAX|DF=|MAPITER", rxUiGet.nonmemEst(.ui("its"))))
})
