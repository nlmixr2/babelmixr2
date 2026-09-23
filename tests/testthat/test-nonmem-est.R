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
  .estRecord <- function(est) {
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(f))
    assign("control", nonmemControl(est=est), envir=.ui)
    rxUiGet.nonmemEst(list(.ui))
  }
  expect_match(.estRecord("its"), "^\\$ESTIMATION METHOD=ITS ")
  expect_match(.estRecord("imp"), "^\\$ESTIMATION METHOD=IMP ")
  expect_match(.estRecord("focei"), "^\\$ESTIMATION METHOD=1 ")
  expect_match(.estRecord("posthoc"), "^\\$ESTIMATION METHOD=0 ")
})
