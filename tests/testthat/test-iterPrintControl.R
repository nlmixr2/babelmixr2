test_that("nlmerControl absorbs print/printNcol/useColor into iterPrintControl", {
  .ctl <- nlmerControl()
  expect_s3_class(.ctl$iterPrintControl, "iterPrintControl")
  expect_equal(.ctl$iterPrintControl$every, 0L)
  expect_null(.ctl[["print"]])
  expect_null(.ctl[["printNcol"]])
  expect_null(.ctl[["useColor"]])

  .ctl <- nlmerControl(print = 2L, printNcol = 3L, useColor = FALSE)
  expect_equal(.ctl$iterPrintControl$every, 2L)
  expect_equal(.ctl$iterPrintControl$ncol, 3L)
  expect_false(.ctl$iterPrintControl$useColor)

  .ctl <- nlmerControl(print = nlmixr2est::iterPrintControl(every = 7L))
  expect_equal(.ctl$iterPrintControl$every, 7L)

  # round-trip through do.call() as getValidNlmixrCtl.nlmer does
  .ctl2 <- do.call(nlmerControl, .ctl)
  expect_equal(.ctl2$iterPrintControl, .ctl$iterPrintControl)
})

test_that("saemixControl absorbs print/printNcol/useColor into iterPrintControl", {
  .ctl <- saemixControl()
  expect_s3_class(.ctl$iterPrintControl, "iterPrintControl")
  expect_equal(.ctl$iterPrintControl$every, 0L)
  expect_null(.ctl[["print"]])
  expect_null(.ctl[["printNcol"]])
  expect_null(.ctl[["useColor"]])

  # legacy logical print still works
  .ctl <- saemixControl(print = TRUE)
  expect_equal(.ctl$iterPrintControl$every, 1L)

  .ctl <- saemixControl(print = nlmixr2est::iterPrintControl(every = 5L, ncol = 2L))
  expect_equal(.ctl$iterPrintControl$every, 5L)
  expect_equal(.ctl$iterPrintControl$ncol, 2L)

  # round-trip through do.call() as getValidNlmixrCtl.saemix does
  .ctl2 <- do.call(saemixControl, .ctl)
  expect_equal(.ctl2$iterPrintControl, .ctl$iterPrintControl)
})
