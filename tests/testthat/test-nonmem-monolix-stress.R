# Edge cases of the nlmixr2 -> NONMEM and nlmixr2 -> Monolix translations
# (translation only; see inst/stress/README.md to run NONMEM/Monolix end
# to end with the same cases)
test_that("NONMEM/Monolix stress test cases translate or are refused", {
  skip_on_cran()
  skip_if_not("linCmtMicro" %in% getNamespaceExports("rxode2"),
              "rxode2 does not have linCmtMicro() for closed-form linCmt()")
  source(system.file("stress", "stress.R", package="babelmixr2"), local=TRUE)
  withr::with_tempdir({
    res <- stressRun(stressCases(), mode="translate", dir=getwd())
  })
  failed <- res[stressFailed(res), ]
  expect_equal(nrow(failed), 0L,
               info=paste(sprintf("[%s] %s: %s %s", failed$engine, failed$case,
                                  failed$message, failed$problems),
                          collapse="\n"))
  # every case ran for both engines
  expect_equal(sort(unique(res$engine)), c("monolix", "nonmem"))
  expect_true(all(res$status %in% c("ok", "refused")))
})
