# Translate models from nlmixr2lib to NONMEM and Monolix; each must
# translate or be refused with a documented error.  By default every
# linCmt() model and a sample of the others are used; set
# BABELMIXR2_STRESS_ALL=true for every model.
test_that("nlmixr2lib models translate to NONMEM/Monolix or are refused", {
  skip_on_cran()
  skip_if_not_installed("nlmixr2lib")
  skip_if_not("linCmtMicro" %in% getNamespaceExports("rxode2"),
              "rxode2 does not have linCmtMicro() for closed-form linCmt()")
  source(system.file("stress", "stress.R", package="babelmixr2"), local=TRUE)
  .all <- identical(tolower(Sys.getenv("BABELMIXR2_STRESS_ALL")), "true")
  withr::with_tempdir({
    res <- stressRun(stressLibCases(all=.all), mode="translate", dir=getwd())
  })
  failed <- res[stressFailed(res), ]
  expect_equal(nrow(failed), 0L,
               info=paste(sprintf("[%s] %s: %s %s", failed$engine, failed$case,
                                  failed$message, failed$problems),
                          collapse="\n"))
})
