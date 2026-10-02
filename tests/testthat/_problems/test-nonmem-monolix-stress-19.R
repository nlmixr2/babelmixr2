# Extracted from test-nonmem-monolix-stress.R:19

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "babelmixr2", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
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
expect_equal(sort(unique(res$engine)), c("monolix", "nonmem"))
expect_true(all(res$status %in% c("ok", "refused")))
