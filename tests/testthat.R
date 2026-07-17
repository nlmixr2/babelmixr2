library(testthat)
library(babelmixr2)

# Test thread policy: cap rxode2 threads to 2 on CRAN (two-core policy) and
# force serial testthat workers on CI/CRAN (a single parallel
# worker only adds the callr message pipe, which base64-serializes large
# objects inlined in an erroring test's backtrace and can blow R's 2GB string
# cap).  Locally, everything is managed normally.
.onCran <- !identical(Sys.getenv("NOT_CRAN"), "true")
.onCI   <- isTRUE(as.logical(Sys.getenv("CI", "false")))

if (.onCI || .onCran) {
  options(Ncpus = 1L)
  Sys.setenv(TESTTHAT_CPUS = "1")
  Sys.setenv(TESTTHAT_PARALLEL = "FALSE")
}
if (.onCran && requireNamespace("rxode2", quietly = TRUE)) {
  rxode2::setRxThreads(2L)
}

# -------------------------------------------------------------------------
# Fast / slow test partitioning (mirrors nlmixr2est's tests/testthat.R)
#
# The estimation tests (full SAEM/nlmer/pseudoOptim/PopED fits and the
# nonmem2rx/monolix2rx file conversions) dominate the suite wall time.  On
# push/PR CI, R-CMD-check runs only the "essential" (fast) subset -- every
# test file EXCEPT the slow ones in .slowBatches below.  The slow files run
# separately, one batch at a time, in a slow-tests workflow that sets
# BABELMIXR2_TEST_BATCH=<n>.
#
# Names are the test file basename with the leading "test-" and trailing
# ".R" removed (what testthat's `filter` matches).  Locally and on CRAN the
# filter is NULL so everything runs (on CRAN the estimation files skip
# themselves via skip_on_cran()).
# -------------------------------------------------------------------------
.slowBatches <- list(
  # batch 1 -- the SAEM comparison fits (direct saemix + nlmixr2 saemix),
  # by far the heaviest file.
  c("saemix"),
  # batch 2 -- nonmem2rx/monolix2rx conversions + PopED design fits.
  c("nonmem-read", "monolix-read", "poped", "as-nlmixr2")
)
.slowAll <- unlist(.slowBatches)

.batch <- Sys.getenv("BABELMIXR2_TEST_BATCH")

.filter <- NULL
if (nzchar(.batch)) {
  # Slow-batch mode: run ONLY this batch's slow files.
  .b <- suppressWarnings(as.integer(.batch))
  if (is.na(.b) || .b < 1L || .b > length(.slowBatches)) {
    stop(sprintf("BABELMIXR2_TEST_BATCH=%s out of range (1..%d)",
                 .batch, length(.slowBatches)))
  }
  .files <- .slowBatches[[.b]]
  .filter <- if (length(.files) == 0L) {
    "^$"  # empty batch: match nothing
  } else {
    paste0("^(", paste(.files, collapse = "|"), ")$")
  }
} else if (.onCI && !.onCran && length(.slowAll) > 0L) {
  # Essential subset on push/PR CI: everything EXCEPT the slow files.
  .filter <- paste0("^(?!(", paste(.slowAll, collapse = "|"), ")$)")
}
# Locally (and on CRAN) .filter stays NULL -> run everything.

test_check("babelmixr2", stop_on_failure = FALSE, filter = .filter, perl = TRUE)
