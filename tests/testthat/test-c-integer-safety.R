# tests/testthat/test-c-integer-safety.R
#
# Tests for integer type safety and bounds-checking fixes in C++ source.
# All tests use only exported C++ functions via RcppExports.R; no PopED or
# other suggested packages are required.

# ---------------------------------------------------------------------------
# popedGetMultipleEndpointModelingTimes / popedMultipleEndpointIndexDataFrame
# ---------------------------------------------------------------------------

test_that("popedGetMultipleEndpointModelingTimes returns correct sorted unique times", {
  times       <- c(1.1, 1.2, 1.3, 2.1, 2.2, 3.1)
  modelSwitch <- c(1L,  1L,  1L,  2L,  2L,  3L)

  result <- popedGetMultipleEndpointModelingTimes(times, modelSwitch, sorted = TRUE)

  expect_type(result, "double")
  expect_equal(result, sort(unique(times)))

  popedMultipleEndpointResetTimeIndex()
})

test_that("popedGetMultipleEndpointModelingTimes handles overlapping times across endpoints", {
  times       <- c(1.1, 1.2, 1.3, 0.5, 2.2, 1.1, 0.75, 0.75)
  modelSwitch <- c(1L,  1L,  1L,  2L,  2L,  2L,  3L,   3L)

  result <- popedGetMultipleEndpointModelingTimes(times, modelSwitch, sorted = TRUE)

  expect_type(result, "double")
  expect_equal(result, sort(unique(times)))

  popedMultipleEndpointResetTimeIndex()
})

test_that("popedGetMultipleEndpointModelingTimes unsorted returns unique times", {
  times       <- c(3.0, 1.0, 2.0)
  modelSwitch <- c(1L, 1L, 1L)

  result_unsorted <- popedGetMultipleEndpointModelingTimes(times, modelSwitch, sorted = FALSE)
  result_sorted   <- popedGetMultipleEndpointModelingTimes(times, modelSwitch, sorted = TRUE)

  # Both should contain the same values
  expect_setequal(result_unsorted, unique(times))
  # Sorted version should be in order
  expect_equal(result_sorted, sort(unique(times)))

  popedMultipleEndpointResetTimeIndex()
})

test_that("popedMultipleEndpointIndexDataFrame returns correct data.frame structure", {
  times       <- c(1.1, 1.2, 1.3, 2.1, 2.2, 3.1)
  modelSwitch <- c(1L,  1L,  1L,  2L,  2L,  3L)

  popedGetMultipleEndpointModelingTimes(times, modelSwitch, sorted = TRUE)
  df <- popedMultipleEndpointIndexDataFrame()

  expect_s3_class(df, "data.frame")
  expect_true("time" %in% names(df))
  # Rows equal number of sorted unique times
  expect_equal(nrow(df), length(sort(unique(times))))
  # 3 model switches → 7 columns: time, MS:1, N:1, MS:2, N:2, MS:3, N:3
  expect_equal(ncol(df), 7L)

  popedMultipleEndpointResetTimeIndex()
})

test_that("popedMultipleEndpointIndexDataFrame stops on modelSwitch value of 0", {
  # info.id = 0 violates the 'info.id <= 0' guard
  times_bad       <- c(1.0, 2.0)
  modelSwitch_bad <- c(0L, 1L)

  popedGetMultipleEndpointModelingTimes(times_bad, modelSwitch_bad, sorted = FALSE)
  expect_error(
    popedMultipleEndpointIndexDataFrame(),
    regexp = "modelSwitch need to be sequential"
  )

  popedMultipleEndpointResetTimeIndex()
})

test_that("popedMultipleEndpointResetTimeIndex clears the global indexer", {
  times       <- c(1.0, 2.0)
  modelSwitch <- c(1L, 2L)

  popedGetMultipleEndpointModelingTimes(times, modelSwitch, sorted = FALSE)
  popedMultipleEndpointResetTimeIndex()

  # After reset, the indexer is not initialized → error
  expect_error(
    popedMultipleEndpointIndexDataFrame(),
    regexp = "has not been initialized"
  )
})

# ---------------------------------------------------------------------------
# popedMultipleEndpointParam
# ---------------------------------------------------------------------------

test_that("popedMultipleEndpointParam returns numeric vector of correct length", {
  p           <- c(1.0, 2.0, 3.0)   # first element is id and is dropped
  times       <- c(0.5, 1.5, 2.5)
  modelSwitch <- c(1L, 2L, 3L)
  maxMT       <- 5L

  result <- popedMultipleEndpointParam(p, times, modelSwitch, maxMT, optTime = FALSE)

  expect_type(result, "double")
  # Output length = (length(p) - 1) + maxMT = 2 + 5 = 7
  expect_length(result, length(p) - 1L + maxMT)

  popedMultipleEndpointResetTimeIndex()
})

test_that("popedMultipleEndpointParam fills excess time slots with max observed time", {
  p           <- c(99.0, 2.0, 3.0)
  times       <- c(0.5, 1.5)
  modelSwitch <- c(1L, 2L)
  maxMT       <- 5L

  result <- popedMultipleEndpointParam(p, times, modelSwitch, maxMT, optTime = FALSE)

  # Slots beyond the unique time count are padded with max(times) = 1.5
  n_params   <- length(p) - 1L  # 2
  n_utimes   <- length(unique(times))  # 2
  fill_start <- n_params + n_utimes + 1L
  if (fill_start <= length(result)) {
    expect_true(all(result[fill_start:length(result)] == max(times)))
  }

  popedMultipleEndpointResetTimeIndex()
})

# ---------------------------------------------------------------------------
# convertDataBack
# ---------------------------------------------------------------------------

test_that("convertDataBack handles evid=0 (observations) and returns correct fields", {
  result <- convertDataBack(
    id   = c(1L, 1L, 1L),
    time = c(0.0, 1.0, 2.0),
    amt  = c(0.0, 0.0, 0.0),
    ii   = c(0.0, 0.0, 0.0),
    evid = c(0L, 0L, 0L),
    cmt  = c(1L, 1L, 1L),
    cmtDvid  = integer(0),
    dvidDvid = integer(0)
  )

  expect_type(result, "list")
  expected_fields <- c("df", "adm", "turnOffCmt", "hasTinf", "hasRate",
                        "hasPhantom", "hasReplace", "hasMult",
                        "hasSs", "hasSs2", "hasSsRate", "nobs")
  expect_true(all(expected_fields %in% names(result)))
  expect_equal(result$nobs, 3L)
  expect_true(all(result$df$.nlmixrKeep))
})

test_that("convertDataBack handles evid=2 (non-observation) rows", {
  result <- convertDataBack(
    id   = c(1L, 1L),
    time = c(0.0, 1.0),
    amt  = c(0.0, 0.0),
    ii   = c(0.0, 0.0),
    evid = c(2L, 0L),
    cmt  = c(1L, 1L),
    cmtDvid  = integer(0),
    dvidDvid = integer(0)
  )

  df <- result$df
  expect_true(all(df$.nlmixrKeep))
  expect_equal(result$nobs, 1L)
})

test_that("convertDataBack handles evid=3 (ODE reset) rows", {
  result <- convertDataBack(
    id   = c(1L, 1L),
    time = c(0.0, 1.0),
    amt  = c(0.0, 0.0),
    ii   = c(0.0, 0.0),
    evid = c(3L, 0L),
    cmt  = c(1L, 1L),
    cmtDvid  = integer(0),
    dvidDvid = integer(0)
  )

  expect_equal(result$df$EVID[1], 3L)
  expect_true(result$df$.nlmixrKeep[1])
})

test_that("convertDataBack bolus dose (cmt=1) creates an adm entry", {
  # EVID encoding for cmt=1 bolus dose (EVIDF_NORMAL, EVID0_REGULAR):
  #   evid = wh100*1e5 + whI*1e4 + wh0*100 + xx
  # For cmt=1: wh100=0, whI=0, wh0=1, xx=1 → evid = 0 + 0 + 100 + 1 = 101
  result <- convertDataBack(
    id   = 1L,
    time = 0.0,
    amt  = 10.0,
    ii   = 0.0,
    evid = 101L,
    cmt  = 1L,
    cmtDvid  = integer(0),
    dvidDvid = integer(0)
  )

  expect_s3_class(result$adm, "data.frame")
  expect_equal(nrow(result$adm), 1L)
  expect_equal(result$df$EVID[1], 1L)
  expect_false(result$turnOffCmt)
})

test_that("convertDataBack reports correct boolean flags", {
  result <- convertDataBack(
    id   = 1L,
    time = 0.0,
    amt  = 0.0,
    ii   = 0.0,
    evid = 0L,
    cmt  = 1L,
    cmtDvid  = integer(0),
    dvidDvid = integer(0)
  )

  expect_false(result$turnOffCmt)
  expect_false(result$hasTinf)
  expect_false(result$hasRate)
  expect_false(result$hasPhantom)
  expect_false(result$hasReplace)
  expect_false(result$hasMult)
  expect_false(result$hasSs)
  expect_false(result$hasSs2)
  expect_false(result$hasSsRate)
})

# ---------------------------------------------------------------------------
# transDv
# ---------------------------------------------------------------------------

test_that("transDv applies Box-Cox transformation (lambda=1 is identity)", {
  inDv    <- c(1.0, 2.0, 3.0)
  inCmt   <- c(1L, 1L, 1L)
  cmtTrans <- 1L
  lambda   <- 1.0
  yj       <- 1L
  low      <- -Inf
  high     <- Inf

  result <- transDv(
    inDv = inDv, inCmt = inCmt,
    cmtTrans = cmtTrans, lambda = lambda,
    yj = yj, low = low, high = high
  )

  expect_type(result, "list")
  expect_true(all(c("dv", "dvid", "cmt", "nCmt", "likAdj") %in% names(result)))
  expect_length(result$dv, length(inDv))
  # All observations match cmt=1, so dvid should all be 1
  expect_true(all(result$dvid == 1L, na.rm = TRUE))
  # nCmt should equal total number of matching observations
  expect_equal(result$nCmt[1], length(inDv))
})

test_that("transDv returns NA for compartments not in cmtTrans", {
  inDv     <- c(1.0, 2.0)
  inCmt    <- c(99L, 99L)
  cmtTrans <- 1L
  lambda   <- 1.0
  yj       <- 1L
  low      <- -Inf
  high     <- Inf

  result <- transDv(
    inDv = inDv, inCmt = inCmt,
    cmtTrans = cmtTrans, lambda = lambda,
    yj = yj, low = low, high = high
  )

  expect_true(all(is.na(result$dv)))
  expect_true(all(is.na(result$dvid)))
  expect_equal(result$cmt, inCmt)
  expect_equal(result$nCmt[1], 0L)
})

test_that("transDv handles empty input gracefully", {
  result <- transDv(
    inDv = numeric(0), inCmt = integer(0),
    cmtTrans = 1L, lambda = 1.0,
    yj = 1L, low = -Inf, high = Inf
  )

  expect_length(result$dv, 0L)
  expect_length(result$dvid, 0L)
  expect_equal(result$nCmt[1], 0L)
})

test_that("transDv handles multiple compartments", {
  inDv     <- c(1.0, 2.0, 3.0, 4.0)
  inCmt    <- c(1L, 2L, 1L, 2L)
  cmtTrans <- c(1L, 2L)
  lambda   <- c(1.0, 1.0)
  yj       <- c(1L, 1L)
  low      <- c(-Inf, -Inf)
  high     <- c(Inf, Inf)

  result <- transDv(
    inDv = inDv, inCmt = inCmt,
    cmtTrans = cmtTrans, lambda = lambda,
    yj = yj, low = low, high = high
  )

  expect_equal(result$nCmt[1], 2L)  # cmt=1 appears twice
  expect_equal(result$nCmt[2], 2L)  # cmt=2 appears twice
  expect_true(all(!is.na(result$dvid)))
})

# ---------------------------------------------------------------------------
# Large-vector test: integer overflow safety (R_xlen_t fix)
# ---------------------------------------------------------------------------
#
# Pre-fix crash mechanism for getDvid() in classicReverse.cpp:
#
#   for (int j = cmtDvid.size(); j--;)
#
#   When cmtDvid.size() = 2^31 (INT_MAX + 1 = 2,147,483,648):
#     int j = 2^31 → INT_MIN = -2,147,483,648  (signed overflow, UB)
#     j-- returns INT_MIN (nonzero → truthy), then sets j = INT_MAX = 2,147,483,647
#     Body: cmtDvid[INT_MAX] → out-of-bounds read → SEGFAULT
#
#   Post-fix: R_xlen_t j = 2^31 holds the value correctly, counts down
#   from 2^31 - 1 = 2,147,483,647 to 0 without overflow.
#
# Memory requirements:
#   Triggering the actual crash requires cmtDvid.size() >= 2^31.
#   One integer vector of 2^31 elements = 4 bytes × 2,147,483,648 = 8 GB.
#   That exceeds typical machine limits; the test below uses n = 6e8
#   (600 million, ~2.4 GB for a single integer vector) which is below INT_MAX
#   and therefore demonstrates correct R_xlen_t handling but does NOT trigger
#   the pre-fix crash.  To reproduce the crash with the unfixed code, set
#   n <- 2^31 on a machine with > 8 GB free RAM, then comment out skip().
#
test_that("R_xlen_t loop variable: convertDataBack handles large cmtDvid correctly", {
  skip("Requires >8 GB RAM: cmtDvid vector of 2^31 integers (8 GB) needed to trigger pre-fix segfault; n=6e8 used here fits in ~2.4 GB and verifies correct R_xlen_t behavior for sub-INT_MAX sizes")

  # Using n = 6e8 (~2.4 GB for one integer vector).
  # Allocate only cmtDvid as a large vector; dvidDvid can stay small because
  # cmt=1L will not match any entry in an all-zero cmtDvid, so dvidDvid is
  # never accessed.  This caps total RAM at ~2.4 GB.
  #
  # To test the actual pre-fix crash (requires >8 GB free RAM):
  #   Replace 6e8 with 2^31 and comment out skip().
  n <- 6e8  # change to 2^31 for the actual crash test (requires >8 GB RAM)
  big_cmtDvid <- integer(n)   # all zeros; ~2.4 GB

  result <- convertDataBack(
    id   = 1L,
    time = 0.0,
    amt  = 0.0,
    ii   = 0.0,
    evid = 0L,    # observation
    cmt  = 1L,    # cmt=1 won't match any 0 in big_cmtDvid → dvidDvid never accessed
    cmtDvid  = big_cmtDvid,
    dvidDvid = integer(0)   # safe: never accessed since cmt=1 ≠ 0
  )

  # Post-fix: processes the 1 observation correctly, dvid = 0 (no match)
  expect_equal(result$nobs, 1L)
  expect_true(result$df$.nlmixrKeep[1])
  expect_equal(result$df$DVID[1], 0L)
  rm(big_cmtDvid)
  gc()
})
