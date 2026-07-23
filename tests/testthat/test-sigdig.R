# Tolerance / standard-option alignment with nlmixr2est.  These assert the
# babelmixr2 in-R estimation controls key their ODE + optimizer tolerances to
# `sigdig` exactly the way the nlmixr2est methods do, and that the standard
# options (indTolRelax, eventSens, sigdigTable) reach the generated
# foceiControl.  Pure control construction -- no fit, so they run on CRAN.

test_that("sigdig scales the ODE tolerances the nlmixr2est way", {
  # rtol = 10^-sigdig, atol = 10^(-sigdig-3); sens matches the main solve,
  # steady-state one order looser
  for (.f in list(nlmerControl, fmeMcmcControl, pseudoOptimControl, saemixControl)) {
    .c4 <- .f(sigdig = 4)$rxControl
    expect_equal(.c4$rtol, 1e-4)
    expect_equal(.c4$atol, 1e-7)
    expect_equal(.c4$rtolSens, 1e-4)
    expect_equal(.c4$atolSens, 1e-7)
    expect_equal(.c4$ssRtol, 1e-3)
    expect_equal(.c4$ssAtol, 1e-6)

    .c6 <- .f(sigdig = 6)$rxControl
    expect_equal(.c6$rtol, 1e-6)
    expect_equal(.c6$atol, 1e-9)
    expect_equal(.c6$rtolSens, 1e-6)
    expect_equal(.c6$ssRtol, 1e-5)
  }
})

test_that("babelmixr2 sigdig tolerances match nlmixr2est's own controls", {
  # the babelmixr2 helpers mirror the (unexported) nlmixr2est ones, so a
  # default-built rxControl must be byte-identical to nlmixr2est's
  .fields <- c("rtol", "atol", "rtolSens", "atolSens", "ssRtol", "ssAtol")
  for (.sig in c(3, 4, 6)) {
    .ref <- nlmixr2est::bobyqaControl(sigdig = .sig)$rxControl[.fields]
    for (.f in list(nlmerControl, fmeMcmcControl, pseudoOptimControl, saemixControl)) {
      expect_equal(.f(sigdig = .sig)$rxControl[.fields], .ref)
    }
  }
})

test_that("an explicit rxControl list overrides sigdig for the fields it names", {
  # skip mechanism: a user atol wins, everything else still follows sigdig
  .c <- fmeMcmcControl(sigdig = 4, rxControl = list(atol = 1e-12))$rxControl
  expect_equal(.c$atol, 1e-12)   # user value preserved
  expect_equal(.c$rtol, 1e-4)    # still sigdig-derived
})

test_that("sigdigTable follows sigdig unless given", {
  expect_equal(nlmerControl(sigdig = 6)$sigdigTable, 6)
  expect_equal(fmeMcmcControl(sigdig = 5)$sigdigTable, 5)
  expect_equal(pseudoOptimControl(sigdig = 4)$sigdigTable, 4)
  expect_equal(saemixControl(sigdig = 5)$sigdigTable, 5)
  # explicit value wins
  expect_equal(saemixControl(sigdig = 5, sigdigTable = 2)$sigdigTable, 2)
})

test_that("tolPwrss and varleft are derived from sigdig", {
  # tuned default preserved at sigdig=4, one order per significant digit
  expect_equal(nlmerControl(sigdig = 4)$tolPwrss, 1e-7)
  expect_equal(nlmerControl(sigdig = 6)$tolPwrss, 1e-9)
  expect_equal(pseudoOptimControl(sigdig = 4)$varleft, 1e-8)
  expect_equal(pseudoOptimControl(sigdig = 6)$varleft, 1e-10)
  # an explicit value still wins
  expect_equal(nlmerControl(tolPwrss = 1e-5)$tolPwrss, 1e-5)
  expect_equal(pseudoOptimControl(varleft = 1e-3)$varleft, 1e-3)
})

test_that("indTolRelax and eventSens are real options with the expected defaults", {
  for (.f in list(nlmerControl, fmeMcmcControl, pseudoOptimControl, saemixControl)) {
    .c <- .f()
    expect_true(.c$indTolRelax)
    expect_identical(.c$eventSens, "jump")
    expect_false(.f(indTolRelax = FALSE)$indTolRelax)
    expect_identical(.f(eventSens = "fd")$eventSens, "fd")
  }
})

test_that("indTolRelax and eventSens propagate into the generated foceiControl", {
  # the final table build runs through a foceiControl; the standard options
  # must reach it rather than being dropped/hard-coded
  .mkEnv <- function(ctlName, control) {
    .e <- new.env(parent = emptyenv())
    .e[[ctlName]] <- control
    .e$ui <- list()
    .e
  }

  .fc <- babelmixr2:::.fmeMcmcControlToFoceiControl(
    .mkEnv("fmeMcmcControl", fmeMcmcControl(indTolRelax = FALSE, eventSens = "fd")),
    assign = FALSE)
  expect_false(.fc$indTolRelax)
  expect_identical(.fc$eventSens, "fd")

  .fc <- babelmixr2:::.pseudoOptimControlToFoceiControl(
    .mkEnv("pseudoOptimControl", pseudoOptimControl(indTolRelax = FALSE, eventSens = "fd")),
    assign = FALSE)
  expect_false(.fc$indTolRelax)
  expect_identical(.fc$eventSens, "fd")
})
