.nonmem2rx <- function(...) suppressWarnings(suppressMessages(nonmem2rx::nonmem2rx(...)))

.as.nonmem2rx <- function(...) suppressWarnings(suppressMessages(nonmem2rx::as.nonmem2rx(...)))

test_that("nonmem2rx imports solve with rxControl(nonmem = TRUE) (#252)", {
  .ode <- rxode2::rxode2({
    cl <- 3 * (1 + 0.05 * time)
    d / dt(central) <- -cl / 30 * central
  })
  expect_equal(.nonmem2rxUseNonmemSolve(.ode), .nonmem2rxHasNonmemSolve())
  # rxode2 treats a delay() outside of d/dt() as a $PK-type statement, so
  # delay models keep the continuous time
  .dde <- rxode2::rxode2({
    cdel <- delay(central, 4)
    d / dt(central) <- -0.1 * central
    d / dt(resp) <- 10 * (1 - cdel / (1 + cdel)) - 0.5 * resp
    resp(0) <- 20
    past(central, 4) <- 0
  })
  expect_false(.nonmem2rxUseNonmemSolve(.dde))

  skip_on_cran()
  mod <- .nonmem2rx(system.file("mods/cpt/runODE032.ctl", package = "nonmem2rx"),
                    determineError = FALSE, lst = ".res", save = FALSE)
  .ctl <- .nonmem2rxToFoceiControl(new.env(), rxode2::rxUiDecompress(mod))
  expect_equal(isTRUE(.ctl$rxControl$nonmem), .nonmem2rxHasNonmemSolve())
  expect_equal(.ctl$rxControl$covsInterpolation,
               rxode2::rxControl(covsInterpolation = "nocb")$covsInterpolation)
  expect_false(.ctl$rxControl$addlKeepsCov)
})

test_that("as.nlmixr2() matches NONMEM for TIME in $PK (#252)", {
  skip_on_cran()
  skip_if_not(.nonmem2rxHasNonmemSolve(), "rxode2/nlmixr2est without the $PK record time")
  skip_if_not(file.exists(test_path("nonmem-pktime.zip")))
  .zip <- normalizePath(test_path("nonmem-pktime.zip"))
  withr::with_tempdir({
    unzip(.zip)
    # nonmem2rx stress kit case `time-in-pk` run with NONMEM 7.4:
    # CL = THETA(1)*(1 + THETA(4)*(1 - EXP(-THETA(5)*TIME)))*EXP(ETA(1))
    mod <- .nonmem2rx("nonmem-pktime/run.ctl", lst = ".lst", save = FALSE,
                      nonmemData = TRUE)
    # drop the untranslated `IF (W .EQ. 0) W = 1`
    .fun <- deparse(mod$fun)
    .fun <- .fun[!grepl("W == 0|W <- 1", .fun)]
    mod2 <- eval(str2lang(paste(.fun, collapse = "\n")))
    mod <- .as.nonmem2rx(mod2, mod)
    fit <- suppressWarnings(suppressMessages(as.nlmixr2(mod)))
    expect_true(inherits(fit, "nlmixr2FitData"))
    expect_true(fit$control$rxControl$nonmem)
    assign("nonmemControl", list(ci = 0.95), fit$env)
    .rel <- .nonmemMergePredsAndCalcRelativeErr(fit)$individualRel
    # with the continuous time CL is off by up to 21.8%
    expect_lt(.rel[["100%"]], 0.01)
  })
})
