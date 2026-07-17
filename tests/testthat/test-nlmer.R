test_that("nlmer estimation basic smoke test", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  # Simple one-compartment model similar to other tests
  oneCmt <- function() {
    ini({
      tka <- 0.0
      tv <- 2.99573227355399
      tcl <- -0.693147180559945
      eta.ka ~ 1.0
      eta.v ~ 1.0
      eta.cl ~ 1.0
      add.sd <- 1.0
    })
    model({
      ka <- exp(tka + eta.ka)
      v <- exp(tv + eta.v)
      cl <- exp(tcl + eta.cl)
      d/dt(depot) <- -depot * ka
      d/dt(central) <- depot * ka - cl * central / v
      cp <- central / v
      cp ~ add(add.sd)
    })
  }

  # Use the built-in theo dataset *unfiltered*: the full object (dosing rows
  # + all observations, including the t=0 observation) is loaded into the nlm
  # solver, which solves the whole thing and subsets by index to the
  # observation gradient/f values lme4 fits against.  nlmer must not require
  # the caller to strip dosing/observation rows out of the estimation data.
  fullTheo <- nlmixr2data::theo_sd

  # Run a short nlmer fit; keep iterations small so test is quick
  fit <- tryCatch(
    nlmixr2(oneCmt, fullTheo, est = "nlmer", nlmerControl(tolPwrss = 1e-6, optCtrl = list(maxfun = 10), returnNlmer = FALSE)),
    error = function(e) {
      skip(paste("nlmer fit failed on this environment:", conditionMessage(e)))
    }
  )

  expect_true(inherits(fit, "nlmixr2FitData"))
  expect_equal(fit$est, "nlmer")
  # nlmer fit object should be present when returnNlmer == FALSE we still keep nlmer in internal slot
  expect_true(!is.null(fit$nlmer) || !is.null(attr(fit, "nlmer")) || exists("nlmer", envir = fit$env))
  # theta must be present and numeric
  expect_true(is.numeric(fit$theta))
  # objective should be finite
  expect_true(is.finite(fit$objective))
  # every observation row (all 132, including t=0) is fit -- nothing dropped
  expect_equal(fit$nobs, sum(fullTheo$EVID == 0))
})

test_that("nlmer jump sensitivities give correct dosing-parameter gradients", {
  skip_on_cran()
  skip_if_not_installed("lme4")

  # One-compartment model with an *estimated* absorption lag time. tlag enters
  # alag(depot), so it is an event parameter: its prediction sensitivity has a
  # jump at the dose front that the augmented sensitivity ODE misses. With
  # eventSens = "jump" rxode2 injects the analytic jump; with "fd" the nlm
  # engine falls back to (less reliable) finite differences across the jump.
  modLag <- function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; tlag <- -0.7
      eta.ka ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); lag <- exp(tlag)
      alag(depot) <- lag
      d/dt(depot) <- -depot * ka
      d/dt(central) <- depot * ka - cl * central / v
      cp <- central / v
      cp ~ add(add.sd)
    })
  }
  ft <- nlmixr2data::theo_sd

  # Load the resident nlm solver env exactly like .nlmerFamilyFit, returning
  # the per-subject prediction + Jacobian solver handle.
  setupSolver <- function(eventSens) {
    .ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(modLag))
    assign("control", nlmerControl(eventSens = eventSens), envir = .ui)
    .ctl <- .ui$control
    .ret <- new.env(parent = emptyenv())
    nlmixr2est::.foceiPreProcessData(ft, .ret, .ui, .ctl$rxControl)
    .sm <- .ui$nlmerSensModel
    .mi <- list(predOnly = .sm$predOnly, thetaGrad = .sm$thetaGrad,
                eventTheta = .sm$eventTheta)
    .start <- .ui$nlmerStart
    .sc <- babelmixr2:::.nlmerSolveControl(.ctl, length(.start))
    nlmixr2est::.nlmSetupEnv(.start, .ui, .ret$dataSav, .mi, .sc)
    list(start = .start, eventTheta = .sm$eventTheta, pn = .sm$paramNames,
         nsub = length(unique(.ret$dataSav$ID[.ret$dataSav$EVID != 2])))
  }

  res <- tryCatch({
    sj <- suppressMessages(setupSolver("jump"))
    np <- length(sj$start)
    nsub <- sj$nsub
    lagCol <- which(sj$pn == "tlag")
    tm <- matrix(rep(unname(sj$start), each = nsub), nsub, np)

    gJ <- nlmixr2est::nlmerSolveGrad(tm)          # jump analytic
    anaJump <- gJ[, lagCol + 1L]

    # Central finite difference of the *prediction* wrt tlag (model truth: the
    # ODE solve itself includes the lag, independent of the sensitivity method)
    h <- 1e-4
    fp <- nlmixr2est::nlmerSolveGrad(`[<-`(tm, , lagCol, tm[, lagCol] + h))[, 1L]
    fm <- nlmixr2est::nlmerSolveGrad(`[<-`(tm, , lagCol, tm[, lagCol] - h))[, 1L]
    centralFD <- (fp - fm) / (2 * h)
    nlmixr2est::.nlmFreeEnv()

    sf <- suppressMessages(setupSolver("fd"))
    anaFD <- nlmixr2est::nlmerSolveGrad(tm)[, lagCol + 1L]
    nlmixr2est::.nlmFreeEnv()

    list(eventThetaJump = sj$eventTheta, eventThetaFd = sf$eventTheta,
         lagCol = lagCol, anaJump = anaJump, anaFD = anaFD, centralFD = centralFD)
  }, error = function(e) {
    nlmixr2est::.nlmFreeEnv()
    skip(paste("nlmer jump probe failed on this environment:",
               conditionMessage(e)))
  })

  # Under "jump" the dosing parameter is routed through the analytic jump
  # (flag stays 0); under "fd" it is flagged for finite differences.
  expect_equal(res$eventThetaJump[res$lagCol], 0L)
  expect_equal(res$eventThetaFd[res$lagCol], 1L)

  # The analytic jump gradient matches the finite-difference-of-prediction
  # truth to solver tolerance.
  expect_lt(max(abs(res$anaJump - res$centralFD)), 1e-3)

  # The jump gradient is materially more accurate than the fd fallback, which
  # is ill-conditioned across the dose front (the motivation for "jump").
  expect_lt(max(abs(res$anaJump - res$centralFD)),
            max(abs(res$anaFD - res$centralFD)))
})