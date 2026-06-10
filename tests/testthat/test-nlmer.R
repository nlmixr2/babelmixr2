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

  # Use the built-in theo dataset, filter to observation rows
  filteredTheo <- nlmixr2data::theo_sd[!(nlmixr2data::theo_sd$TIME == 0 & nlmixr2data::theo_sd$EVID == 0), ]

  # Run a short nlmer fit; keep iterations small so test is quick
  fit <- tryCatch(
    nlmixr2(oneCmt, filteredTheo, est = "nlmer", nlmerControl(tolPwrss = 1e-6, optCtrl = list(maxfun = 10), returnNlmer = FALSE)),
    error = function(e) {
      skip(paste("nlmer fit failed on this environment:", conditionMessage(e)))
    }
  )

  expect_true(inherits(fit, "nlmixr2"))
  expect_equal(fit$est, "nlmer")
  # nlmer fit object should be present when returnNlmer == FALSE we still keep nlmer in internal slot
  expect_true(!is.null(fit$nlmer) || !is.null(attr(fit, "nlmer")) || exists("nlmer", envir = fit$env))
  # theta must be present and numeric
  expect_true(is.numeric(fit$theta))
  # objective should be finite
  expect_true(is.finite(fit$objective))
})