# est = "poped" walks each subject's records itself, so it has to advance
# newind the way rxSolve() does: 1 on a subject's first record, 2 after it.
# A model that reads newind is otherwise designed as a different model.

newindEvents <- function() {
  rxode2::et(amt = 100) |>
    rxode2::et(c(0.25, 1.2, 2, 4, 10, 20)) |>
    rxode2::et(id = 1:2) |>
    as.data.frame()
}

# six times per group: maxn = 15 takes the modeled-times walk, maxn = 3 the
# full walk
newindMaxn <- c(modeledTimes = 15L, full = 3L)

newindProbe <- function() {
  ini({
    tcl <- log(0.5)
    tv <- log(10)
    eta.cl ~ 0.1
    add.err <- 0.1
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    d/dt(central) <- -cl / v * central
    cp <- newind + 0 * central
    cp ~ add(add.err)
  })
}

newindE0 <- function() {
  ini({
    tcl <- log(0.5)
    tv <- log(10)
    eta.cl ~ 0.1
    add.err <- 0.1
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    d/dt(central) <- -cl / v * central
    e0 <- 1
    if (newind > 1) e0 <- 1.2
    cp <- central / v * e0
    cp ~ add(add.err)
  })
}

newindConstant <- function() {
  ini({
    tcl <- log(0.5)
    tv <- log(10)
    eta.cl ~ 0.1
    add.err <- 0.1
  })
  model({
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    d/dt(central) <- -cl / v * central
    cp <- central / v * 1.2
    cp ~ add(add.err)
  })
}

newindDesign <- function(ui, maxn) {
  suppressMessages(
    nlmixr2est::nlmixr2(ui, newindEvents(), "poped", popedControl(maxn = maxn, groupsize = 20))
  )
}

test_that("est='poped' predicts newind as rxSolve() does, on both walks", {
  skip_on_cran()
  skip_if_not_installed("PopED")
  solved <- as.data.frame(rxode2::rxSolve(newindProbe, newindEvents()))
  # every observation follows its subject's dose record
  expect_equal(solved$cp, rep(2, 12))
  for (walk in names(newindMaxn)) {
    pred <- PopED::model_prediction(newindDesign(newindProbe, newindMaxn[[walk]]))
    expect_equal(pred$PRED, solved$cp, info = walk)
  }
})

test_that("est='poped' designs a model that reads newind as the model it is", {
  skip_on_cran()
  skip_if_not_installed("PopED")
  for (walk in names(newindMaxn)) {
    # after the dose record e0 is 1.2 at every observation
    expect_equal(
      PopED::evaluate_design(newindDesign(newindE0, newindMaxn[[walk]]))$ofv,
      PopED::evaluate_design(newindDesign(newindConstant, newindMaxn[[walk]]))$ofv,
      tolerance = 1e-6,
      info = walk
    )
  }
})
