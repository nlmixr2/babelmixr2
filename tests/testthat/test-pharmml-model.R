# Shared fixtures for the ModelDefinition writers.

.pharmmlTestUiOneCmt <- function() {
  .f <- function() {
    ini({
      tka <- log(1.57); label("Ka")
      tcl <- log(2.72); label("Cl")
      tv  <- log(31.5); label("V")
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v  ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      vc <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / vc * center
      cp <- center / vc
      cp ~ add(add.sd)
    })
  }
  rxode2::rxUiDecompress(.f())
}

# Same model with a correlated eta block, to exercise <Correlation>.
.pharmmlTestUiCorr <- function() {
  .f <- function() {
    ini({
      tka <- log(1.57)
      tcl <- log(2.72)
      tv  <- log(31.5)
      eta.ka + eta.cl ~ c(0.6, 0.01, 0.3)
      eta.v  ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      vc <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / vc * center
      cp <- center / vc
      cp ~ add(add.sd)
    })
  }
  rxode2::rxUiDecompress(.f())
}

test_that("the variability model declares an id level and a residual level", {
  .x <- .pharmmlVariabilityModel(.pharmmlTestUiOneCmt())

  expect_match(.x, 'blkId="vm1"')
  expect_match(.x, 'type="parameterVariability"')
  expect_match(.x, 'referenceLevel="true"')
  expect_match(.x, 'symbId="id"')

  expect_match(.x, 'blkId="vm2"')
  expect_match(.x, 'type="residualError"')
  expect_match(.x, 'symbId="residual"')
})

test_that("the variability model omits vm1 when the model has no etas", {
  .f <- function() {
    ini({ tcl <- log(2.72); add.sd <- 0.7 })
    model({ cl <- exp(tcl); d/dt(center) <- -cl * center
            cp <- center; cp ~ add(add.sd) })
  }
  .x <- .pharmmlVariabilityModel(rxode2::rxUiDecompress(.f()))
  expect_false(grepl('type="parameterVariability"', .x, fixed = TRUE))
  expect_match(.x, 'type="residualError"')
})
