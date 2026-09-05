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

# Wrap a ModelDefinition fragment in the smallest valid document.
.pharmmlWrapMdef <- function(x) {
  paste0(
    '<?xml version="1.0" encoding="UTF-8"?>\n',
    '<PharmML xmlns="http://www.pharmml.org/pharmml/0.9/PharmML"\n',
    '    xmlns:ct="http://www.pharmml.org/pharmml/0.9/CommonTypes"\n',
    '    xmlns:math="http://www.pharmml.org/pharmml/0.9/Maths"\n',
    '    xmlns:mdef="http://www.pharmml.org/pharmml/0.9/ModelDefinition"\n',
    '    xmlns:po="http://www.pharmml.org/probonto/ProbOnto"\n',
    '    writtenVersion="0.9" id="i1">\n',
    '  <ct:Name>fixture</ct:Name>\n',
    '  <IndependentVariable symbId="t"/>\n',
    '  <mdef:ModelDefinition>\n', x, '\n  </mdef:ModelDefinition>\n',
    '</PharmML>\n')
}

test_that("the parameter model declares thetas, omegas, etas and individual parameters", {
  .ui <- .pharmmlTestUiOneCmt()
  .x <- .pharmmlParameterModel(.ui)

  expect_match(.x, 'blkId="pm1"')

  # population parameters: one per theta, one per omega
  expect_match(.x, '<mdef:PopulationParameter symbId="tka"/>')
  expect_match(.x, '<mdef:PopulationParameter symbId="omega_eta.ka"/>')

  # random variable referencing the id variability level
  expect_match(.x, '<mdef:RandomVariable symbId="eta.ka">')
  expect_match(.x, '<ct:SymbRef blkIdRef="vm1" symbIdRef="id"/>')
  expect_match(.x, '<po:ProbOnto name="Normal2">')
  expect_match(.x, '<po:Parameter name="var">')

  # individual parameter, log-transformed because the model used exp()
  expect_match(.x, '<mdef:IndividualParameter symbId="ka">')
  expect_match(.x, '<mdef:Transformation type="log"/>')
})

test_that("the parameter model maps every mu-referenced parameter", {
  .x <- .pharmmlParameterModel(.pharmmlTestUiOneCmt())
  for (.v in c("ka", "cl", "vc")) {
    expect_match(.x, paste0('<mdef:IndividualParameter symbId="', .v, '">'), info = .v)
  }
})

test_that("the parameter model emits a correlation block for correlated etas", {
  .x <- .pharmmlParameterModel(.pharmmlTestUiCorr())
  expect_match(.x, "<mdef:Correlation>")
  expect_match(.x, 'symbIdRef="eta.ka"')
  expect_match(.x, 'symbIdRef="eta.cl"')
  expect_match(.x, "cov_eta.ka_eta.cl")
})

test_that("the parameter model omits the correlation block when etas are diagonal", {
  .x <- .pharmmlParameterModel(.pharmmlTestUiOneCmt())
  expect_false(grepl("<mdef:Correlation>", .x, fixed = TRUE))
})

test_that("the parameter model is schema-valid", {
  for (.ui in list(.pharmmlTestUiOneCmt(), .pharmmlTestUiCorr())) {
    .doc <- .pharmmlWrapMdef(paste(.pharmmlVariabilityModel(.ui),
                                   .pharmmlParameterModel(.ui), sep = "\n"))
    expect_true(pharmmlValidate(.doc))
  }
})
