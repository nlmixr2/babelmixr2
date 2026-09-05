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

test_that("the structural model emits derivative and assignment variables", {
  .x <- .pharmmlStructuralModel(.pharmmlTestUiOneCmt())

  expect_match(.x, 'blkId="sm1"')
  expect_match(.x, '<ct:DerivativeVariable symbId="depot" symbolType="real">')
  expect_match(.x, '<ct:DerivativeVariable symbId="center" symbolType="real">')
  expect_match(.x, '<ct:Variable symbId="cp" symbolType="real">')

  # every derivative declares its independent variable and initial condition
  expect_match(.x, "<ct:IndependentVariable>")
  expect_match(.x, "<ct:InitialCondition>")
  expect_match(.x, "<ct:InitialValue>")
})

test_that("the structural model does not emit the error-model helper lines", {
  .x <- .pharmmlStructuralModel(.pharmmlTestUiOneCmt())
  expect_false(grepl("rx_pred_", .x, fixed = TRUE))
})

test_that("the structural model qualifies parameter references with their block", {
  .x <- .pharmmlStructuralModel(.pharmmlTestUiOneCmt())
  # ka/cl/vc are individual parameters, so they live in pm1
  expect_match(.x, '<ct:SymbRef blkIdRef="pm1" symbIdRef="ka"/>')
  expect_match(.x, '<ct:SymbRef blkIdRef="pm1" symbIdRef="cl"/>')
  # states are local to the structural model, so they carry no blkIdRef
  expect_match(.x, '<ct:SymbRef symbIdRef="depot"/>')
})

test_that("the structural model honours explicit initial conditions", {
  .f <- function() {
    ini({ tcl <- log(2.72); eta.cl ~ 0.3; add.sd <- 0.7 })
    model({
      cl <- exp(tcl + eta.cl)
      center(0) <- 100
      d/dt(center) <- -cl * center
      cp <- center
      cp ~ add(add.sd)
    })
  }
  .x <- .pharmmlStructuralModel(rxode2::rxUiDecompress(.f()))
  expect_match(.x, "<ct:Real>100</ct:Real>")
})

test_that("the structural model is schema-valid", {
  for (.ui in list(.pharmmlTestUiOneCmt(), .pharmmlTestUiCorr())) {
    .doc <- .pharmmlWrapMdef(paste(.pharmmlVariabilityModel(.ui),
                                   .pharmmlParameterModel(.ui),
                                   .pharmmlStructuralModel(.ui), sep = "\n"))
    expect_true(pharmmlValidate(.doc))
  }
})

# Build a one-endpoint model with the given residual line.  `iniExtra` must
# declare exactly the error parameters the line uses -- nlmixr2 rejects an ini
# entry that the model body never references.
.pharmmlTestUiErr <- function(errLine, iniExtra) {
  .f <- eval(bquote(function() {
    ini({
      tcl <- log(2.72)
      eta.cl ~ 0.3
      .(iniExtra)
    })
    model({
      cl <- exp(tcl + eta.cl)
      d/dt(center) <- -cl * center
      cp <- center
      .(errLine)
    })
  }))
  rxode2::rxUiDecompress(.f())
}

.pharmmlErrAdd  <- quote(add.sd <- 0.7)
.pharmmlErrProp <- quote(prop.sd <- 0.1)
.pharmmlErrBoth <- quote({ add.sd <- 0.7; prop.sd <- 0.1 })

test_that("the observation model wires output, error model and residual", {
  .x <- .pharmmlObservationModel(.pharmmlTestUiOneCmt())

  expect_match(.x, 'blkId="om1"')
  expect_match(.x, "<mdef:ContinuousData>")
  expect_match(.x, '<mdef:PopulationParameter symbId="add.sd"/>')
  expect_match(.x, '<ct:SymbRef blkIdRef="vm2" symbIdRef="residual"/>')
  expect_match(.x, '<po:ProbOnto name="StandardNormal1"/>')
  expect_match(.x, '<ct:SymbRef blkIdRef="sm1" symbIdRef="cp"/>')
  expect_match(.x, "<mdef:ErrorModel>")
  expect_match(.x, "<mdef:ResidualError>")
})

test_that("the additive error model is a bare constant", {
  .x <- .pharmmlObservationModel(.pharmmlTestUiErr(quote(cp ~ add(add.sd)), .pharmmlErrAdd))
  expect_match(.x, '<ct:SymbRef symbIdRef="add.sd"/>')
  expect_false(grepl('op="times"', .x, fixed = TRUE))
})

test_that("the proportional error model multiplies by the prediction", {
  .x <- .pharmmlObservationModel(.pharmmlTestUiErr(quote(cp ~ prop(prop.sd)), .pharmmlErrProp))
  expect_match(.x, 'math:Binop op="times"')
  expect_match(.x, '<ct:SymbRef symbIdRef="prop.sd"/>')
})

test_that("combined1 is a + b*f and combined2 is sqrt(a^2 + (b*f)^2)", {
  .ui <- .pharmmlTestUiErr(quote(cp ~ add(add.sd) + prop(prop.sd)), .pharmmlErrBoth)

  rxode2::rxAssignControlValue(.ui, "addProp", "combined1")
  .x1 <- .pharmmlObservationModel(.ui)
  expect_match(.x1, 'math:Binop op="plus"')
  expect_false(grepl('Uniop op="sqrt"', .x1, fixed = TRUE))

  rxode2::rxAssignControlValue(.ui, "addProp", "combined2")
  .x2 <- .pharmmlObservationModel(.ui)
  expect_match(.x2, 'Uniop op="sqrt"')
  expect_match(.x2, 'math:Binop op="power"')
})

test_that("unsupported residual distributions error by name", {
  expect_error(.pharmmlObservationModel(
    .pharmmlTestUiErr(quote(cp ~ pow(prop.sd, add.sd)), .pharmmlErrBoth)), "pow")
})

test_that("the observation model is schema-valid", {
  .cases <- list(list(quote(cp ~ add(add.sd)), .pharmmlErrAdd),
                 list(quote(cp ~ prop(prop.sd)), .pharmmlErrProp),
                 list(quote(cp ~ add(add.sd) + prop(prop.sd)), .pharmmlErrBoth))
  for (.case in .cases) {
    .e <- .case[[1]]
    .ui <- .pharmmlTestUiErr(.e, .case[[2]])
    .doc <- .pharmmlWrapMdef(paste(.pharmmlVariabilityModel(.ui),
                                   .pharmmlParameterModel(.ui),
                                   .pharmmlStructuralModel(.ui),
                                   .pharmmlObservationModel(.ui), sep = "\n"))
    expect_true(pharmmlValidate(.doc), info = deparse1(.e))
  }
})
