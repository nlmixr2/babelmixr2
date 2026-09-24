.pharmmlTestData <- function() {
  .d <- nlmixr2data::theo_sd
  .d$logWT70 <- log(.d$WT / 70)
  .d$SEX <- ifelse(.d$ID %% 2 == 0, "M", "F")
  .d
}

# Model with one continuous and one categorical covariate.
.pharmmlTestUiCovBoth <- function() {
  .f <- function() {
    ini({
      tka <- log(1.57)
      tcl <- log(2.72)
      tv <- log(31.5)
      wt.cl <- 0.75
      sex.cl <- 0.1
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + wt.cl * logWT70 + sex.cl * SEX + eta.cl)
      vc <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / vc * center
      cp <- center / vc
      cp ~ add(add.sd)
    })
  }
  rxode2::rxUiDecompress(.f())
}

test_that("covariate classification separates continuous from categorical", {
  .info <- .pharmmlCovariateInfo(.pharmmlTestUiCovBoth(), .pharmmlTestData())

  expect_equal(.info$logWT70$type, "continuous")
  expect_equal(.info$SEX$type, "categorical")
  # levels come from the original data, alphabetically as factor() orders them
  expect_equal(.info$SEX$levels, c("F", "M"))
  expect_equal(.info$SEX$codes, c(1L, 2L))
})

test_that("covariate classification treats numeric covariates as continuous", {
  .d <- .pharmmlTestData()
  .d$SEX <- ifelse(.d$SEX == "M", 1, 0) # already numerically coded
  .info <- .pharmmlCovariateInfo(.pharmmlTestUiCovBoth(), .d)
  expect_equal(.info$SEX$type, "continuous")
})

test_that("covariate classification without data leaves everything continuous", {
  .info <- .pharmmlCovariateInfo(.pharmmlTestUiCovBoth(), NULL)
  expect_equal(.info$logWT70$type, "continuous")
  expect_equal(.info$SEX$type, "continuous")
})

test_that("the covariate model emits Categorical with a Category per level", {
  .x <- .pharmmlCovariateModel(.pharmmlTestUiCovBoth(), .pharmmlTestData())

  expect_match(.x, '<mdef:Covariate symbId="logWT70">')
  expect_match(.x, "<mdef:Continuous/>")
  expect_match(.x, '<mdef:Covariate symbId="SEX">')
  expect_match(.x, "<mdef:Categorical>")
  expect_match(.x, '<mdef:Category catId="F"/>')
  expect_match(.x, '<mdef:Category catId="M"/>')
})

test_that("the covariate model still emits Continuous when given no data", {
  .x <- .pharmmlCovariateModel(.pharmmlTestUiCovBoth(), NULL)
  expect_match(.x, "<mdef:Continuous/>")
  expect_false(grepl("<mdef:Categorical>", .x, fixed = TRUE))
})

# --- trial design ---------------------------------------------------------

test_that("the trial design maps each dataset column to its model symbol", {
  .ui <- .pharmmlTestUiCovBoth()
  .x <- .pharmmlTrialDesign(.ui, .pharmmlTestData(), dataFile = "theo.csv")

  expect_match(.x, '<design:ExternalDataSet toolName="NONMEM"')
  # TIME maps to the independent variable
  expect_match(.x, '<ds:ColumnRef columnIdRef="TIME"/>')
  expect_match(.x, '<ct:SymbRef symbIdRef="t"/>')
  # a covariate maps into the covariate block, under its *model* name
  expect_match(.x, '<ds:ColumnRef columnIdRef="LOGWT70"/>')
  expect_match(.x, '<ct:SymbRef blkIdRef="cm1" symbIdRef="logWT70"/>')
  # the observation column maps to the observation model
  expect_match(.x, '<ds:ColumnRef columnIdRef="DV"/>')
  expect_match(.x, '<ct:SymbRef blkIdRef="om1" symbIdRef="cp_obs"/>')
})

test_that("the trial design declares every column with a type and a number", {
  .x <- .pharmmlTrialDesign(
    .pharmmlTestUiCovBoth(),
    .pharmmlTestData(),
    dataFile = "theo.csv"
  )
  expect_match(
    .x,
    '<ds:Column columnId="ID" columnType="id" valueType="int" columnNum="1"/>'
  )
  expect_match(.x, 'columnId="TIME" columnType="idv" valueType="real"')
  expect_match(.x, 'columnId="DV" columnType="dv" valueType="real"')
  expect_match(.x, 'columnId="AMT" columnType="dose" valueType="real"')
  expect_match(.x, 'columnId="EVID" columnType="evid" valueType="int"')
  expect_match(.x, 'columnId="CMT" columnType="cmt" valueType="int"')
  expect_match(.x, 'columnId="LOGWT70" columnType="covariate" valueType="real"')
})

test_that("the trial design does not export the internal row-number column", {
  .x <- .pharmmlTrialDesign(
    .pharmmlTestUiCovBoth(),
    .pharmmlTestData(),
    dataFile = "theo.csv"
  )
  expect_false(grepl("nlmixrRowNums", .x, fixed = TRUE))
})

test_that("a categorical covariate column carries a CategoricalMapping", {
  .x <- .pharmmlTrialDesign(
    .pharmmlTestUiCovBoth(),
    .pharmmlTestData(),
    dataFile = "theo.csv"
  )
  expect_match(.x, "<ds:CategoryMapping>")
  expect_match(.x, '<ds:Map dataSymbol="1" modelSymbol="F"/>')
  expect_match(.x, '<ds:Map dataSymbol="2" modelSymbol="M"/>')
})

test_that("the trial design names the external data file", {
  .x <- .pharmmlTrialDesign(
    .pharmmlTestUiCovBoth(),
    .pharmmlTestData(),
    dataFile = "theo.csv"
  )
  expect_match(.x, "<ds:path>theo.csv</ds:path>")
  expect_match(.x, "<ds:format>CSV</ds:format>")
  expect_match(.x, "<ds:delimiter>COMMA</ds:delimiter>")
})

# --- modelling steps ------------------------------------------------------

test_that("the estimation step lists every estimated parameter", {
  .x <- .pharmmlModellingSteps(.pharmmlTestUiCovBoth())

  expect_match(.x, "<mstep:EstimationStep")
  expect_match(.x, "<mstep:ParametersToEstimate>")
  for (.p in c(
    "tka",
    "tcl",
    "tv",
    "wt.cl",
    "sex.cl",
    "add.sd",
    "omega_eta.ka",
    "omega_eta.cl",
    "omega_eta.v"
  )) {
    expect_match(.x, paste0('symbIdRef="', .p, '"'), info = .p)
  }
})

test_that("the estimation step carries initial values and bounds", {
  .x <- .pharmmlModellingSteps(.pharmmlTestUiCovBoth())
  expect_match(.x, "<mstep:InitialEstimate>")
  # add.sd is bounded below at 0 in the iniDf
  expect_match(.x, "<mstep:LowerBound>")
})

test_that("a fixed parameter is declared fixed rather than estimated", {
  .f <- function() {
    ini({
      tcl <- fix(log(2.72))
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      cl <- exp(tcl + eta.cl)
      d / dt(center) <- -cl * center
      cp <- center
      cp ~ add(add.sd)
    })
  }
  .x <- .pharmmlModellingSteps(rxode2::rxUiDecompress(.f()))
  expect_match(.x, 'fixed="true"')
})

test_that("the estimation step names the block each parameter is declared in", {
  .x <- .pharmmlModellingSteps(.pharmmlTestUiCovBoth())
  # residual-error parameters live in the ObservationModel, not the
  # ParameterModel
  expect_match(.x, '<ct:SymbRef blkIdRef="om1" symbIdRef="add.sd"/>')
  expect_match(.x, '<ct:SymbRef blkIdRef="pm1" symbIdRef="tcl"/>')
  expect_equal(
    .pharmmlDanglingRefs(as.pharmml(
      .pharmmlTestUiCovBoth(),
      .pharmmlTestData()
    )),
    character(0)
  )
})

# --- whole document -------------------------------------------------------

test_that("a complete PharmML document is schema-valid", {
  .ui <- .pharmmlTestUiCovBoth()
  .doc <- paste0(
    '<?xml version="1.0" encoding="UTF-8"?>\n',
    '<PharmML xmlns="http://www.pharmml.org/pharmml/0.9/PharmML"\n',
    '    xmlns:ct="http://www.pharmml.org/pharmml/0.9/CommonTypes"\n',
    '    xmlns:math="http://www.pharmml.org/pharmml/0.9/Maths"\n',
    '    xmlns:mdef="http://www.pharmml.org/pharmml/0.9/ModelDefinition"\n',
    '    xmlns:design="http://www.pharmml.org/pharmml/0.9/TrialDesign"\n',
    '    xmlns:mstep="http://www.pharmml.org/pharmml/0.9/ModellingSteps"\n',
    '    xmlns:ds="http://www.pharmml.org/pharmml/0.9/Dataset"\n',
    '    xmlns:po="http://www.pharmml.org/probonto/ProbOnto"\n',
    '    writtenVersion="0.9" id="i1">\n',
    '  <ct:Name>full document</ct:Name>\n',
    '  <IndependentVariable symbId="t"/>\n',
    .pharmmlModelDefinition(.ui, data = .pharmmlTestData()),
    "\n",
    .pharmmlTrialDesign(.ui, .pharmmlTestData(), dataFile = "theo.csv"),
    "\n",
    .pharmmlModellingSteps(.ui),
    "\n",
    '</PharmML>\n'
  )
  expect_true(pharmmlValidate(.doc))
})

test_that("a category-code mismatch is an error, not a silent mis-mapping", {
  # Guard the derived level -> code mapping.  If bblDatToNonmem() ever coded
  # categories differently, the mapping would be transposed rather than wrong
  # in an obvious way, so the check must actually fire.
  .d <- .pharmmlTestData()
  .nmBad <- .pharmmlNonmemData(.pharmmlTestUiCovBoth(), .d)
  .nmBad$SEX <- 3 - .nmBad$SEX # swap the two codes
  expect_error(
    .pharmmlCovariateInfo(.pharmmlTestUiCovBoth(), .d, .nmBad),
    "do not match the levels derived"
  )
})

test_that("the category-code check passes on the real conversion", {
  .d <- .pharmmlTestData()
  .nm <- .pharmmlNonmemData(.pharmmlTestUiCovBoth(), .d)
  .info <- .pharmmlCovariateInfo(.pharmmlTestUiCovBoth(), .d, .nm)
  expect_equal(.info$SEX$levels, c("F", "M"))
})
