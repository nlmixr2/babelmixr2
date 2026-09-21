.pharmmlApiModel <- function() {
  one.cmt <- function() {
    ini({
      tka <- log(1.57); tcl <- log(2.72); tv <- log(31.5)
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
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
  one.cmt
}

test_that("pharmmlControl validates its arguments", {
  expect_s3_class(pharmmlControl(), "pharmmlControl")
  expect_true(pharmmlControl()$validate)
  expect_false(pharmmlControl(validate = FALSE)$validate)
  expect_error(pharmmlControl(version = "0.7"), "unsupported PharmML version")
  expect_error(pharmmlControl(validate = "yes"))
})

test_that("rxToPharmml translates an expression", {
  expect_match(rxToPharmml("ka * depot"), 'math:Binop op="times"')
  expect_match(rxToPharmml("exp(tcl)"), 'Uniop op="exp"')
})

test_that("rxToPharmml qualifies symbols when given a ui", {
  .ui <- rxode2::rxUiDecompress(.pharmmlApiModel()())
  expect_match(rxToPharmml("ka * depot", .ui), 'blkIdRef="pm1" symbIdRef="ka"')
  # without a ui there is nothing to resolve against
  expect_false(grepl("blkIdRef", rxToPharmml("ka * depot"), fixed = TRUE))
})

test_that("as.pharmml writes a schema-valid model-only document", {
  .x <- as.pharmml(.pharmmlApiModel())
  expect_s3_class(.x, "babelmixr2pharmml")
  expect_match(.x, "<PharmML")
  expect_match(.x, 'writtenVersion="0.9"')
  expect_match(.x, "<mdef:ModelDefinition>")
  # no data, so no trial design or estimation step.  Assert on the elements,
  # not the bare words -- the xmlns declarations mention both namespaces.
  expect_false(grepl("<design:TrialDesign>", .x, fixed = TRUE))
  expect_false(grepl("<mstep:ModellingSteps>", .x, fixed = TRUE))
  expect_true(pharmmlValidate(unclass(.x)))
})

test_that("as.pharmml with data writes the full document", {
  .x <- as.pharmml(.pharmmlApiModel(), nlmixr2data::theo_sd)
  expect_match(.x, "<design:TrialDesign>")
  expect_match(.x, "<mstep:ModellingSteps>")
  expect_true(pharmmlValidate(unclass(.x)))
})

test_that("as.pharmml names the document after the model", {
  # modelName depends on how the model reaches assertRxUi(), so take it from
  # the ui rather than hard-coding what the fixture happens to be called.
  .ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(.pharmmlApiModel()))
  .x <- as.pharmml(.pharmmlApiModel())
  expect_match(.x, paste0("<ct:Name>", .ui$modelName, "</ct:Name>"), fixed = TRUE)
  expect_match(.x, "translated to PharmML 0.9 by babelmixr2")
})

test_that("as.pharmml honours a supplied description", {
  .x <- as.pharmml(.pharmmlApiModel(),
                   control = pharmmlControl(description = "for the archive"))
  expect_match(.x, "<ct:Description>for the archive</ct:Description>")
})

test_that("as.pharmml writes the model and its dataset to disk", {
  .dir <- withr::local_tempdir()
  .f <- file.path(.dir, "theo.xml")
  as.pharmml(.pharmmlApiModel(), nlmixr2data::theo_sd, file = .f)

  .ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(.pharmmlApiModel()))
  .csvFile <- file.path(.dir, paste0(.ui$modelName, ".csv"))

  expect_true(file.exists(.f))
  expect_true(file.exists(.csvFile))
  # the written document is the one that was validated
  expect_true(pharmmlValidate(.f))

  .csv <- utils::read.csv(.csvFile)
  expect_true(all(c("ID", "TIME", "DV", "AMT", "EVID", "CMT") %in% names(.csv)))
  expect_false("nlmixrRowNums" %in% names(.csv))
})

test_that("as.pharmml can skip writing the dataset", {
  .dir <- withr::local_tempdir()
  .f <- file.path(.dir, "theo.xml")
  as.pharmml(.pharmmlApiModel(), nlmixr2data::theo_sd, file = .f,
             control = pharmmlControl(writeData = FALSE))
  .ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(.pharmmlApiModel()))
  expect_true(file.exists(.f))
  expect_false(file.exists(file.path(.dir, paste0(.ui$modelName, ".csv"))))
})

test_that("as.pharmml refuses a model PharmML cannot express", {
  # A power residual builds fine as an nlmixr2 model but has no PharmML
  # Standard error-model form, so it must be refused by name rather than
  # emitted as something plausible-looking.
  .f <- function() {
    ini({ tcl <- log(2.72); eta.cl ~ 0.3; prop.sd <- 0.1; pw <- 0.8 })
    model({ cl <- exp(tcl + eta.cl); d/dt(center) <- -cl * center
            cp <- center; cp ~ pow(prop.sd, pw) })
  }
  expect_error(as.pharmml(.f), "not supported")
})

test_that("as.pharmml rejects a control object it did not make", {
  expect_error(as.pharmml(.pharmmlApiModel(), control = list(version = "0.9")),
               "pharmmlControl")
})

test_that("the document id is a valid NCName even for an awkward model name", {
  expect_equal(.pharmmlId("one.cmt"), "one.cmt")
  expect_equal(.pharmmlId("2cmt model"), "i2cmt_model")
  expect_equal(.pharmmlId("a/b"), "a_b")
})

test_that("printing a document does not mangle it", {
  .x <- as.pharmml(.pharmmlApiModel())
  expect_output(print(.x), "<PharmML")
})
