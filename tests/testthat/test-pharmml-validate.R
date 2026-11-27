test_that("pharmmlValidate accepts a valid 0.9 document", {
  .f <- testthat::test_path("pharmml", "minimal-0.9.xml")
  expect_true(pharmmlValidate(.f))
})

test_that("pharmmlValidate rejects an invalid 0.9 document and reports why", {
  .f <- withr::local_tempfile(fileext = ".xml")
  .x <- readLines(testthat::test_path("pharmml", "minimal-0.9.xml"), warn = FALSE)
  .x <- sub("<IndependentVariable ", "<IndependentVariableTYPO ", .x, fixed = TRUE)
  writeLines(.x, .f)

  .v <- pharmmlValidate(.f, error = FALSE)
  expect_false(.v)
  expect_match(paste(attr(.v, "errors"), collapse = " "), "IndependentVariableTYPO")
})

test_that("pharmmlValidate errors by default on an invalid document", {
  .f <- withr::local_tempfile(fileext = ".xml")
  writeLines("<PharmML xmlns='http://www.pharmml.org/pharmml/0.9/PharmML'/>", .f)
  expect_error(pharmmlValidate(.f), "does not validate")
})

test_that("pharmmlValidate rejects an unsupported version", {
  .f <- testthat::test_path("pharmml", "minimal-0.9.xml")
  expect_error(pharmmlValidate(.f, version = "0.7"),
               "unsupported PharmML version")
})

test_that("pharmmlValidate accepts a character document, not just a path", {
  .x <- paste(readLines(testthat::test_path("pharmml", "minimal-0.9.xml"),
                        warn = FALSE), collapse = "\n")
  expect_true(pharmmlValidate(.x))
})

test_that("pharmmlVersions lists the supported versions", {
  expect_true("0.9" %in% pharmmlVersions())
})
