# A model with a repeated (same()) omega block: nlmixr2est refuses it to an
# estimation method that does not say it keeps the tie.  pknca and poped take
# omega as given ("sameBlock" = "fixed"); the methods that estimate omega apart
# are still refused.

sameModel <- function() {
  ini({
    tvka <- 0.45
    lcl <- 1
    lvc <- 3.45
    prop.err <- 0.5
    eta.cl ~ 0.3
    eta.vc ~ same()
  })
  model({
    ka <- tvka
    cl <- exp(lcl + eta.cl)
    vc <- exp(lvc + eta.vc)
    cp <- linCmt()
    cp ~ prop(prop.err)
  })
}

refusal <- "does not honour a repeated (same()) omega block"

skipWithoutGuard <- function() {
  testthat::skip_if_not(
    exists(".nlmixr0preProcessAssertSameBlock", envir = asNamespace("nlmixr2est"), inherits = FALSE),
    "nlmixr2est does not refuse a repeated omega block"
  )
}

test_that("the methods that take omega as given declare it", {
  expect_identical(attr(nlmixr2Est.pknca, "sameBlock"), "fixed")
  expect_identical(attr(nlmixr2Est.poped, "sameBlock"), "fixed")
  # the methods that estimate omega say nothing
  for (.est in c("nonmem", "monolix", "saemix", "nlmer")) {
    expect_null(attr(get(paste0("nlmixr2Est.", .est)), "sameBlock"), info = .est)
  }
})

test_that("a model with a same() block reaches est='pknca'", {
  skipWithoutGuard()
  skip_if_not_installed("PKNCA")
  suppressMessages(expect_s3_class(
    nlmixr2est::nlmixr(object = sameModel, data = nlmixr2data::theo_sd, est = "pknca"),
    "pkncaEst"
  ))
})

test_that("a model with a same() block is not refused by est='poped'", {
  skipWithoutGuard()
  skip_if_not_installed("PopED")
  e <- rxode2::et(amt = 25)
  e <- rxode2::et(e, c(0.25, 1.2, 2, 4, 10, 20))
  e <- rxode2::et(e, id = 1:2)
  e <- as.data.frame(e)
  e$low <- NA_real_
  e$high <- NA_real_
  # what the design does with the block is the design's affair: it is not refused
  # before it gets there
  msg <- tryCatch(
    {
      suppressMessages(nlmixr2est::nlmixr2(sameModel, e, "poped", popedControl(maxn = 15)))
      ""
    },
    error = conditionMessage
  )
  expect_false(grepl(refusal, msg, fixed = TRUE))
})

test_that("a model with a same() block is still refused by the methods that estimate omega", {
  skipWithoutGuard()
  for (.est in c("nonmem", "monolix", "saemix", "nlmer")) {
    expect_error(
      suppressMessages(nlmixr2est::nlmixr(object = sameModel, data = nlmixr2data::theo_sd, est = .est)),
      paste0("'", .est, "' ", refusal),
      fixed = TRUE,
      info = .est
    )
  }
})
