# Edge cases of the nlmixr2 -> NONMEM and nlmixr2 -> Monolix translations
# (translation only; see inst/stress/README.md to run NONMEM/Monolix end
# to end with the same cases)
test_that("NONMEM/Monolix stress test cases translate or are refused", {
  skip_on_cran()
  skip_if_not("linCmtMicro" %in% getNamespaceExports("rxode2"),
              "rxode2 does not have linCmtMicro() for closed-form linCmt()")
  source(system.file("stress", "stress.R", package="babelmixr2"), local=TRUE)
  withr::with_tempdir({
    res <- stressRun(stressCases(), mode="translate", dir=getwd())
  })
  failed <- res[stressFailed(res), ]
  expect_equal(nrow(failed), 0L,
               info=paste(sprintf("[%s] %s: %s %s", failed$engine, failed$case,
                                  failed$message, failed$problems),
                          collapse="\n"))
  # every case ran for both engines
  expect_equal(sort(unique(res$engine)), c("monolix", "nonmem"))
  expect_true(all(res$status %in% c("ok", "refused")))
})

test_that("stress run mode checks the fit, the predictions and a rerun", {
  skip_on_cran()
  .e <- new.env()
  source(system.file("stress", "stress.R", package = "babelmixr2"), local = .e)
  .case <- .e$stressCases()[["rerun reads saved output"]]
  .translate <- .e$.stressFit
  # write the real files, but return a made up fit (no NONMEM here)
  .e$.fits <- list()
  .e$.stressFit <- function(case, engine, mode, dir, modelName, runCommand) {
    .translate(case, engine, "translate", dir, modelName, runCommand)
    .ret <- .e$.fits[[1]]
    .e$.fits <- .e$.fits[-1]
    .ret
  }
  .fit <- function(objf, theta = c(tka = 0.45)) {
    structure(list(objf = objf, theta = theta), class = "nlmixr2FitData")
  }
  .e$.ipred <- 0.1
  .e$.stressPredDiff <- function(fit, engine) c(.e$.ipred, 0.2)
  .run <- function() {
    withr::with_tempdir(.e$stressRunCase(
      .case,
      "nonmem",
      mode = "run",
      dir = getwd()
    ))
  }

  .e$.fits <- list(.fit(100), .fit(100))
  .r <- .run()
  expect_equal(.r$status, "ok", info = .r$problems)
  expect_equal(.r$objf, 100)
  expect_equal(.r$ipredRelDiff, 0.1)
  expect_equal(.r$predRelDiff, 0.2)
  expect_false(is.na(.r$rerunSeconds))

  .e$.fits <- list(.fit(100), .fit(101))
  .r <- .run()
  expect_match(.r$problems, "rerun gave a different objective function")

  .e$.fits <- list(.fit(100), simpleError("no output"))
  .r <- .run()
  expect_match(.r$problems, "rerun failed: no output")

  .e$.fits <- list(.fit(NA_real_, c(tka = Inf)), .fit(NA_real_, c(tka = Inf)))
  .r <- .run()
  expect_equal(.r$status, "problem")
  expect_match(.r$problems, "objective function is not finite")
  expect_match(.r$problems, "non-finite estimates: tka")

  .e$.ipred <- 12
  .e$.fits <- list(.fit(100), .fit(100))
  .r <- .run()
  expect_match(.r$problems, "IPRED differs from nonmem by 12.00%")
  .e$.ipred <- NA_real_
  .e$.fits <- list(.fit(100), .fit(100))
  .r <- .run()
  expect_match(.r$problems, "cannot compare the predictions")
})

test_that("engine specific stress cases only run for their engine", {
  skip_on_cran()
  .e <- new.env()
  source(system.file("stress", "stress.R", package = "babelmixr2"), local = .e)
  .c <- .e$stressCases()
  .c <- .c[c("NONMEM est=imp", "Monolix stiff ODEs")]
  .res <- withr::with_tempdir(.e$stressRun(
    .c,
    mode = "translate",
    dir = getwd()
  ))
  expect_equal(.res$engine, c("nonmem", "monolix"))
  expect_equal(
    .res$status,
    c("ok", "ok"),
    info = paste(.res$message, .res$problems)
  )
})
