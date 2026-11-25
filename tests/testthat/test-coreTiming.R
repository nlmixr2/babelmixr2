.hasCoreTiming <- function() {
  exists("coreTimingRegister", envir = asNamespace("nlmixr2extra"),
         mode = "function", inherits = FALSE)
}

test_that("babelmixr2 registers its methods with coreTiming()", {
  skip_if_not(.hasCoreTiming(), "nlmixr2extra without coreTimingRegister()")
  .m <- nlmixr2extra::coreTimingMethods()
  expect_true(all(c("nonmemControl", "monolixControl", "saemixControl",
                    "nlmerControl", "fmeMcmcControl", "pseudoOptimControl")
                  %in% .m$name))
  expect_equal(.m$args[.m$name == "nonmemControl"], "maxeval, niter")
  expect_true(.m$cores[.m$name == "nonmemControl"])
  expect_true(.m$cores[.m$name == "monolixControl"])
})

test_that("the NONMEM and Monolix truncations cap every phase", {
  skip_if_not(.hasCoreTiming(), "nlmixr2extra without coreTimingRegister()")
  .tr <- loadNamespace("nlmixr2extra")$.ctTruncateControl
  .ctl <- .tr(nonmemControl(), 3L, "nonmem")
  expect_equal(.ctl$maxeval, 3L)
  expect_equal(.ctl$niter, 3L)
  .ctl <- .tr(monolixControl(), 3L, "monolix")
  expect_equal(.ctl$exploratoryIterations, 3L)
  expect_equal(.ctl$smoothingIterations, 3L)
  expect_equal(.ctl$burnInIterations, 3L)
  expect_equal(.ctl$simulatedAnnealingIterations, 3L)
  expect_false(.ctl$exploratoryAutoStop)
  expect_false(.ctl$smoothingAutoStop)
  # the external engines warn that only R-side threads change
  expect_warning(loadNamespace("nlmixr2extra")$.ctSetCores(nonmemControl(), 2L, "nonmem"),
                 "runCommand")
})

test_that("coreTiming() times saemix and nlmer with their caps", {
  skip_on_cran()
  skip_if_not(.hasCoreTiming(), "nlmixr2extra without coreTimingRegister()")
  skip_if_not_installed("saemix")
  skip_if_not_installed("lme4")
  oneCmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  .ui <- rxode2::as.rxUi(oneCmt)
  .ct <- suppressWarnings(nlmixr2extra::coreTiming(
    .ui, data = nlmixr2data::theo_sd, est = "saemix", cores = 1L,
    estControl = saemixControl(print = 0),
    control = list(nIter = 3, reps = 1, warmup = FALSE)))
  expect_true(all(is.na(.ct$data$error)))
  expect_true("nbiter.saemix = 3, 3" %in% .ct$limits)
  .ct <- suppressWarnings(nlmixr2extra::coreTiming(
    .ui, data = nlmixr2data::theo_sd, est = "nlmer", cores = 1L,
    estControl = nlmerControl(print = 0),
    control = list(nIter = 3, reps = 1, warmup = FALSE)))
  expect_true(all(is.na(.ct$data$error)))
  expect_true("optCtrl$maxfun = 50" %in% .ct$limits)
})

test_that("coreTiming() times pseudoOptim and fmeMcmc with their caps", {
  skip_on_cran()
  skip_if_not(.hasCoreTiming(), "nlmixr2extra without coreTimingRegister()")
  skip_if_not_installed("FME")
  rxode2::rxWithSeed(42, {
    dsn <- data.frame(i = 1:200)
    dsn$time <- exp(rnorm(200))
    dsn$DV <- rbinom(200, 1, exp(-1 + dsn$time) / (1 + exp(-1 + dsn$time)))
  })
  mod <- function() {
    ini({
      E0 <- 0.5
      Em <- 0.5
      E50 <- 2
      g <- fix(2)
    })
    model({
      v <- E0 + Em * time^g / (E50^g + time^g)
      ll(bin) ~ DV * v - log(1 + exp(v))
    })
  }
  .ui <- rxode2::as.rxUi(mod)
  .ct <- suppressWarnings(nlmixr2extra::coreTiming(
    .ui, data = dsn, est = "fmeMcmc", cores = 1L,
    estControl = fmeMcmcControl(print = 0),
    control = list(nIter = 20, reps = 1, warmup = FALSE)))
  expect_true(all(is.na(.ct$data$error)))
  expect_true("niter = 20" %in% .ct$limits)
  # pseudoOptim needs finite bounds on every estimated parameter
  .ui <- rxode2::ini(.ui, E0 = c(-5, 0.5, 5), Em = c(-5, 0.5, 5), E50 = c(0.01, 2, 20))
  .ct <- suppressWarnings(nlmixr2extra::coreTiming(
    .ui, data = dsn, est = "pseudoOptim", cores = 1L,
    estControl = pseudoOptimControl(print = 0),
    control = list(nIter = 20, reps = 1, warmup = FALSE)))
  expect_true(all(is.na(.ct$data$error)))
  expect_true("numiter = 20" %in% .ct$limits)
})
