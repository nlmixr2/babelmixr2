## Monolix 2024R1 runs of the same theophylline model (made with the
## issue #207 check script), trimmed to the files babelmixr2 reads:
##
## - mle:          no prior
## - ka_weak:      prior(tka) ~ dnorm(log(3), 100), so MAP ~ the MLE
## - beta_strong:  prior(cl.wt) ~ dnorm(0.75, 0.001), so MAP ~ 0.75
## - logit_strong: prior(tfr) ~ dnorm(logit(0.6), 0.01), so fr ~ 0.6
## - two_priors:   tka ~ dnorm(log(3), 0.01) and tcl ~ dnorm(log(2), 0.01)
## - addsd_strong: prior(add.sd) ~ dnorm(1.5, 0.001), written by an earlier
##                 build; Monolix ignored it, which is why babelmixr2 now
##                 refuses priors on residual error parameters
##
## They check that babelmixr2 reads a 2024R1 export with a mu-referenced
## covariate, and that Monolix used the MAP prior the way nlmixr2 means it.

.mlx2024Model <- function() {
  ini({
    tka <- log(1.5)
    tcl <- log(2.7)
    tv <- log(32)
    tfr <- logit(0.9)
    cl.wt <- 0
    eta.ka ~ 0.4
    eta.cl ~ 0.1
    eta.v ~ 0.02
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl + cl.wt * lWT)
    v <- exp(tv + eta.v)
    fr <- expit(tfr)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- fr * center / v
    cp ~ add(add.sd)
  })
}

.mlx2024Data <- function() {
  .dat <- nlmixr2data::theo_sd
  .dat$lWT <- log(.dat$WT / 70)
  .dat
}

.mlx2024Ui <- function(case) {
  .u <- rxode2::rxode2(.mlx2024Model)
  if (case == "mle") return(.u)
  .p <- try(ini(.u, prior(tka) ~ dnorm(0, 1)), silent=TRUE)
  if (inherits(.p, "try-error") || !any(names(.p$iniDf) == "prior")) {
    skip("this lotri/rxode2 does not support priors")
  }
  switch(case,
         ka_weak=ini(.u, prior(tka) ~ dnorm(log(3), 100)),
         beta_strong=ini(.u, prior(cl.wt) ~ dnorm(0.75, 0.001)),
         logit_strong=ini(.u, prior(tfr) ~ dnorm(0.4054651, 0.01)),
         two_priors=ini(ini(.u, prior(tka) ~ dnorm(log(3), 0.01)),
                        prior(tcl) ~ dnorm(log(2), 0.01)))
}

.mlx2024Fit <- function(case) {
  .zip <- normalizePath(test_path("monolix2024-priors.zip"))
  .u <- .mlx2024Ui(case)
  withr::with_tempdir({
    utils::unzip(.zip)
    setwd(file.path("monolix2024", case))
    # the current translation still writes the prior lines Monolix 2024R1
    # accepted (the [POPULATION] definitions and which parameters are MAP)
    .prior <- function(l) {
      l[grepl("typical=[-0-9.e]+, sd=|method=MAP", l)]
    }
    expect_equal(.prior(strsplit(.u$mlxtran, "\n")[[1]]),
                 .prior(readLines(paste0(case, "-monolix.mlxtran"))))
    suppressWarnings(suppressMessages(
      nlmixr2(.u, .mlx2024Data(), "monolix", monolixControl(modelName=case))))
  })
}

test_that("a Monolix 2024R1 export with a mu-referenced covariate is read", {
  skip_if_not(file.exists(test_path("monolix2024-priors.zip")))
  .f <- .mlx2024Fit("mle")
  expect_true(inherits(.f, "nlmixr2FitData"))
  # the covariance used to look up 'NA_pop' for the covariate effect
  expect_true(all(c("beta_cl_lWT", "ka_pop") %in% dimnames(.f$cov)[[1]]) ||
                all(c("cl.wt", "tka") %in% dimnames(.f$cov)[[1]]))
  expect_equal(.f$theta[["cl.wt"]], 0.5653053, tolerance=1e-5)
  expect_equal(exp(.f$theta[["tka"]]), 1.565491, tolerance=1e-5)
})

test_that("Monolix 2024R1 MAP estimates follow the ini() priors", {
  skip_if_not(file.exists(test_path("monolix2024-priors.zip")))
  .mle <- .mlx2024Fit("mle")

  # a tight prior pins the covariate effect at its mean
  .f <- .mlx2024Fit("beta_strong")
  expect_equal(.f$theta[["cl.wt"]], 0.75, tolerance=1e-3)
  expect_gt(abs(.mle$theta[["cl.wt"]] - 0.75), 0.1)

  # a vague prior leaves the estimate at the MLE
  .f <- .mlx2024Fit("ka_weak")
  expect_equal(exp(.f$theta[["tka"]]), exp(.mle$theta[["tka"]]), tolerance=0.05)
})

test_that("Monolix 2024R1 MAP estimates follow logit and joint priors", {
  skip_if_not(file.exists(test_path("monolix2024-priors.zip")))
  # logitNormal prior: the mean is back-transformed, the sd is on the logit scale
  .f <- .mlx2024Fit("logit_strong")
  expect_equal(expit(.f$theta[["tfr"]]), 0.6, tolerance=1e-3)

  .f <- .mlx2024Fit("two_priors")
  expect_equal(exp(.f$theta[["tka"]]), 3, tolerance=0.01)
  expect_equal(exp(.f$theta[["tcl"]]), 2, tolerance=0.01)
})

test_that("Monolix 2024R1 ignores a MAP prior on a residual error parameter", {
  skip_if_not(file.exists(test_path("monolix2024-priors.zip")))
  .zip <- normalizePath(test_path("monolix2024-priors.zip"))
  withr::with_tempdir({
    utils::unzip(.zip)
    .pop <- function(case) {
      read.csv(file.path("monolix2024", case, paste0(case, "-monolix"),
                         "populationParameters.txt"))[, c("parameter", "value")]
    }
    .mlx <- readLines(file.path("monolix2024", "addsd_strong", "addsd_strong-monolix.mlxtran"))
    # the prior was written ...
    expect_true(any(.mlx == "add__sd={value=0.7, method=MAP}"))
    expect_true(any(.mlx == "add__sd = {distribution=normal, typical=1.5, sd=0.001}"))
    # ... and every estimate is identical to the run without it
    expect_equal(.pop("addsd_strong"), .pop("mle"))
  })
  # so babelmixr2 refuses it rather than let the prior silently vanish
  .u <- ini(.mlx2024Ui("ka_weak"), prior(add.sd) ~ dnorm(1.5, 0.001))
  expect_error(.u$mlxtranParameter, "residual error parameter 'add.sd'")
})
