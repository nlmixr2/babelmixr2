test_that("fmeMcmc works", {
  skip_if_not_installed("FME")
  skip_if_not_installed("coda")
  skip_on_cran()

  rxode2::rxWithSeed(42,{
    dsn <- data.frame(i=1:1000)
    dsn$time <- exp(rnorm(1000))
    dsn$DV=rbinom(1000,1,exp(-1+dsn$time)/(1+exp(-1+dsn$time)))

    mod <- function() {
      ini({
        E0 <- 0.5
        Em <- 0.5
        E50 <- 2
        g <- fix(2)
      })
      model({
        v <- E0+Em*time^g/(E50^g+time^g)
        ll(bin) ~ DV * v - log(1 + exp(v))
      })
    }


    fit1 <- suppressMessages(nlmixr(mod, dsn, est="fmeMcmc",
                                    control=fmeMcmcControl(print=0)))
    expect_s3_class(fit1, "nlmixr2.fmeMcmc")

    fit2 <- suppressMessages(nlmixr(mod, dsn, est="fmeMcmc",
                                    control=fmeMcmcControl(scaleType="nlmixr2",
                                                           print=0)))

    expect_s3_class(fit2, "nlmixr2.fmeMcmc")

    expect_s3_class(coda::as.mcmc(fit2), "mcmc")

    # issue #178: when the chain is explored in a scaled space, `$pars` must be
    # back-transformed to the natural parameter space.
    .chain <- fit2$fmeMcmc$pars

    # the chain columns carry the real parameter names, not FME's generic
    # `p1`/`p2`/... labels
    expect_setequal(colnames(.chain), c("E0", "Em", "E50"))

    # the natural-scale fixed-effect estimates lie within the explored chain.
    # `Estimate` is `bestpar`, which is by definition a member of the chain, so
    # this holds exactly for a back-transformed chain -- but fails for the old
    # scaled chain, whose samples cluster near the ~1 scaled starting point.
    .est <- setNames(fit2$parFixedDf$Estimate, rownames(fit2$parFixedDf))
    .est <- .est[colnames(.chain)]
    expect_true(all(.est >= apply(.chain, 2L, min) &
                      .est <= apply(.chain, 2L, max)))

    # the covariance (used for the standard errors) is the sample covariance of
    # the back-transformed chain
    expect_equal(fit2$fmeMcmc$cov, stats::cov(.chain))

  })


})

test_that("fmeMcmc chain is determined by control's seed", {
  skip_if_not_installed("FME")
  skip_on_cran()

  dsn <- rxode2::rxWithSeed(42, {
    .d <- data.frame(i=1:1000)
    .d$time <- exp(rnorm(1000))
    .d$DV <- rbinom(1000, 1, exp(-1+.d$time)/(1+exp(-1+.d$time)))
    .d
  })

  mod <- function() {
    ini({
      E0 <- 0.5
      Em <- 0.5
      E50 <- 2
      g <- fix(2)
    })
    model({
      v <- E0+Em*time^g/(E50^g+time^g)
      ll(bin) ~ DV * v - log(1 + exp(v))
    })
  }

  .fit <- function(seed) {
    suppressMessages(nlmixr(mod, dsn, est="fmeMcmc",
                            control=fmeMcmcControl(print=0, niter=50,
                                                   seed=seed)))$fmeMcmc$pars
  }

  # the run is wrapped in rxWithSeed(), so the same seed reproduces the chain
  # exactly even though the two fits run back-to-back from different RNG states
  expect_equal(.fit(42), .fit(42))

  # ... and a different seed explores a different chain
  expect_false(isTRUE(all.equal(.fit(42), .fit(7))))
})

test_that("fmeMcmc samples with ini({}) priors (#208)", {
  skip_if_not_installed("FME")
  skip_on_cran()
  skip_if_not(exists("rxPriorLogDensity", envir=asNamespace("rxode2")),
              "rxode2 does not support ini({}) priors")

  dsn <- rxode2::rxWithSeed(42, {
    .d <- data.frame(i=1:1000)
    .d$time <- exp(rnorm(1000))
    .d$DV <- rbinom(1000, 1, exp(-1+.d$time)/(1+exp(-1+.d$time)))
    .d
  })

  # a tight prior far from where the data put E0 (~ -0.5)
  mod <- function() {
    ini({
      E0 <- 0.5
      Em <- 0.5
      E50 <- 2
      g <- fix(2)
      prior(E0) ~ dnorm(5, 0.01)
      prior(Em) ~ dnorm(0.5, 2)
    })
    model({
      v <- E0+Em*time^g/(E50^g+time^g)
      ll(bin) ~ DV * v - log(1 + exp(v))
    })
  }

  # the prior is evaluated on the natural scale, so it pulls the chain to the
  # same place whether or not FME samples a rescaled space
  for (.scaleType in c("none", "nlmixr2")) {
    fit <- suppressMessages(nlmixr(mod, dsn, est="fmeMcmc",
                                   control=fmeMcmcControl(print=0,
                                                          scaleType=.scaleType)))
    expect_s3_class(fit, "nlmixr2.fmeMcmc")
    .pars <- fit$fmeMcmc$pars
    .e0 <- .pars[seq(nrow(.pars) %/% 2L, nrow(.pars)), "E0"]
    expect_equal(mean(.e0), 5, tolerance=0.02)
    # FME stores -2*log(prior) for each sample; it is the joint declared prior
    .n <- nrow(.pars)
    expect_equal(fit$fmeMcmc$prior[.n],
                 -2 * (stats::dnorm(.pars[[.n, "E0"]], 5, 0.01, log=TRUE) +
                         stats::dnorm(.pars[[.n, "Em"]], 0.5, 2, log=TRUE)),
                 tolerance=1e-6)
  }

  # supplying both an ini({}) prior and fmeMcmcControl(prior=) is an error
  expect_error(nlmixr(mod, dsn, est="fmeMcmc",
                      control=fmeMcmcControl(print=0, prior=function(p) 0)),
               "fmeMcmcControl\\(prior=\\)")

  # an omega prior never reaches the prior kernel: the model is refused as
  # not population-only first
  modEta <- function() {
    ini({
      E0 <- 0.5
      Em <- 0.5
      E50 <- 2
      g <- fix(2)
      eta.E0 ~ 0.1
      prior(E0) ~ dnorm(0, 10)
      prior(eta.E0) ~ dnorm(0.1, 0.01)
    })
    model({
      v <- E0+eta.E0+Em*time^g/(E50^g+time^g)
      ll(bin) ~ DV * v - log(1 + exp(v))
    })
  }
  expect_error(nlmixr(modEta, dsn, est="fmeMcmc",
                      control=fmeMcmcControl(print=0)),
               "population estimates")
})

test_that("fmeMcmc fits linCmt() models", {
  skip_if_not_installed("FME")
  skip_if_not_installed("coda")
  skip_on_cran()
  lin <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  fit <- suppressMessages(nlmixr(lin, nlmixr2data::theo_sd, est="fmeMcmc",
                                 control=fmeMcmcControl(print=0, niter=50)))
  expect_s3_class(fit, "nlmixr2.fmeMcmc")
})

test_that("fmeMcmc honors jump and moves zero initial estimates (#276)", {
  skip_if_not_installed("FME")
  skip_on_cran()

  d <- rxode2::rxWithSeed(1, {
    .d <- data.frame(id = 1, time = 1:20)
    .d$dv <- 1 + 0.5 * .d$time + rnorm(20, sd = 0.3)
    .d
  })

  mod <- function() {
    ini({
      a <- 0
      b <- 0.4
      add.sd <- fix(0.3)
    })
    model({
      y <- a + b * time
      y ~ add(add.sd)
    })
  }

  # Default jump = NULL must still move a parameter initialized at 0
  fitDefault <- suppressMessages(nlmixr(
    mod, d, est = "fmeMcmc",
    control = fmeMcmcControl(niter = 400, print = 0, seed = 99)
  ))
  expect_gt(length(unique(fitDefault$fmeMcmc$pars[, "a"])), 10L)

  # User-supplied vector jump is honored and changes the sampled chain
  fitVec <- suppressMessages(nlmixr(
    mod, d, est = "fmeMcmc",
    control = fmeMcmcControl(niter = 400, jump = c(0.2, 0.05), print = 0, seed = 99)
  ))
  expect_gt(length(unique(fitVec$fmeMcmc$pars[, "a"])), 10L)
  expect_false(isTRUE(all.equal(fitDefault$fmeMcmc$pars, fitVec$fmeMcmc$pars)))

  # Natural-scale vector jump maps through parameter scaling to produce the
  # same natural-scale chain under scaleType = "nlmixr2"
  fitVecScaled <- suppressMessages(nlmixr(
    mod, d, est = "fmeMcmc",
    control = fmeMcmcControl(
      niter = 400, jump = c(0.2, 0.05), scaleType = "nlmixr2", print = 0, seed = 99
    )
  ))
  expect_equal(fitVecScaled$fmeMcmc$pars, fitVec$fmeMcmc$pars, tolerance = 1e-10)

  # Covariance matrix jump is honored and invariant across scaleType
  .covJump <- matrix(c(0.04, -0.002, -0.002, 0.0025), nrow = 2)
  fitMat <- suppressMessages(nlmixr(
    mod, d, est = "fmeMcmc",
    control = fmeMcmcControl(niter = 400, jump = .covJump, print = 0, seed = 99)
  ))
  fitMatScaled <- suppressMessages(nlmixr(
    mod, d, est = "fmeMcmc",
    control = fmeMcmcControl(
      niter = 400, jump = .covJump, scaleType = "nlmixr2", print = 0, seed = 99
    )
  ))
  expect_gt(length(unique(fitMat$fmeMcmc$pars[, "a"])), 10L)
  expect_equal(fitMatScaled$fmeMcmc$pars, fitMat$fmeMcmc$pars, tolerance = 1e-10)

  # Scalar jump applies to all parameters
  fitScalar <- suppressMessages(nlmixr(
    mod, d, est = "fmeMcmc",
    control = fmeMcmcControl(niter = 100, jump = 0.05, print = 0, seed = 99)
  ))
  expect_gt(length(unique(fitScalar$fmeMcmc$pars[, "a"])), 5L)

  # Invalid jump dimensions error clearly
  expect_error(
    suppressMessages(nlmixr(
      mod, d, est = "fmeMcmc",
      control = fmeMcmcControl(niter = 50, jump = c(0.1, 0.2, 0.3), print = 0)
    )),
    "jump"
  )
})
