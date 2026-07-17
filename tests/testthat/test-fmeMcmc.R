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
