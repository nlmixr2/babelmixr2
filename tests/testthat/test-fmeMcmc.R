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

  })


})
