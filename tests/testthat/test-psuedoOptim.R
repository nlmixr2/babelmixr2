test_that("fmeMcmc works", {
  skip_if_not_installed("FME")
  skip_on_cran()

  rxode2::rxWithSeed(42,{

    dsn <- data.frame(i=1:1000)
    dsn$time <- exp(rnorm(1000))
    dsn$DV=rbinom(1000,1,exp(-1+dsn$time)/(1+exp(-1+dsn$time)))

    mod <- function() {
      ini({
        # This estimation method requires all parameters
        # to be bounded:
        E0 <- c(-100, 0.5, 100)
        Em <- c(0, 0.5, 10)
        E50 <- c(0, 2, 20)
        g <- fix(c(0.1, 2, 10))
      })
      model({
        v <- E0+Em*time^g/(E50^g+time^g)
        ll(bin) ~ DV * v - log(1 + exp(v))
      })
    }

    fit2 <- suppressMessages(nlmixr(mod, dsn, est="pseudoOptim"))

    expect_s3_class(fit2, "nlmixr2.pseudoOptim")

    mod <- function() {
      ini({
        # This estimation method requires all parameters
        # to be bounded:
        E0 <- c(-100, 0.5)
        Em <- c(0, 0.5, 10)
        E50 <- c(0, 2, 20)
        g <- fix(c(0.1, 2, 10))
      })
      model({
        v <- E0+Em*time^g/(E50^g+time^g)
        ll(bin) ~ DV * v - log(1 + exp(v))
      })
    }

    expect_error(suppressMessages(nlmixr(mod, dsn, est="pseudoOptim")))

    mod <- function() {
      ini({
        # This estimation method requires all parameters
        # to be bounded:
        E0 <- c(-Inf, 0.5, 100)
        Em <- c(0, 0.5, 10)
        E50 <- c(0, 2, 20)
        g <- fix(c(0.1, 2, 10))
      })
      model({
        v <- E0+Em*time^g/(E50^g+time^g)
        ll(bin) ~ DV * v - log(1 + exp(v))
      })
    }

    expect_error(suppressMessages(nlmixr(mod, dsn, est="pseudoOptim")))

  })

})
