test_that("nonmem $pred model", {

  dat <- nlmixr2data::Wang2007
  dat$DV <- dat$Y

  ipredWang2007 <- function() {
    ini({
      tke <- 0.5
      eta.ke ~ 0.04
      prop.sd <- sqrt(0.1)
    })
    model({
      ke <- tke * exp(eta.ke)
      ipre <- 10 * exp(-ke * t)
      f2 <- ipre / (ipre + 5)
      f3 <- f2 * 3
      lipre <- log(ipre)
      ipre ~ prop(prop.sd)
    })
  }

  skip_if_not(file.exists("ipredWang2007.zip"))
  .path <- normalizePath("ipredWang2007.zip")
  withr::with_tempdir({
    unzip(.path)
    expect_error(nlmixr(ipredWang2007, dat, "nonmem"), NA)
  })
})
