test_that("bioavailability makes sense (issue #89)", {
  skip_on_cran()

  one.cmt <- function() {
    ini({
      lfdepotHigh <- 0.1
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
      prop.sd <- 0.01
    })
    model({
      fdepot <- exp(lfdepotHigh*HIGH)
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      kel <- cl/v
      d/dt(depot) <- -ka*depot
      d/dt(central) <- ka*depot - kel*central
      cp <- central/v
      f(depot) <- fdepot
      cp ~ add(add.sd) + prop(prop.sd)
    })
  }


  f <- one.cmt()

  f2 <- strsplit(f$nonmemModel, "\n")[[1]]

  f2 <- gsub(" *", "",
             gsub(" *;.*", "",
                  f2[which(regexpr("^ *F1=", f2) != -1)]))

  expect_equal(f2, "F1=DEXP(HIGH*THETA(1))")
  
})
