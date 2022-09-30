test_that("nonmem amt=0 evid=1 conversion test", {
  
  one.compartment <- function() {
    ini({
      tka <- 0.45 # Log Ka
      tcl <- 1 # Log Cl
      tv <- 3.45    # Log V
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl / v * center
      cp = center / v
      cp ~ add(add.sd)
    })
  }

  et <- rxode2::et(amt=0) %>% rxode2::et(1)
  et$DV <- 100

  conv <- bblDatToNonmem(one.compartment, et)
  expect_equal(conv$EVID, c(2L, 0L))
  
})
