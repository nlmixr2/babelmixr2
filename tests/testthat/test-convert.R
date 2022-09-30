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

test_that("pknca conversion keeps extra columns", {
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

  et <- rxode2::et(amt=10) %>% rxode2::et(1)
  et$DV <- 100
  bblDatToPknca(one.compartment, et)

  et <- rxode2::et(amt=0) %>% rxode2::et(1)
  et$DV <- 100
  suppressMessages(expect_error(
    bblDatToPknca(one.compartment, et),
    regexp="no dosing rows (EVID = 1 or 4) detected",
    fixed=TRUE
  ))
})

test_that("getStandardColNames", {
  expect_equal(
    getStandardColNames(data.frame(ID=1, DV=2, Time=3, CmT=4)),
    c(id = "ID", time = "Time", amt = NA, rate = NA, dur = NA, evid = NA,
      cmt = "CmT", ss = NA, ii = NA, addl = NA, dv = "DV", mdv = NA, dvid = NA,
      cens = NA, limit = NA)
  )
  expect_error(
    getStandardColNames(data.frame(ID=1, DV=2, Time=3, CmT=4, cmt=5)),
    regexp = "Multiple data columns match cmt when converted to lower case"
  )
})
