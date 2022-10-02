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

  # Normal
  et <- rxode2::et(amt=10) %>% rxode2::et(1)
  et$DV <- 100
  suppressMessages(dClean <- bblDatToPknca(one.compartment, et))
  expect_named(dClean, c("obs", "dose"))
  expect_equal(nrow(dClean$obs), 1)
  expect_equal(nrow(dClean$dose), 1)

  # Only dosing and only observation rows for some subjects
  et <-
    dplyr::bind_rows(
      as.data.frame(rxode2::et(amt=10, id=1:2)),
      as.data.frame(rxode2::et(time=1, id=2:3))
    )
  et$DV <- 1
  suppressWarnings(suppressMessages(
    expect_message(
      dClean <- bblDatToPknca(one.compartment, et),
      regexp = "Dropping 1 observation rows with no doses for the subject with PKNCA estimation"
    )
  ))
  expect_named(dClean, c("obs", "dose"))
  expect_equal(nrow(dClean$obs), 1)
  expect_equal(nrow(dClean$dose), 1)
  expect_equal(dClean$obs$id, 2)
  expect_equal(dClean$dose$id, 2)

  # Drop a whole subject if they use ADDL
  et <-
    dplyr::bind_rows(
      as.data.frame(rxode2::et(amt=10, id=1, addl=1, ii=1)),
      as.data.frame(rxode2::et(amt=10, id=2)),
      as.data.frame(rxode2::et(time=1, id=1:2))
    )
  et$DV <- 1
  suppressWarnings(suppressMessages(
    expect_message(expect_message(
      dClean <- bblDatToPknca(one.compartment, et),
      regexp = "ADDL dosing not supported with PKNCA estimation, dropping subjects using ADDL: 1 rows"),
      regexp = "Dropping 1 observation rows with no doses for the subject with PKNCA estimation"
    )
  ))
  expect_named(dClean, c("obs", "dose"))
  expect_equal(nrow(dClean$obs), 1)
  expect_equal(nrow(dClean$dose), 1)
  expect_equal(dClean$obs$id, 2)
  expect_equal(dClean$dose$id, 2)

  # No dosing
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
