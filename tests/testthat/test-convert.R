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
    data.frame(
      id = c(1, 2, 2, 3),
      amt = c(10, 10, 0, 0),
      evid = c(1, 1, 0, 0),
      time = c(0, 0, 1, 1),
      DV = 1
    )
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
    data.frame(
      id = c(1, 2, 1, 2),
      time = c(0, 0, 1, 1),
      amt=c(10, 10, NA, NA),
      ii=c(1, NA, NA, NA),
      addl=c(1, NA, NA, NA),
      evid =c(1, 1, 0, 0),
      DV = 1
  )
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

  # Drop right censored subject
  et <-
    data.frame(
      amt=c(10, 10, NA, NA),
      id=c(1L, 2L, 1L, 2L),
      evid=c(1, 1, NA, NA),
      time=c(0, 0, 1, 1),
      cens=c(NA, NA, -1, NA),
      DV=c(NA, NA, 1, 1)
    )
  suppressMessages(
    expect_message(expect_message(
      dClean <- bblDatToPknca(one.compartment, et),
      regexp = "Right censoring and left censoring with a value above zero is not supported with PKNCA estimation, dropping subjects with those censoring types: 1 rows"),
      regexp = "Dropping 1 dosing rows with no observations for the subject with PKNCA estimation"
    )
  )
  expect_named(dClean, c("obs", "dose"))
  expect_equal(nrow(dClean$obs), 1)
  expect_equal(nrow(dClean$dose), 1)
  expect_equal(dClean$obs$id, 2)
  expect_equal(dClean$dose$id, 2)

  # Drop left censored subject with LIMIT > 0
  et <-
    data.frame(
      amt=c(10, 10, NA, NA),
      id=c(1L, 2L, 1L, 2L),
      evid=c(1, 1, NA, NA),
      time=c(0, 0, 1, 1),
      cens=c(NA, NA, 1, NA),
      limit=c(NA, NA, 0.5, NA),
      DV=c(NA, NA, 1, 1)
    )
  suppressMessages(
    expect_message(expect_message(
      dClean <- bblDatToPknca(one.compartment, et),
      regexp = "Right censoring and left censoring with a value above zero is not supported with PKNCA estimation, dropping subjects with those censoring types: 1 rows"),
      regexp = "Dropping 1 dosing rows with no observations for the subject with PKNCA estimation"
    )
  )
  expect_named(dClean, c("obs", "dose"))
  expect_equal(nrow(dClean$obs), 1)
  expect_equal(nrow(dClean$dose), 1)
  expect_equal(dClean$obs$id, 2)
  expect_equal(dClean$dose$id, 2)

  # Keep left censored subject with no LIMIT, convert DV to 0
  et <-
    data.frame(
      amt=c(10, 10, NA, NA),
      id=c(1L, 2L, 1L, 2L),
      evid=c(1, 1, NA, NA),
      time=c(0, 0, 1, 1),
      cens=c(NA, NA, 1, NA),
      DV=c(NA, NA, 1, 1)
    )
  suppressMessages(
    expect_message(
      dClean <- bblDatToPknca(one.compartment, et),
      regexp = "Setting DV to zero for PKNCA estimation with left censoring: 1 rows"
    )
  )
  expect_named(dClean, c("obs", "dose"))
  expect_equal(nrow(dClean$obs), 2)
  expect_equal(nrow(dClean$dose), 2)
  expect_equal(dClean$obs$id, 1:2)
  expect_equal(dClean$dose$id, 1:2)
  expect_equal(dClean$obs$DV, 0:1)

  # Keep left censored subject with zero LIMIT, convert DV to 0
  et <-
    data.frame(
      amt=c(10, 10, NA, NA),
      id=c(1L, 2L, 1L, 2L),
      evid=c(1, 1, NA, NA),
      time=c(0, 0, 1, 1),
      cens=c(NA, NA, 1, NA),
      limit=c(NA, NA, 0, NA),
      DV=c(NA, NA, 1, 1)
    )
  suppressMessages(
    expect_message(
      dClean <- bblDatToPknca(one.compartment, et),
      regexp = "Setting DV to zero for PKNCA estimation with left censoring: 1 rows"
    )
  )
  expect_named(dClean, c("obs", "dose"))
  expect_equal(nrow(dClean$obs), 2)
  expect_equal(nrow(dClean$dose), 2)
  expect_equal(dClean$obs$id, 1:2)
  expect_equal(dClean$dose$id, 1:2)
  expect_equal(dClean$obs$DV, 0:1)

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

test_that("invalid nonmem conversion", {
  
  skip_if_not(file.exists("bad-nonmem-data-convert.qs"))
  
  d <- qs::qread("bad-nonmem-data-convert.qs")

  f <-
    function() {
      ini({
        
        tpm1 <- c(log(0.0001),log(1))
        tpm2 <- c(log(0.0001),log(5))
        tpm3 <- c(log(0.0001),log(0.1))
        tpm4 <- c(log(0.0001),log(5))
        tpm5 <- c(log(0.0001), log(1), log(10))
        tpm6 <- c(log(0.0001), log(0.1), log(1))
        
        eta.pm1 ~ 0.1
        eta.pm2 ~ 0.1
        eta.pm6 ~ 0.1
        eta.pm3 ~ 0.1
        eta.pm4 ~ 0.1
        
        eps.prop <- c(0,1)
        
        
      })
      
      model({
        ipm1 <- exp(tpm1 + eta.pm1)
        pm2 <- exp(tpm2 + eta.pm2)
        pm5 <- exp(tpm5)
        pm6 <- exp(tpm6 + eta.pm6)
        pm3 <- exp(tpm3 + eta.pm3)
        pm4 <- exp(tpm4 + eta.pm4)
        
        tmevent <- tevent - 7
        if(tmevent<0){tmevent = 0}
        
        if(time>tmevent){
          pm5=pm5
          pm6=pm6
        } else{
          pm5=1
          pm6=1
        }
        
        pm1 <- ipm1*pm5
        pm8 <- pm1/pm2
        pm9 <- pm3/pm2
        pm99 <- pm3/pm4
        
        pm7 = pm6
        cp = (cent/pm2)*pm7
        
        d/dt(cent) = -pm8*cp - pm9*cp + pm99*periph
        d/dt(periph) = pm9*cp - pm99*periph
        
        IPRED = log(cp)
        IPRED ~ prop(eps.prop) | abc
        
      })
    }

  d$tevent <- 0.5

  expect_error(bblDatToNonmem(f, d), NA)
  
})
