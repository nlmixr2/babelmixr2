test_that("NONMEM failures are classified from the output (#46)", {
  expect_null(.nonmemClassifyFailure(NULL))
  expect_null(.nonmemClassifyFailure(c("1NONLINEAR MIXED EFFECTS MODEL PROGRAM",
                                       " License Registered to: Someone",
                                       " Expiration Date:    14 JUN 2099",
                                       "0MINIMIZATION SUCCESSFUL")))
  expect_equal(.nonmemClassifyFailure(c("", "NONMEM license has expired",
                                        "contact ICON"))$cause,
               "license")
  .data <- .nonmemClassifyFailure(c("",
                                    " (DATA ERROR) RECORD         3, DATA ITEM   6, CONTENTS: 1",
                                    " ITEM IS OUT OF RANGE."))
  expect_equal(.data$cause, "data")
  expect_equal(.data$lines,
               c(" (DATA ERROR) RECORD         3, DATA ITEM   6, CONTENTS: 1",
                 " ITEM IS OUT OF RANGE."))
  .ctl <- .nonmemClassifyFailure(c(" AN ERROR WAS FOUND IN THE CONTROL STATEMENTS.",
                                   "",
                                   " AT LINE 12: 208 UNDEFINED VARIABLE: RXQ"))
  expect_equal(.ctl$cause, "controlStream")
  expect_equal(.ctl$lines[2], " AT LINE 12: 208 UNDEFINED VARIABLE: RXQ")
  expect_equal(.nonmemClassifyFailure(c("1NONLINEAR MIXED EFFECTS MODEL PROGRAM",
                                        "0PROGRAM TERMINATED BY OBJ",
                                        " MESSAGE ISSUED FROM ESTIMATION STEP"))$cause,
               "crash")
  # NM-TRAN stopping before NONMEM starts is not an estimation crash
  expect_equal(.nonmemClassifyFailure(c(" WARNING: THE NUMBER OF WARNINGS EXCEEDS THE MAXIMUM.",
                                        " PROGRAM TERMINATED."))$cause,
               "nmtran")
  # the registration line is not a license failure
  expect_null(.nonmemClassifyFailure(c(" License Registered to: Missing Data Solutions",
                                       "1NONLINEAR MIXED EFFECTS MODEL PROGRAM")))
  # a license warning does not hide a crash
  .crash <- .nonmemClassifyFailure(c(" WARNING: LICENSE EXPIRED, GRACE PERIOD",
                                     "1NONLINEAR MIXED EFFECTS MODEL PROGRAM",
                                     "0PROGRAM TERMINATED BY OBJ"))
  expect_equal(.crash$cause, "crash")
  expect_true("0PROGRAM TERMINATED BY OBJ" %in% .crash$lines)
  # the model name is not read as a license message
  expect_null(.nonmemClassifyFailure(c(" PROBLEM NO.:  1  license_missing translated from babelmixr2"),
                                     modelName="license_missing"))
})

withr::with_tempdir({

  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -depot*ka
      d/dt(central) <- depot*ka - cl*central/v
      cp <- central/v
      cp ~ add(add.sd)
    })
  }

  # a runCommand that writes `lines` as NONMEM's output
  .fakeNonmem <- function(lines) {
    force(lines)
    function(ctl, directory, ui) {
      if (!is.null(lines)) {
        writeLines(lines, file.path(directory, ui$nonmemNmlst))
      }
      invisible()
    }
  }

  .fit <- function(runCommand, modelName) {
    nlmixr2(one.cmt, nlmixr2data::theo_sd, "nonmem",
            nonmemControl(runCommand=runCommand, modelName=modelName))
  }

  # the error message from fitting with a runCommand
  .failure <- function(runCommand, modelName) {
    .e <- tryCatch(suppressMessages(.fit(runCommand, modelName)),
                   error=function(e) e)
    expect_s3_class(.e, "error")
    conditionMessage(.e)
  }

  test_that("a missing NONMEM command is reported (#46)", {
    skip_on_os("windows")
    .msg <- .failure("babelmixr2-no-such-nmfe", "fail_cmd")
    expect_match(.msg, "did not create its output file")
    expect_match(.msg, "exit status: 127", fixed=TRUE)
    expect_match(.msg, "nonmemControl(runCommand=)", fixed=TRUE)
  })

  test_that("a runCommand function without output is reported (#46)", {
    .msg <- .failure(.fakeNonmem(NULL), "fail_none")
    expect_match(.msg, "did not create its output file")
    expect_no_match(.msg, "exit status")
  })

  test_that("NONMEM license failures are reported (#46)", {
    .msg <- .failure(.fakeNonmem(c("License file nonmem.lic has expired",
                         "Please contact ICON")), "fail_lic")
    expect_match(.msg, "license problem")
    expect_match(.msg, "nonmem.lic has expired", fixed=TRUE)
  })

  test_that("NM-TRAN control stream errors are reported (#46)", {
    .msg <- .failure(.fakeNonmem(c(" AN ERROR WAS FOUND IN THE CONTROL STATEMENTS.",
                         " AT LINE 12: 208 UNDEFINED VARIABLE: RXQ")), "fail_ctl")
    expect_match(.msg, "UNDEFINED VARIABLE: RXQ", fixed=TRUE)
    expect_match(.msg, "fail_ctl.nmctl", fixed=TRUE)
    expect_match(.msg, "babelmixr2/issues", fixed=TRUE)
  })

  test_that("NM-TRAN data errors are reported (#46)", {
    .msg <- .failure(.fakeNonmem(c(" (DATA ERROR) RECORD         3, DATA ITEM   6, CONTENTS: 1",
                         " ITEM IS OUT OF RANGE.")), "fail_data")
    expect_match(.msg, "ITEM IS OUT OF RANGE", fixed=TRUE)
    expect_match(.msg, "fail_data.csv", fixed=TRUE)
  })

  test_that("NONMEM that never starts estimation is reported (#46)", {
    .msg <- .failure(.fakeNonmem(c("$PROBLEM translated from babelmixr2",
                         "gfortran: error: cannot compile")), "fail_start")
    expect_match(.msg, "did not start estimation")
    expect_match(.msg, "cannot compile", fixed=TRUE)
  })

  test_that("NONMEM crashes during estimation are reported (#46)", {
    .msg <- .failure(.fakeNonmem(c("1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1",
                         " ITERATION NO.:   15    OBJECTIVE VALUE:   120.5")), "fail_crash")
    expect_match(.msg, "started but did not finish")
    expect_match(.msg, "ITERATION NO.:   15", fixed=TRUE)

    .msg <- .failure(.fakeNonmem(c("1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1",
                         "0PROGRAM TERMINATED BY OBJ",
                         " ERROR IN NCONTR WHILE COMPUTING OBJECTIVE")), "fail_obj")
    expect_match(.msg, "PROGRAM TERMINATED BY OBJ", fixed=TRUE)
  })

  test_that("a crash after NONMEM read the data is reported (#46)", {
    # with the data summary, nonmem2rx reads the crash as NONMEM's
    # termination message
    .msg <- .failure(.fakeNonmem(c("1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1",
                                   " TOT. NO. OF OBS RECS:      132",
                                   " TOT. NO. OF INDIVIDUALS:       12",
                                   " #TERM:",
                                   "0PROGRAM TERMINATED BY OBJ",
                                   " ERROR IN NCONTR WHILE COMPUTING OBJECTIVE")), "fail_obj2")
    expect_match(.msg, "NONMEM stopped during the run")
    expect_match(.msg, "PROGRAM TERMINATED BY OBJ", fixed=TRUE)
    expect_no_match(.msg, "minimization not successful")
  })

  test_that("a finished NONMEM run that did not converge is reported (#46)", {
    .msg <- .failure(.fakeNonmem(c("1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1",
                                   " TOT. NO. OF OBS RECS:      132",
                                   " TOT. NO. OF INDIVIDUALS:       12",
                                   " #TERM:",
                                   "0MINIMIZATION TERMINATED",
                                   " DUE TO MAX. NO. OF FUNCTION EVALUATIONS EXCEEDED",
                                   " #TERE:")), "fail_maxeval")
    expect_match(.msg, "minimization not successful")
  })

  test_that("unreadable NONMEM output is reported (#46)", {
    .msg <- .failure(.fakeNonmem(c("1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1",
                         " #TERM:",
                         " garbled")), "fail_read")
    expect_match(.msg, "could not read NONMEM's output")
  })

  test_that("the echoed control stream is not read as a NONMEM message (#46)", {
    .dir <- file.path(tempfile(), "echo-nonmem")
    dir.create(.dir, recursive=TRUE)
    writeLines(c("$PK", "  LICENSE_MISSING=1", "; PROGRAM TERMINATED"),
               file.path(.dir, "echo.nmctl"))
    writeLines(c("$PK", "  LICENSE_MISSING=1", "; PROGRAM TERMINATED",
                 "gfortran: error: cannot compile"),
               file.path(.dir, "echo.lst"))
    .lines <- .nonmemFailureLines(.dir, "echo.lst", "echo.nmctl")
    expect_equal(.lines, "gfortran: error: cannot compile")
    expect_null(.nonmemClassifyFailure(.lines))
    unlink(dirname(.dir), recursive=TRUE)
  })

  test_that("errors reading a finished NONMEM run are reported (#46)", {
    .ui <- suppressMessages(.fit(NA, "fail_final"))
    writeLines(c("1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1",
                 " #TERM:",
                 "0MINIMIZATION SUCCESSFUL"),
               file.path(.ui$nonmemExportPath, .ui$nonmemNmlst))
    local_mocked_bindings(.nonmemFinalizeEnv=function(env, oldUi) {
      stop("mock post-processing error", call.=FALSE)
    })
    expect_error(.nonmemFinalizeOrExplain(new.env(), .ui),
                 "could not read NONMEM's output(.|\n)*mock post-processing error")
    local_mocked_bindings(.nonmemFinalizeEnv=function(env, oldUi) "fit")
    expect_equal(.nonmemFinalizeOrExplain(new.env(), .ui), "fit")
  })

  test_that("an unset run command says how to set it (#46)", {
    expect_error(suppressMessages(.fit("", "fail_unset")),
                 "nonmemControl(runCommand=)", fixed=TRUE)
  })
})
