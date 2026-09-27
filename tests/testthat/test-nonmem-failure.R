test_that("NONMEM failures are classified from the output (#46)", {
  expect_null(.nonmemClassifyFailure(NULL))
  expect_null(.nonmemClassifyFailure(c(
    "1NONLINEAR MIXED EFFECTS MODEL PROGRAM",
    " License Registered to: Someone",
    " Expiration Date:    14 JUN 2099",
    "0MINIMIZATION SUCCESSFUL"
  )))
  expect_equal(
    .nonmemClassifyFailure(c(
      "",
      "NONMEM license has expired",
      "contact ICON"
    ))$cause,
    "license"
  )
  .data <- .nonmemClassifyFailure(c(
    "",
    " (DATA ERROR) RECORD         3, DATA ITEM   6, CONTENTS: 1",
    " ITEM IS OUT OF RANGE."
  ))
  expect_equal(.data$cause, "data")
  expect_equal(
    .data$lines,
    c(
      " (DATA ERROR) RECORD         3, DATA ITEM   6, CONTENTS: 1",
      " ITEM IS OUT OF RANGE."
    )
  )
  .ctl <- .nonmemClassifyFailure(c(
    " AN ERROR WAS FOUND IN THE CONTROL STATEMENTS.",
    "",
    " AT LINE 12: 208 UNDEFINED VARIABLE: RXQ"
  ))
  expect_equal(.ctl$cause, "controlStream")
  expect_equal(.ctl$lines[2], " AT LINE 12: 208 UNDEFINED VARIABLE: RXQ")
  expect_equal(
    .nonmemClassifyFailure(c(
      "1NONLINEAR MIXED EFFECTS MODEL PROGRAM",
      "0PROGRAM TERMINATED BY OBJ",
      " MESSAGE ISSUED FROM ESTIMATION STEP"
    ))$cause,
    "crash"
  )
  # NM-TRAN stopping before NONMEM starts is not an estimation crash
  expect_equal(
    .nonmemClassifyFailure(c(
      " WARNING: THE NUMBER OF WARNINGS EXCEEDS THE MAXIMUM.",
      " PROGRAM TERMINATED."
    ))$cause,
    "nmtran"
  )
  # the registration line is not a license failure
  expect_null(.nonmemClassifyFailure(c(
    " License Registered to: Missing Data Solutions",
    "1NONLINEAR MIXED EFFECTS MODEL PROGRAM"
  )))
  # a license warning does not hide a crash
  .crash <- .nonmemClassifyFailure(c(
    " WARNING: LICENSE EXPIRED, GRACE PERIOD",
    "1NONLINEAR MIXED EFFECTS MODEL PROGRAM",
    "0PROGRAM TERMINATED BY OBJ"
  ))
  expect_equal(.crash$cause, "crash")
  expect_true("0PROGRAM TERMINATED BY OBJ" %in% .crash$lines)
  # NONMEM stopping after estimation finished is not an estimation crash
  expect_null(.nonmemClassifyFailure(c(
    "1NONLINEAR MIXED EFFECTS MODEL PROGRAM",
    " #TERM:",
    "0MINIMIZATION SUCCESSFUL",
    " #TERE:",
    "0PROGRAM TERMINATED BY FNLETA",
    " MESSAGE ISSUED FROM TABLE STEP"
  )))
  # the model name is not read as a license message
  expect_null(.nonmemClassifyFailure(
    c(
      " PROBLEM NO.:  1  license_missing translated from babelmixr2",
      " DATA FILE: license_missing.csv"
    ),
    modelName = "license_missing"
  ))
  # nor does a model named "license" (real NONMEM license messages)
  for (.l in c(
    " ERROR reading license file /opt/nm730/license/nonmem.lic",
    "  **** NONMEM LICENSE HAS EXPIRED ****",
    "License file has expired"
  )) {
    expect_equal(
      .nonmemClassifyFailure(.l, modelName = "license")$cause,
      "license",
      info = .l
    )
  }
  # but a short model name does not hide a license message
  expect_equal(
    .nonmemClassifyFailure("License file has expired", modelName = "a")$cause,
    "license"
  )
  expect_equal(
    .nonmemDropModelName(c("a.csv a-nonmem/a.lst", "data", "pk.1+a"), "a"),
    c(".csv -nonmem/.lst", "data", "pk.1+")
  )
  expect_equal(.nonmemDropModelName("x (m.1) m.1.csv", "m.1"), "x () .csv")
  # directories are not license messages, but license files are
  expect_null(.nonmemClassifyFailure(
    "gfortran: error: /home/u/missing_license/FSUBS.f90: failed"
  ))
  expect_equal(
    .nonmemClassifyFailure("Cannot find /opt/nm/license/nonmem.lic")$cause,
    "license"
  )
})

test_that("NONMEM exiting with an error is warned about (#46)", {
  expect_no_warning(.nonmemWarnStatus(NULL))
  expect_no_warning(.nonmemWarnStatus(0L))
  expect_warning(.nonmemWarnStatus(137L), "exited with status 137")
})

test_that("real NONMEM output is classified (#46)", {
  # PsN's collection of NONMEM output, shipped with nonmem2rx
  .zip <- system.file("PsN.zip", package = "nonmem2rx")
  skip_if(.zip == "")
  withr::with_tempdir({
    unzip(.zip)
    .cause <- function(f) {
      .l <- .nonmemFailureLines(dirname(f), basename(f), "none.nmctl")
      .c <- .nonmemClassifyFailure(.l, "run1")
      if (is.null(.c)) "none" else .c$cause
    }
    .d <- file.path("PsN", "test_files", "output")
    .expected <- c(
      "special_mod/license_missing.lst" = "license",
      "special_mod/license_dummy.lst" = "license",
      "special_mod/license_expired.lst" = "license",
      "special_mod/data_missing.lst" = "data",
      "onePROB/oneEST/noSIM/hessian_error.lst" = "crash",
      "onePROB/oneEST/noSIM/nm710_fail_negV.lst" = "crash",
      # the first of two estimations crashed
      "onePROB/multEST/firstEstTerm.lst" = "crash",
      # NONMEM stopped in the covariance or table step, after estimation
      "onePROB/oneEST/noSIM/large_s_matrix_cov_fail.lst" = "none",
      "special_mod/interrupted_at_eigen.lst" = "none",
      "nm73/UseCase7.lst" = "none",
      "nm73/mox_fail_nonp.lst" = "none"
    )
    for (.f in names(.expected)) {
      expect_equal(.cause(file.path(.d, .f)), .expected[[.f]], info = .f)
    }
    # the control stream line NM-TRAN shows with an error is kept
    .e <- file.path("PsN", "test_files", "modelfit", "diagnose_lst_errors")
    .ctl <- .nonmemClassifyFailure(
      .nonmemFailureLines(.e, "psn.lst", "psn.mod"),
      "psn"
    )
    expect_equal(.ctl$cause, "controlStream")
    expect_true(any(
      grepl("HEJSAN", .ctl$lines, fixed = TRUE) &
        grepl("$ESTIMATION", .ctl$lines, fixed = TRUE)
    ))
    # no finished run in the collection is read as a failure
    .fs <- list.files(
      .d,
      pattern = "[.]lst$",
      recursive = TRUE,
      full.names = TRUE
    )
    .fs <- .fs[!(.fs %in% file.path(.d, names(.expected)))]
    .done <- vapply(
      .fs,
      function(f) {
        any(grepl("#TERE:", .nonmemFailureReadLines(f), fixed = TRUE))
      },
      logical(1)
    )
    for (.f in .fs[.done]) {
      expect_equal(.cause(.f), "none", info = .f)
    }
  })
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
      d / dt(depot) <- -depot * ka
      d / dt(central) <- depot * ka - cl * central / v
      cp <- central / v
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

  .fit <- function(runCommand, modelName, ...) {
    nlmixr2(
      one.cmt,
      nlmixr2data::theo_sd,
      "nonmem",
      nonmemControl(runCommand = runCommand, modelName = modelName, ...)
    )
  }

  # the error message from fitting with a runCommand
  .failure <- function(runCommand, modelName, ctl = list()) {
    .e <- tryCatch(
      suppressMessages(do.call(.fit, c(list(runCommand, modelName), ctl))),
      error = function(e) e
    )
    expect_s3_class(.e, "error")
    conditionMessage(.e)
  }

  test_that("a missing NONMEM command is reported (#46)", {
    skip_on_os("windows")
    .msg <- .failure("babelmixr2-no-such-nmfe", "fail_cmd")
    expect_match(.msg, "did not create its output file")
    expect_match(.msg, "exit status: 127", fixed = TRUE)
    expect_match(.msg, "nonmemControl(runCommand=)", fixed = TRUE)
  })

  test_that("a runCommand function without output is reported (#46)", {
    .msg <- .failure(.fakeNonmem(NULL), "fail_none")
    expect_match(.msg, "did not create its output file")
    expect_no_match(.msg, "exit status")
  })

  test_that("NONMEM license failures are reported (#46)", {
    .msg <- .failure(
      .fakeNonmem(c(
        "License file nonmem.lic has expired",
        "Please contact ICON"
      )),
      "fail_lic"
    )
    expect_match(.msg, "license problem")
    expect_match(.msg, "nonmem.lic has expired", fixed = TRUE)
  })

  test_that("NM-TRAN control stream errors are reported (#46)", {
    .msg <- .failure(
      .fakeNonmem(c(
        " AN ERROR WAS FOUND IN THE CONTROL STATEMENTS.",
        " AT LINE 12: 208 UNDEFINED VARIABLE: RXQ"
      )),
      "fail_ctl"
    )
    expect_match(.msg, "UNDEFINED VARIABLE: RXQ", fixed = TRUE)
    expect_match(.msg, "fail_ctl.nmctl", fixed = TRUE)
    expect_match(.msg, "babelmixr2/issues", fixed = TRUE)
  })

  test_that("NM-TRAN data errors are reported (#46)", {
    .msg <- .failure(
      .fakeNonmem(c(
        " (DATA ERROR) RECORD         3, DATA ITEM   6, CONTENTS: 1",
        " ITEM IS OUT OF RANGE."
      )),
      "fail_data"
    )
    expect_match(.msg, "ITEM IS OUT OF RANGE", fixed = TRUE)
    expect_match(.msg, "fail_data.csv", fixed = TRUE)
  })

  test_that("NONMEM that never starts estimation is reported (#46)", {
    .msg <- .failure(
      .fakeNonmem(c(
        "$PROBLEM translated from babelmixr2",
        "gfortran: error: cannot compile"
      )),
      "fail_start"
    )
    expect_match(.msg, "did not start estimation")
    expect_match(.msg, "cannot compile", fixed = TRUE)
  })

  test_that("NONMEM crashes during estimation are reported (#46)", {
    .msg <- .failure(
      .fakeNonmem(c(
        "1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1",
        " ITERATION NO.:   15    OBJECTIVE VALUE:   120.5"
      )),
      "fail_crash"
    )
    expect_match(.msg, "started but did not finish")
    expect_match(.msg, "ITERATION NO.:   15", fixed = TRUE)

    .msg <- .failure(
      .fakeNonmem(c(
        "1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1",
        "0PROGRAM TERMINATED BY OBJ",
        " ERROR IN NCONTR WHILE COMPUTING OBJECTIVE"
      )),
      "fail_obj"
    )
    expect_match(.msg, "PROGRAM TERMINATED BY OBJ", fixed = TRUE)
  })

  test_that("a crash after NONMEM read the data is reported (#46)", {
    # as NONMEM writes it: the crash, then an empty termination block
    .msg <- .failure(
      .fakeNonmem(c(
        "1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1",
        " TOT. NO. OF OBS RECS:      132",
        " TOT. NO. OF INDIVIDUALS:       12",
        "0PROGRAM TERMINATED BY OBJ",
        " ERROR IN NCONTR WHILE COMPUTING OBJECTIVE",
        " MESSAGE ISSUED FROM ESTIMATION STEP",
        " #TERM:",
        " #TERE:"
      )),
      "fail_obj2"
    )
    expect_match(.msg, "NONMEM stopped during the run")
    expect_match(.msg, "PROGRAM TERMINATED BY OBJ", fixed = TRUE)
    expect_no_match(.msg, "minimization not successful")
  })

  test_that("a finished NONMEM run that did not converge is reported (#46)", {
    .msg <- .failure(
      .fakeNonmem(c(
        "1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1",
        " TOT. NO. OF OBS RECS:      132",
        " TOT. NO. OF INDIVIDUALS:       12",
        " #TERM:",
        "0MINIMIZATION TERMINATED",
        " DUE TO MAX. NO. OF FUNCTION EVALUATIONS EXCEEDED",
        " #TERE:"
      )),
      "fail_maxeval"
    )
    expect_match(.msg, "minimization not successful")
  })

  test_that("a crash points to NONMEM's solving errors (#46)", {
    .withPrderr <- function(ctl, directory, ui) {
      writeLines(
        c(
          "1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1",
          "0PROGRAM TERMINATED BY PRED"
        ),
        file.path(directory, ui$nonmemNmlst)
      )
      writeLines("ERROR IN LSODA", file.path(directory, "PRDERR"))
    }
    .msg <- .failure(.withPrderr, "fail_prderr")
    expect_match(
      .msg,
      "solving errors: 'fail_prderr-nonmem/PRDERR'",
      fixed = TRUE
    )
  })

  test_that("NM-TRAN stopping before NONMEM runs is reported (#46)", {
    .msg <- .failure(
      .fakeNonmem(c(
        " WARNING: THE NUMBER OF WARNINGS EXCEEDS THE MAXIMUM.",
        " PROGRAM TERMINATED."
      )),
      "fail_nmtran"
    )
    expect_match(.msg, "NM-TRAN stopped before NONMEM could run")
    expect_match(.msg, "fail_nmtran.nmctl", fixed = TRUE)
  })

  test_that("a NONMEM evaluation (MAXEVALS=0) is not a crash (#46)", {
    # evaluations have no termination message, only the omitted
    # estimation step
    .eval <- c(
      "1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1",
      " TOT. NO. OF INDIVIDUALS:       12",
      " ESTIMATION STEP OMITTED:                 YES ",
      " #OBJV:*******************      742.051       *******"
    )
    .msg <- .failure(.fakeNonmem(.eval), "eval_stop")
    expect_match(.msg, "minimization not successful")
    .msg <- .failure(.fakeNonmem(.eval), "eval_read", list(readBadOpt = TRUE))
    expect_match(.msg, "could not read or use NONMEM's output")
    # nor is a model with every parameter fixed
    .msg <- .failure(
      .fakeNonmem(c(
        .eval[1:2],
        " ESTIMATION STEP IMPLEMENTED BUT THE NUMBER OF PARAMETERS TO BE ESTIMATED IS 0"
      )),
      "allfix_stop"
    )
    expect_match(.msg, "minimization not successful")
  })

  test_that("NONMEM output that is not UTF-8 is still explained (#46)", {
    .dir <- file.path(tempfile(), "enc-nonmem")
    dir.create(.dir, recursive = TRUE)
    writeBin(
      c(
        charToRaw("gfortran: error: C:\\Users\\M"),
        as.raw(0xfc),
        charToRaw("ller\\FSUBS.f90\nLicense file has expired\n")
      ),
      file.path(.dir, "enc.lst")
    )
    .lines <- .nonmemFailureLines(.dir, "enc.lst", "enc.nmctl")
    expect_true(all(validUTF8(.lines)))
    expect_equal(.nonmemClassifyFailure(.lines, "enc")$cause, "license")
    unlink(dirname(.dir), recursive = TRUE)
  })

  test_that("the echoed control stream is not read as a NONMEM message (#46)", {
    .dir <- file.path(tempfile(), "echo-nonmem")
    dir.create(.dir, recursive = TRUE)
    writeLines(
      c("$PK", "  LICENSE_MISSING=1", "; PROGRAM TERMINATED"),
      file.path(.dir, "echo.nmctl")
    )
    writeLines(
      c(
        "$PK",
        "  LICENSE_MISSING=1",
        "; PROGRAM TERMINATED",
        "gfortran: error: cannot compile"
      ),
      file.path(.dir, "echo.lst")
    )
    .lines <- .nonmemFailureLines(.dir, "echo.lst", "echo.nmctl")
    expect_equal(.lines, "gfortran: error: cannot compile")
    expect_null(.nonmemClassifyFailure(.lines))
    unlink(dirname(.dir), recursive = TRUE)
  })

  test_that("errors reading a finished NONMEM run are reported (#46)", {
    .ui <- suppressMessages(.fit(NA, "fail_final"))
    writeLines(
      c(
        "1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1",
        " #TERM:",
        "0MINIMIZATION SUCCESSFUL"
      ),
      file.path(.ui$nonmemExportPath, .ui$nonmemNmlst)
    )
    local_mocked_bindings(.nonmemFinalizeEnv = function(env, oldUi) {
      stop("mock post-processing error", call. = FALSE)
    })
    expect_error(
      .nonmemFinalizeOrExplain(new.env(), .ui),
      "could not read or use NONMEM's output(.|\n)*mock post-processing error"
    )
    local_mocked_bindings(.nonmemFinalizeEnv = function(env, oldUi) "fit")
    expect_equal(.nonmemFinalizeOrExplain(new.env(), .ui), "fit")
  })

  test_that("an earlier failed run does not explain a new failure (#46)", {
    skip_on_os("windows")
    .msg <- .failure(
      .fakeNonmem(c(
        "1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1",
        "0PROGRAM TERMINATED BY OBJ"
      )),
      "fail_stale"
    )
    expect_match(.msg, "NONMEM stopped during the run")
    writeLines(" (DATA ERROR) RECORD 1", file.path("fail_stale-nonmem", "FMSG"))
    .msg <- .failure("babelmixr2-no-such-nmfe", "fail_stale")
    expect_match(.msg, "exit status: 127", fixed = TRUE)
    expect_no_match(.msg, "PROGRAM TERMINATED|DATA ERROR")
  })

  test_that("an NM-TRAN error without NONMEM output is reported (#46)", {
    .nmtranOnly <- function(ctl, directory, ui) {
      writeLines(
        c(
          " (DATA ERROR) RECORD         3, DATA ITEM   6, CONTENTS: 1",
          " ITEM IS OUT OF RANGE."
        ),
        file.path(directory, "FMSG")
      )
    }
    .msg <- .failure(.nmtranOnly, "fail_fmsg")
    expect_match(.msg, "NM-TRAN found an error in the NONMEM data")
    expect_match(.msg, "ITEM IS OUT OF RANGE", fixed = TRUE)
  })

  test_that("finished runs that did not converge can still be read (#46)", {
    .head <- c(
      "1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1",
      " TOT. NO. OF OBS RECS:      132",
      " TOT. NO. OF INDIVIDUALS:       12",
      " #TERM:"
    )
    # a covariance step failure after the termination block
    .cov <- c(
      "0PROGRAM TERMINATED BY OBJ",
      " MESSAGE ISSUED FROM COVARIANCE STEP"
    )
    .cases <- list(
      list(
        name = "read_round_cov",
        ctl = list(readRounding = TRUE),
        lines = c(
          "0MINIMIZATION TERMINATED",
          " DUE TO ROUNDING ERRORS (ERROR=134)",
          " #TERE:",
          .cov
        )
      ),
      list(
        name = "read_round",
        ctl = list(readRounding = TRUE),
        lines = c(
          "0MINIMIZATION TERMINATED",
          " DUE TO ROUNDING ERRORS (ERROR=134)"
        )
      ),
      list(
        name = "read_badopt",
        ctl = list(readBadOpt = TRUE),
        lines = c(
          "0MINIMIZATION TERMINATED",
          " DUE TO MAX. NO. OF FUNCTION EVALUATIONS EXCEEDED"
        )
      ),
      list(
        name = "read_its",
        ctl = list(readBadOpt = TRUE),
        lines = " OPTIMIZATION WAS COMPLETED"
      ),
      list(name = "read_empty", ctl = list(readBadOpt = TRUE), lines = "")
    )
    for (.c in .cases) {
      # the fake output has no estimates, so with the flag the run gets
      # past the failure checks and only reading the estimates fails
      .msg <- .failure(
        .fakeNonmem(c(.head, .c$lines, " #TERE:")),
        .c$name,
        .c$ctl
      )
      expect_match(
        .msg,
        "could not read or use NONMEM's output",
        info = .c$name
      )
      # without the flag the run stops as not successful
      .msg <- .failure(
        .fakeNonmem(c(.head, .c$lines, " #TERE:")),
        paste0(.c$name, "_stop")
      )
      expect_match(.msg, "minimization not successful", info = .c$name)
    }
  })

  test_that("NONMEM exiting abnormally after estimation is reported (#46)", {
    skip_on_os("windows")
    .sh <- normalizePath(tempfile(fileext = ".sh"), mustWork = FALSE)
    writeLines(
      c(
        "printf '%s\\n' '1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.1' \\",
        "  ' TOT. NO. OF INDIVIDUALS:       12' ' #TERM:' '0MINIMIZATION SUCCESSFUL' > \"$2\"",
        "exit 137"
      ),
      .sh
    )
    .msg <- .failure(paste("sh", shQuote(.sh)), "fail_killed")
    expect_match(.msg, "exited with an error after estimation")
    expect_match(.msg, "exit status: 137", fixed = TRUE)
    unlink(.sh)
  })

  test_that("output from a manual run is kept without a run command (#46)", {
    expect_error(suppressMessages(.fit("", "keep_lst")), "runCommand")
    .lst <- file.path("keep_lst-nonmem", "keep_lst.lst")
    writeLines("manual run", .lst)
    expect_error(suppressMessages(.fit("", "keep_lst")), "runCommand")
    expect_true(file.exists(.lst))
  })

  test_that("an unset run command says how to set it (#46)", {
    expect_error(
      suppressMessages(.fit("", "fail_unset")),
      "nonmemControl(runCommand=)",
      fixed = TRUE
    )
  })
})
