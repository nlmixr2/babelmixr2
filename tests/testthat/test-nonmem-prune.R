test_that("nonmemControl(prune=) option", {
  expect_false(nonmemControl()$prune)
  expect_true(nonmemControl(prune=TRUE)$prune)
  expect_error(nonmemControl(prune="a"))
  expect_error(nonmemControl(prune=NA))
})

test_that("NONMEM logical expressions used as numbers use indicator variables (#11)", {
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
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl/v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }
  ui <- rxode2::rxode2(one.cmt)
  withr::with_options(list(babelmixr2.protectZeros=FALSE), {
    expect_equal(rxToNonmem("z <- (wt > 1)*k1 + (1-(wt > 1))*k2", ui),
                 paste(c("  RXL001=0",
                         "  IF (WT.GT.1) RXL001=1",
                         "  Z=(RXL001)*K1+(1-(RXL001))*K2 ; z <- (wt > 1) * k1 + (1 - (wt > 1)) * k2"),
                       collapse="\n"))
    expect_equal(rxToNonmem("z <- !(wt > 1 && sex <= 2)", ui),
                 paste(c("  RXL001=0",
                         "  IF (.NOT. ((WT.GT.1.AND.SEX.LE.2))) RXL001=1",
                         "  Z=RXL001 ; z <- !(wt > 1 && sex <= 2)"),
                       collapse="\n"))
    # a numeric condition is true when it is not zero
    expect_equal(rxToNonmem("if (wt) {z=1}", ui),
                 paste(c("  IF ((WT).NE.0) THEN",
                         "    Z=1 ; z = 1",
                         "  END IF",
                         ""),
                       collapse="\n"))
  })
})

test_that("NONMEM models with nested if/else can be pruned (#11)", {
  f <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
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
      if (WT > 70) {
        if (SEX == 1) {
          cl2 <- cl * 1.2
        } else {
          cl2 <- cl * 1.1
        }
      } else if (WT > 50) {
        cl2 <- cl
      } else {
        cl2 <- cl * 0.8
      }
      if (SEX == 1) {
        fd <- 0.8
      } else {
        fd <- 1
      }
      d/dt(depot) <- -ka * depot
      f(depot) <- fd
      d/dt(central) <- ka * depot - cl2 / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }
  d <- nlmixr2data::theo_sd
  d$WT <- ifelse(d$ID %% 2 == 0, 80, 60)
  d$SEX <- as.integer(d$ID %% 3 == 0)

  withr::with_tempdir({
    expect_error(
      suppressMessages(nlmixr2(f, d, "nonmem",
                               nonmemControl(runCommand=NA, modelName="noprune"))),
      "nonmemControl\\(prune=TRUE\\)")

    suppressMessages(
      expect_error(nlmixr2(f, d, "nonmem",
                           nonmemControl(runCommand=NA, modelName="prune",
                                         prune=TRUE)),
                   NA))
    .ctl <- readLines(file.path("prune-nonmem", "prune.nmctl"))
    # no ELSE statements are written
    expect_false(any(grepl("ELSE", .ctl)))
    # the conditions are written as indicator variables
    expect_true(any(grepl("^ +IF \\(WT\\.GT\\.70\\) RXL[0-9]+=1$", .ctl)))
    expect_true(any(grepl("^ +IF \\(SEX\\.EQ\\.1\\) RXL[0-9]+=1$", .ctl)))
    # logical expressions are not used as numbers
    .noComment <- sub(";.*$", "", .ctl)
    .noIf <- .noComment[!grepl("^ *IF \\(", .noComment)]
    expect_false(any(grepl("\\.(GT|EQ|LT|GE|LE|NE)\\.", .noIf)))
    # the bioavailability indicator is defined in $PK before F1
    .pk <- which(.ctl == "$PK")
    .des <- which(.ctl == "$DES")
    .f1 <- grep("^  F1=", .ctl)
    expect_length(.f1, 1L)
    expect_true(.f1 > .pk && .f1 < .des)
    .rxl <- regmatches(.ctl[.f1], regexpr("RXL[0-9]+", .ctl[.f1]))
    .def <- grep(paste0("^  ", .rxl, "=0$"), .ctl)
    expect_length(.def, 1L)
    expect_true(.def > .pk && .def < .f1)
  })
})

test_that("models without if/else are not changed by pruning (#11)", {
  f <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
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
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }
  .env <- new.env(parent=emptyenv())
  .env$ui <- rxode2::rxode2(f)
  .ui <- .env$ui
  .bblPruneIf(.env, "NONMEM")
  expect_identical(.env$ui, .ui)
})
