test_that("monolixControl(prune=) option", {
  expect_equal(monolixControl()$prune, "auto")
  expect_true(monolixControl(prune = TRUE)$prune)
  expect_false(monolixControl(prune = FALSE)$prune)
  expect_error(monolixControl(prune = "a"), "prune")
  expect_error(monolixControl(prune = NA), "prune")
})

test_that("Monolix logical expressions as numbers are indicators (#11)", {
  one.cmt <- function() {
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
  ui <- rxode2::rxode2(one.cmt)
  expect_equal(
    rxToMonolix("z <- (wt > 1)*k1 + (1-(wt > 1))*k2", ui),
    paste(
      c(
        "  rx_l001 = 0",
        "  if wt>1",
        "    rx_l001 = 1",
        "  end",
        "   z = (rx_l001)*k1+(1-(rx_l001))*k2"
      ),
      collapse = "\n"
    )
  )
  expect_equal(
    rxToMonolix("z <- !(wt > 1 && sex <= 2)", ui),
    paste(
      c(
        "  rx_l001 = 0",
        "  if ~(wt>1&&sex<=2)",
        "    rx_l001 = 1",
        "  end",
        "   z = rx_l001"
      ),
      collapse = "\n"
    )
  )
  # conditions are still written directly
  expect_equal(
    rxToMonolix("if (wt > 1 && !(sex == 1)) {z <- 1}", ui),
    "  if wt>1&&~(sex==1)\n   z = 1\n  end\n"
  )
})

test_that("prune=\"auto\" only prunes Monolix models when needed (#11)", {
  .needs <- function(expr) {
    .bblNeedsPrune(list(lstExpr = list(expr)), nested = TRUE)
  }
  # Monolix writes if/elseif/else (and nested if) statements
  expect_false(.needs(quote(
    if (a > 1) {
      b <- 1
    }
  )))
  expect_false(.needs(quote(
    if (a > 1) {
      b <- 1
    } else if (a > 0) {
      b <- 2
    } else {
      b <- 3
    }
  )))
  expect_false(.needs(quote(
    if (a > 1) {
      if (c > 1) {
        b <- 1
      }
    }
  )))
  # but not ifelse()
  expect_true(.needs(quote(b <- ifelse(a > 1, 1, 2))))
  expect_true(.needs(quote(
    if (a > 1) {
      b <- 1
    } else {
      b <- ifelse(a > 0, 1, 2)
    }
  )))

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
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl2 / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }
  fIfelse <- rxode2::rxode2(f)
  fIfelse <- rxode2::model(fIfelse, fd <- ifelse(SEX == 1, 0.8, 1), append = v)
  fIfelse <- rxode2::model(fIfelse, cp <- central / v * fd)
  d <- nlmixr2data::theo_sd
  d$WT <- ifelse(d$ID %% 2 == 0, 80, 60)
  d$SEX <- as.integer(d$ID %% 3 == 0)
  .mod <- function(name) {
    readLines(paste0(name, "-monolix.txt"))
  }
  withr::with_tempdir({
    # nested if/else is written as Monolix if blocks (not pruned)
    suppressMessages(expect_error(
      nlmixr2(
        f,
        d,
        "monolix",
        monolixControl(runCommand = NA, modelName = "nested")
      ),
      NA
    ))
    .m <- .mod("nested")
    expect_true(any(grepl("^ *else$", .m)))
    expect_true(any(grepl("^ *if WT>50$", .m)))
    expect_false(any(grepl("rx_l", .m)))

    # prune=TRUE prunes it
    suppressMessages(expect_error(
      nlmixr2(
        f,
        d,
        "monolix",
        monolixControl(runCommand = NA, modelName = "pruned", prune = TRUE)
      ),
      NA
    ))
    .m <- .mod("pruned")
    expect_false(any(grepl("else", .m)))
    expect_true(any(grepl("^ *if WT>70$", .m)))
    expect_true(any(grepl("^ *rx_l[0-9]+ = 1$", .m)))

    # ifelse() cannot be written without pruning
    expect_error(
      suppressMessages(nlmixr2(
        fIfelse,
        d,
        "monolix",
        monolixControl(runCommand = NA, modelName = "noprune", prune = FALSE)
      )),
      "ifelse"
    )
    # and is pruned by default
    suppressMessages(expect_error(
      nlmixr2(
        fIfelse,
        d,
        "monolix",
        monolixControl(runCommand = NA, modelName = "fdmodel")
      ),
      NA
    ))
    .m <- .mod("fdmodel")
    expect_false(any(grepl("ifelse|else", .m)))
    expect_true(any(grepl("^ *fd = .*rx_l[0-9]+", .m)))
  })
})

test_that("Monolix compartment properties write their indicators first (#11)", {
  one.cmt <- function() {
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
  ui <- rxode2::rxode2(one.cmt)
  expect_equal(
    rxToMonolix("f(depot) <- (sex == 1)*0.8 + (1 - (sex == 1))", ui),
    paste(
      c(
        "  ;f defined in PK section",
        "  rx_l001 = 0",
        "  if sex==1",
        "    rx_l001 = 1",
        "  end",
        "   rx_f_depot = (rx_l001)*0.8+(1-(rx_l001))"
      ),
      collapse = "\n"
    )
  )
})
