.asToolModel <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    wt.cl <- 0.1
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl + wt.cl * WT)
    v <- exp(tv + eta.v)
    linCmt() -> cp
    cp ~ add(add.sd)
  })
}

test_that("as.nonmem() writes the control stream named like nonmemControl(modelName=)", {
  skip_on_cran()
  withr::with_tempdir({
    .f <- suppressMessages(as.nonmem(.asToolModel, modelName = "x"))
    expect_equal(.f, c(nmctl = file.path("x-nonmem", "x.nmctl"),
                       csv = file.path("x-nonmem", "x.csv")))
    expect_true(all(file.exists(.f)))
    .nm <- readLines(.f["nmctl"])
    expect_true(any(grepl("^\\$INPUT .*WT", .nm)))
    # a time-invariant covariate is mu-referenced
    expect_true(any(grepl("MU_2=THETA\\(2\\)\\+WT\\*THETA\\(4\\)", .nm)))
    # the template data has one subject
    .d <- read.csv(.f["csv"])
    expect_equal(unique(.d$ID), 1L)

    # the default name is the model variable
    .f <- suppressMessages(as.nonmem(.asToolModel))
    expect_equal(.f[["nmctl"]], file.path(".asToolModel-nonmem", ".asToolModel.nmctl"))

    # the control's model name is used too
    .f <- suppressMessages(as.nonmem(.asToolModel, control = nonmemControl(modelName = "y")))
    expect_equal(.f[["nmctl"]], file.path("y-nonmem", "y.nmctl"))
  })
})

test_that("as.nonmem() does not mu-reference time-varying covariates", {
  skip_on_cran()
  withr::with_tempdir({
    .f <- suppressMessages(as.nonmem(.asToolModel, modelName = "tv", timeVarying = "WT"))
    .nm <- readLines(.f["nmctl"])
    expect_true(any(grepl("^  MU_2=THETA\\(2\\)$", .nm)))
    expect_true(any(grepl("RX__WT_CL\\*WT", .nm)))
    .d <- read.csv(.f["csv"])
    expect_gt(length(unique(.d$WT)), 1L)
    # covariates are matched ignoring case
    .f2 <- suppressMessages(as.nonmem(.asToolModel, modelName = "tv2", timeVarying = "wt"))
    expect_equal(gsub("tv2", "tv", readLines(.f2["nmctl"]))[-1], .nm[-1])
  })
})

test_that("as.nonmem() data options add the data items", {
  skip_on_cran()
  withr::with_tempdir({
    .f <- suppressMessages(as.nonmem(.asToolModel, modelName = "inf", doseType = "infusion", ss = TRUE))
    .nm <- readLines(.f["nmctl"])
    .input <- .nm[grepl("^\\$INPUT", .nm)]
    expect_true(grepl("\\bRATE\\b", .input))
    expect_true(grepl("\\bSS\\b", .input))
    expect_true(grepl("\\bII\\b", .input))

    .f <- suppressMessages(as.nonmem(.asToolModel, modelName = "cens", cens = TRUE))
    .nm <- readLines(.f["nmctl"])
    expect_true(grepl("\\bCENS\\b", .nm[grepl("^\\$INPUT", .nm)]))
  })
})

test_that("as.nonmem() can translate with data", {
  skip_on_cran()
  .d <- nlmixr2data::theo_sd
  .d$WT <- 70
  withr::with_tempdir({
    .f <- suppressMessages(as.nonmem(.asToolModel, data = .d, modelName = "d"))
    suppressMessages(nlmixr2(.asToolModel, .d, "nonmem",
                             nonmemControl(runCommand = NA, modelName = "e")))
    expect_equal(gsub("\\be\\b", "d", readLines(file.path("e-nonmem", "e.nmctl"))),
                 readLines(.f["nmctl"]))
    expect_equal(nrow(read.csv(.f["csv"])), nrow(.d))
    expect_warning(suppressMessages(as.nonmem(.asToolModel, data = .d, modelName = "w",
                                              timeVarying = "WT")),
                   "ignored")
  })
})

test_that("as.monolix() writes the files named like monolixControl(modelName=)", {
  skip_on_cran()
  withr::with_tempdir({
    .f <- suppressMessages(as.monolix(.asToolModel, modelName = "x"))
    expect_equal(.f, c(mlxtran = "x-monolix.mlxtran", model = "x-monolix.txt",
                       csv = "x-monolix.csv"))
    expect_true(all(file.exists(.f)))
    .mlx <- readLines(.f["mlxtran"])
    expect_true(any(grepl("WT = \\{use=covariate", .mlx)))
    expect_true(any(grepl("covariate = \\{WT\\}", .mlx)))

    .f <- suppressMessages(as.monolix(.asToolModel, modelName = "tv", timeVarying = "WT"))
    .mlx <- readLines(.f["mlxtran"])
    expect_true(any(grepl("WT = \\{use=regressor\\}", .mlx)))
    expect_false(any(grepl("covariate = \\{WT\\}", .mlx)))

    .f <- suppressMessages(as.monolix(.asToolModel, modelName = "inf", doseType = "infusion"))
    .mlx <- readLines(.f["mlxtran"])
    expect_true(any(grepl("use=rate|use=infusiontime", .mlx)))
  })
})

test_that("as.nonmem()/as.monolix() check their arguments", {
  skip_on_cran()
  withr::with_tempdir({
    expect_error(as.nonmem(.asToolModel, timeVarying = "AGE"), "not covariates")
    expect_error(as.monolix(.asToolModel, dose = "periph"), "not in the model")
    expect_error(as.nonmem(.asToolModel, dose = 5))
    expect_error(as.nonmem(.asToolModel, doseType = "oral"))
    expect_error(as.nonmem(.asToolModel, control = monolixControl()))
  })
})
