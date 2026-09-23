.monolixPriorModel <- function() {
  ini({
    tka <- 0.45
    tcl <- log(c(0, 2.7, 100))
    tv <- 3.45
    temax <- logit(0.5)
    cl.wt <- 0
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl + WT * cl.wt)
    v <- exp(tv + eta.v)
    emax <- expit(temax)
    d/dt(depot) <- -depot*ka
    d/dt(central) <- depot*ka - cl*central/v
    cp <- emax * central / v
    cp ~ add(add.sd)
  })
}

# the model, or a skip when this lotri/rxode2 cannot store priors
.monolixPriorUi <- function() {
  .u <- rxode2::rxode2(.monolixPriorModel)
  .p <- try(ini(.u, prior(tka) ~ dnorm(0, 1)), silent=TRUE)
  if (inherits(.p, "try-error") || !any(names(.p$iniDf) == "prior")) {
    skip("this lotri/rxode2 does not support priors")
  }
  .u
}

.mlxLines <- function(txt, pattern) {
  .l <- strsplit(txt, "\n")[[1]]
  .l[grepl(pattern, .l)]
}

test_that("no priors, no [POPULATION] section and no MAP", {
  .u <- rxode2::rxode2(.monolixPriorModel)
  expect_equal(.u$mlxtranModelPopulation, "")
  expect_false(grepl("[POPULATION]", .u$mlxtranModel, fixed=TRUE))
  expect_false(grepl("MAP", .u$mlxtranParameter, fixed=TRUE))
})

test_that("normal priors become Monolix MAP estimation", {
  .u <- ini(.monolixPriorUi(), prior(tka) ~ dnorm(log(1.5), 0.5))
  .u <- ini(.u, prior(temax) ~ dnorm(0, 2))
  .u <- ini(.u, prior(cl.wt) ~ stdNormal())
  .u <- ini(.u, prior(add.sd) ~ dnorm(0.7, 0.1))

  .par <- .u$mlxtranParameter
  expect_equal(.mlxLines(.par, "^ka_pop="),
               paste0("ka_pop={value=", exp(0.45), ", method=MAP}"))
  expect_equal(.mlxLines(.par, "^emax_pop="),
               "emax_pop={value=0.5, method=MAP}")
  expect_equal(.mlxLines(.par, "^beta_cl_WT="),
               "beta_cl_WT={value=0, method=MAP}")
  expect_equal(.mlxLines(.par, "^add__sd="),
               "add__sd={value=0.7, method=MAP}")
  # parameters without a prior are still maximum likelihood
  expect_equal(.mlxLines(.par, "^v_pop="),
               paste0("v_pop={value=", exp(3.45), ", method=MLE}"))
  expect_equal(.mlxLines(.par, "^omega_ka="),
               paste0("omega_ka={value=", sqrt(0.6), ", method=MLE}"))

  # the prior mean is back-transformed like the estimate; the sd is the
  # Gaussian-space sd Monolix uses for its prior distributions
  .pop <- .u$mlxtranModelPopulation
  expect_equal(.mlxLines(.pop, "^ka_pop "),
               "ka_pop = {distribution=logNormal, typical=1.5, sd=0.5}")
  expect_equal(.mlxLines(.pop, "^emax_pop "),
               "emax_pop = {distribution=logitNormal, min=0, max=1, typical=0.5, sd=2}")
  expect_equal(.mlxLines(.pop, "^beta_cl_WT "),
               "beta_cl_WT = {distribution=normal, typical=0, sd=1}")
  expect_equal(.mlxLines(.pop, "^add__sd "),
               "add__sd = {distribution=normal, typical=0.7, sd=0.1}")

  # [POPULATION] comes before [INDIVIDUAL] in <MODEL>
  .mod <- .u$mlxtranModel
  expect_true(regexpr("[POPULATION]", .mod, fixed=TRUE) <
                regexpr("[INDIVIDUAL]", .mod, fixed=TRUE))
})

test_that("a normal prior written with named arguments is read the same", {
  .u <- ini(.monolixPriorUi(), prior(tka) ~ dnorm(mean=log(1.5), sd=0.5))
  expect_equal(.mlxLines(.u$mlxtranModelPopulation, "^ka_pop "),
               "ka_pop = {distribution=logNormal, typical=1.5, sd=0.5}")
})

test_that("a fixed parameter keeps FIXED and has no prior definition", {
  .u <- ini(.monolixPriorUi(), prior(tka) ~ dnorm(log(1.5), 0.5))
  .u <- ini(.u, tka=fix(0.45))
  expect_equal(.mlxLines(.u$mlxtranParameter, "^ka_pop="),
               paste0("ka_pop={value=", exp(0.45), ", method=FIXED}"))
  expect_equal(.u$mlxtranModelPopulation, "")
})

test_that("priors Monolix cannot represent are refused", {
  .u <- ini(.monolixPriorUi(), prior(om.eta.cl) ~ dnorm(0.3, 0.1))
  expect_error(.u$mlxtranParameter, "prior on the omega")

  # a block prior is stored on the block's first diagonal element
  .u <- ini(.monolixPriorUi(), eta.ka + eta.cl ~ c(0.6, 0.01, 0.3))
  .u <- ini(.u, prior(eta.ka, eta.cl) ~ invWishart(4))
  expect_error(.u$mlxtranParameter, "prior on the omega")

  # Monolix's probitNormal is only on (0, 1)
  .p <- model(rxode2::rxode2(.monolixPriorModel), emax <- probitInv(temax, 0, 100))
  .p <- ini(.p, temax=0)
  expect_error(ini(.p, prior(temax) ~ dnorm(0, 1))$mlxtranParameter, "probitInv")
  .p01 <- model(.p, emax <- probitInv(temax))
  expect_equal(.mlxLines(ini(.p01, prior(temax) ~ dnorm(0, 1))$mlxtranModelPopulation,
                         "^emax_pop "),
               "emax_pop = {distribution=probitNormal, typical=0.5, sd=1}")

  .u <- ini(.monolixPriorUi(), prior(tcl) ~ dcauchy(0, 1))
  expect_error(.u$mlxtranParameter, "univariate normal prior")

  .u <- ini(.monolixPriorUi(), prior(tcl, tv) ~ multiNormal(c(1, 3), lotri(tcl + tv ~ c(1, 0.1, 1))))
  expect_error(.u$mlxtranParameter, "univariate normal prior")

  .u <- ini(.monolixPriorUi(), prior(add.sd) ~ dnorm(0.7, 0.1))
  .u2 <- try(model(.u, cp ~ add(add.sd) + var()), silent=TRUE)
  if (!inherits(.u2, "try-error")) {
    expect_error(.u2$mlxtranParameter, "on a variance")
  }
})

test_that("nlmixr2 writes a MAP mlxtran for a model with a normal prior", {
  .u <- ini(.monolixPriorUi(), prior(tka) ~ dnorm(log(1.5), 0.5))
  withr::with_tempdir({
    suppressMessages(try(
      nlmixr2(.u, nlmixr2data::theo_sd, "monolix",
              monolixControl(runCommand=NA, modelName="monolixPrior")),
      silent=TRUE))
    .f <- "monolixPrior-monolix.mlxtran"
    expect_true(file.exists(.f))
    .txt <- paste(readLines(.f), collapse="\n")
    expect_equal(.mlxLines(.txt, "^ka_pop="),
                 paste0("ka_pop={value=", exp(0.45), ", method=MAP}"))
    expect_equal(.mlxLines(.txt, "^ka_pop "),
                 "ka_pop = {distribution=logNormal, typical=1.5, sd=0.5}")
  })
})

test_that("nlmixr2 refuses an omega prior for monolix before writing files", {
  .u <- ini(.monolixPriorUi(), prior(om.eta.cl) ~ dnorm(0.3, 0.1))
  withr::with_tempdir({
    expect_error(
      nlmixr2(.u, nlmixr2data::theo_sd, "monolix",
              monolixControl(runCommand=NA, modelName="monolixPrior")),
      "omega")
    expect_false(file.exists("monolixPrior-monolix.mlxtran"))
  })
})
