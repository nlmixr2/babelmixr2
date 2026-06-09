test_that("saemix continuous PK model comparison", {
  skip_on_cran()
  # 1. Direct saemix run on theo.saemix
  data(theo.saemix, package = "saemix")
  saemix.data <- saemix::saemixData(name.data = theo.saemix, header = TRUE, sep = " ", na = NA,
                                    name.group = c("Id"), name.predictors = c("Dose", "Time"),
                                    name.response = c("Concentration"),
                                    units = list(x = "hr", y = "mg/L"), name.X = "Time",
                                    verbose = FALSE)

  model1cpt <- function(psi, id, xidep) {
    dose <- xidep[, 1]
    tim <- xidep[, 2]
    ka <- psi[id, 1]
    V <- psi[id, 2]
    CL <- psi[id, 3]
    k <- CL/V
    ypred <- dose * ka / (V * (ka - k)) * (exp(-k * tim) - exp(-ka * tim))
    return(ypred)
  }

  saemix.model <- saemix::saemixModel(model = model1cpt,
                                      psi0 = matrix(c(1.0, 20, 0.5), ncol = 3, byrow = TRUE,
                                                    dimnames = list(NULL, c("ka", "V", "CL"))),
                                      transform.par = c(1, 1, 1),
                                      error.model = "constant",
                                      verbose = FALSE)

  saemix.options <- list(seed = 632545, save = FALSE, save.graphs = FALSE, print = FALSE,
                         nbiter.saemix = c(10, 5), fim = FALSE)

  fit_direct <- saemix::saemix(saemix.model, saemix.data, saemix.options)

  # 2. nlmixr2 saemix run
  one.cmt <- function() {
    ini({
      tka <- 0.0
      tv <- 2.99573227355399
      tcl <- -0.693147180559945
      eta.ka ~ 1.0
      eta.v ~ 1.0
      eta.cl ~ 1.0
      add.sd <- 1.0
    })
    model({
      ka <- exp(tka + eta.ka)
      v <- exp(tv + eta.v)
      cl <- exp(tcl + eta.cl)
      d/dt(depot) <- -depot * ka
      d/dt(central) <- depot * ka - cl * central / v
      cp <- central / v
      cp ~ add(add.sd)
    })
  }

  filtered_theo <- nlmixr2data::theo_sd[!(nlmixr2data::theo_sd$TIME == 0 & nlmixr2data::theo_sd$EVID == 0), ]
  fit_nlmixr <- nlmixr2(one.cmt, filtered_theo, est = "saemix",
                        saemixControl(seed = 632545, nbiter.saemix = c(10, 5),
                                      fim = FALSE, warnings = FALSE))

  # Compare population parameter estimates (fixed effects)
  direct_est <- fit_direct@results@fixed.effects
  nlmixr_est <- exp(fit_nlmixr$theta[c("tka", "tv", "tcl")])

  expect_equal(as.numeric(nlmixr_est), as.numeric(direct_est), tolerance = 1e-3)

  # Compare residual error estimates
  expect_equal(as.numeric(fit_nlmixr$theta["add.sd"]), as.numeric(fit_direct@results@respar[1]), tolerance = 1e-3)

  # Compare random effects variances (omega diagonals)
  direct_omega <- diag(fit_direct@results@omega)
  nlmixr_omega <- diag(fit_nlmixr$omega)
  expect_equal(as.numeric(nlmixr_omega[c("eta.ka", "eta.v", "eta.cl")]), as.numeric(direct_omega), tolerance = 1e-3)

  # Verify saemix model is embedded in the fit environment
  expect_s4_class(fit_nlmixr$saemix, "SaemixObject")
})

test_that("saemix discrete likelihood model comparison", {
  # 1. Direct saemix run on toenail.saemix
  skip_on_cran()
  data(toenail.saemix, package = "saemix")
  saemix.data <- saemix::saemixData(name.data = toenail.saemix, name.group = c("id"), name.predictors = c("time", "y"),
                                    name.response = "y", name.X = c("time"),
                                    verbose = FALSE)

  binary.model <- function(psi, id, xidep) {
    tim <- xidep[, 1]
    y <- xidep[, 2]
    inter <- psi[id, 1]
    slope <- psi[id, 2]
    logit <- inter + slope * tim
    pevent <- exp(logit) / (1 + exp(logit))
    pobs = (y == 0) * (1 - pevent) + (y == 1) * pevent
    logpdf <- log(pobs)
    return(logpdf)
  }

  saemix.model <- saemix::saemixModel(model = binary.model, description = "Binary model",
                                      modeltype = "likelihood",
                                      psi0 = matrix(c(-5, -0.1), ncol = 2, byrow = TRUE,
                                                    dimnames = list(NULL, c("theta1", "theta2"))),
                                      transform.par = c(0, 0),
                                      covariance.model = matrix(c(1, 0, 0, 1), ncol = 2),
                                      omega.init = matrix(c(1, 0, 0, 1), ncol = 2),
                                      verbose = FALSE)

  saemix.options <- list(seed = 1234567, save = FALSE, save.graphs = FALSE, print = FALSE,
                         nbiter.saemix = c(10, 5), fim = FALSE, nb.chains = 1)

  fit_direct <- saemix::saemix(saemix.model, saemix.data, saemix.options)

  # 2. nlmixr2 saemix run
  binary_mod_rx <- function() {
    ini({
      theta1 <- -5.0
      theta2 <- -0.1
      eta.theta1 ~ 1.0
      eta.theta2 ~ 1.0
    })
    model({
      inter <- theta1 + eta.theta1
      slope <- theta2 + eta.theta2
      logit <- inter + slope * TIME
      pevent <- exp(logit) / (1 + exp(logit))
      ll(bin) ~ (DV == 0) * log(1 - pevent) + (DV == 1) * log(pevent)
    })
  }

  toenail_df <- toenail.saemix
  colnames(toenail_df)[colnames(toenail_df) == "id"] <- "ID"
  colnames(toenail_df)[colnames(toenail_df) == "time"] <- "TIME"
  colnames(toenail_df)[colnames(toenail_df) == "y"] <- "DV"
  toenail_df$EVID <- 0

  fit_nlmixr <- nlmixr2(binary_mod_rx, toenail_df, est = "saemix",
                        saemixControl(seed = 1234567, nbiter.saemix = c(10, 5),
                                      nb.chains = 1, fim = FALSE, warnings = FALSE))

  # Compare estimates
  direct_est <- fit_direct@results@fixed.effects
  nlmixr_est <- fit_nlmixr$theta[c("theta1", "theta2")]

  expect_equal(as.numeric(nlmixr_est), as.numeric(direct_est), tolerance = 1e-3)
})
