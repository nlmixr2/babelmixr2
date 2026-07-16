test_that("saemix continuous PK model comparison", {
  skip_on_cran()
  # 1. Direct saemix run on theo.saemix
  data(theo.saemix, package = "saemix")
  saemixData <- saemix::saemixData(name.data = theo.saemix, header = TRUE, sep = " ", na = NA,
                                   name.group = c("Id"), name.predictors = c("Dose", "Time"),
                                   name.response = c("Concentration"),
                                   units = list(x = "hr", y = "mg/L"), name.X = "Time",
                                   verbose = FALSE)

  model1Cpt <- function(psi, id, xidep) {
    dose <- xidep[, 1]
    tim <- xidep[, 2]
    ka <- psi[id, 1]
    V <- psi[id, 2]
    CL <- psi[id, 3]
    k <- CL/V
    yPred <- dose * ka / (V * (ka - k)) * (exp(-k * tim) - exp(-ka * tim))
    return(yPred)
  }

  saemixModel <- saemix::saemixModel(model = model1Cpt,
                                     psi0 = matrix(c(1.0, 20, 0.5), ncol = 3, byrow = TRUE,
                                                   dimnames = list(NULL, c("ka", "V", "CL"))),
                                     transform.par = c(1, 1, 1),
                                     error.model = "constant",
                                     verbose = FALSE)

  saemixOptions <- list(seed = 632545, save = FALSE, save.graphs = FALSE, print = FALSE,
                        nbiter.saemix = c(10, 5), fim = FALSE)

  fitDirect <- saemix::saemix(saemixModel, saemixData, saemixOptions)

  # 2. nlmixr2 saemix run
  oneCmt <- function() {
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

  filteredTheo <- nlmixr2data::theo_sd[!(nlmixr2data::theo_sd$TIME == 0 & nlmixr2data::theo_sd$EVID == 0), ]
  fitNlmixr <- nlmixr2(oneCmt, filteredTheo, est = "saemix",
                       saemixControl(seed = 632545, nbiter.saemix = c(10, 5),
                                     fim = FALSE, warnings = FALSE))

  # Compare population parameter estimates (fixed effects)
  directEst <- fitDirect@results@fixed.effects
  nlmixrEst <- exp(fitNlmixr$theta[c("tka", "tv", "tcl")])

  expect_equal(as.numeric(nlmixrEst), as.numeric(directEst), tolerance = 1e-3)

  # Compare residual error estimates
  expect_equal(as.numeric(fitNlmixr$theta["add.sd"]), as.numeric(fitDirect@results@respar[1]), tolerance = 1e-3)

  # Compare random effects variances (omega diagonals)
  directOmega <- diag(fitDirect@results@omega)
  nlmixrOmega <- diag(fitNlmixr$omega)
  expect_equal(as.numeric(nlmixrOmega[c("eta.ka", "eta.v", "eta.cl")]), as.numeric(directOmega), tolerance = 1e-3)

  # Verify saemix model is embedded in the fit environment
  expect_s4_class(fitNlmixr$saemix, "SaemixObject")
})

test_that("saemix discrete likelihood model comparison", {
  # 1. Direct saemix run on toenail.saemix
  skip_on_cran()
  data(toenail.saemix, package = "saemix")
  saemixData <- saemix::saemixData(name.data = toenail.saemix, name.group = c("id"), name.predictors = c("time", "y"),
                                   name.response = "y", name.X = c("time"),
                                   verbose = FALSE)

  binaryModel <- function(psi, id, xidep) {
    tim <- xidep[, 1]
    y <- xidep[, 2]
    inter <- psi[id, 1]
    slope <- psi[id, 2]
    logit <- inter + slope * tim
    pevent <- exp(logit) / (1 + exp(logit))
    pObs = (y == 0) * (1 - pevent) + (y == 1) * pevent
    logpdf <- log(pObs)
    return(logpdf)
  }

  saemixModel <- saemix::saemixModel(model = binaryModel, description = "Binary model",
                                     modeltype = "likelihood",
                                     psi0 = matrix(c(-5, -0.1), ncol = 2, byrow = TRUE,
                                                   dimnames = list(NULL, c("theta1", "theta2"))),
                                     transform.par = c(0, 0),
                                     covariance.model = matrix(c(1, 0, 0, 1), ncol = 2),
                                     omega.init = matrix(c(1, 0, 0, 1), ncol = 2),
                                     verbose = FALSE)

  saemixOptions <- list(seed = 1234567, save = FALSE, save.graphs = FALSE, print = FALSE,
                        nbiter.saemix = c(10, 5), fim = FALSE, nb.chains = 1)

  fitDirect <- saemix::saemix(saemixModel, saemixData, saemixOptions)

  # 2. nlmixr2 saemix run
  binaryModRx <- function() {
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

  toenailDf <- toenail.saemix
  colnames(toenailDf)[colnames(toenailDf) == "id"] <- "ID"
  colnames(toenailDf)[colnames(toenailDf) == "time"] <- "TIME"
  colnames(toenailDf)[colnames(toenailDf) == "y"] <- "DV"
  toenailDf$EVID <- 0

  fitNlmixr <- nlmixr2(binaryModRx, toenailDf, est = "saemix",
                       saemixControl(seed = 1234567, nbiter.saemix = c(10, 5),
                                     nb.chains = 1, fim = FALSE, warnings = FALSE))

  # Compare estimates
  directEst <- fitDirect@results@fixed.effects
  nlmixrEst <- fitNlmixr$theta[c("theta1", "theta2")]

  expect_equal(as.numeric(nlmixrEst), as.numeric(directEst), tolerance = 1e-3)
})
