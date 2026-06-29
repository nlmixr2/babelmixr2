## Benchmark: nlmer inner nonlinear+Jacobian evaluation
##
## Compares the new approach -- a single multithreaded C solve over the
## resident model (nlmixr2est::nlmerSolveGrad) -- against the previous
## approach of looping rxode2::rxSolve() per subject in R.  Both produce the
## per-observation prediction and d(pred)/d(THETA) Jacobian that lme4::nlmer
## requests on every PWRSS iteration, so this measures the per-iteration cost
## that dominates an nlmer fit.
##
## Run with:  Rscript benchmarks/nlmer-solve-benchmark.R

suppressMessages({
  library(rxode2); library(nlmixr2est); library(microbenchmark)
})
## Load the in-development babelmixr2 (this checkout), not an installed copy.
suppressMessages(pkgload::load_all(".", quiet = TRUE))

oneCmt <- function() {
  ini({ tka <- 0.0; tv <- 2.99573227355399; tcl <- -0.693147180559945
        eta.ka ~ 1.0; eta.v ~ 1.0; eta.cl ~ 1.0; add.sd <- 1.0 })
  model({ ka <- exp(tka + eta.ka); v <- exp(tv + eta.v); cl <- exp(tcl + eta.cl)
          d/dt(depot) <- -depot * ka; d/dt(central) <- depot * ka - cl * central / v
          cp <- central / v; cp ~ add(add.sd) })
}

## Replicate as many subjects as requested by stacking copies of theo_sd
makeData <- function(nrep = 1L) {
  base <- nlmixr2data::theo_sd[!(nlmixr2data::theo_sd$TIME == 0 &
                                   nlmixr2data::theo_sd$EVID == 0), ]
  if (nrep == 1L) return(base)
  nid <- length(unique(base$ID))
  do.call(rbind, lapply(seq_len(nrep) - 1L, function(k) {
    d <- base; d$ID <- d$ID + k * nid; d
  }))
}

benchOne <- function(nrep) {
  ft <- makeData(nrep)
  ui <- rxode2::rxUiDecompress(rxode2::as.rxUi(oneCmt))
  assign("control", nlmerControl(), envir = ui)
  ctl <- ui$control
  ret <- new.env(parent = emptyenv())
  nlmixr2est::.foceiPreProcessData(ft, ret, ui, ctl$rxControl)
  sm <- ui$nlmerSensModel
  start <- ui$nlmerStart
  np <- length(start)

  dsAll <- ret$dataSav[ret$dataSav$EVID != 2, ]
  sid <- unique(dsAll$ID)
  nsub <- length(sid)
  subjEv <- lapply(sid, function(s) dsAll[dsAll$ID == s, ])
  tm <- matrix(rep(unname(start), each = nsub), nsub, np)
  theta1 <- setNames(unname(start), sprintf("THETA[%d]", seq_len(np)))
  rxCtl <- as.list(ctl$rxControl)

  ## ---- old: per-subject R rxSolve loop (benchmarked FIRST, since its
  ##      standalone rxSolve()s reset the rxode2 global solve state that the
  ##      resident solver below depends on) ----
  oldSolve <- function() {
    lapply(subjEv, function(ev) {
      s <- do.call(rxode2::rxSolve,
                   c(list(object = sm$thetaGrad, params = theta1, events = ev),
                     rxCtl))
      list(pred = s$rx_pred_,
           grad = as.matrix(s[, paste0("rx__sens_rx_pred__BY_THETA_",
                                       seq_len(np), "___"), drop = FALSE]))
    })
  }
  mbOld <- microbenchmark(old = oldSolve(), times = 25L)

  ## ---- new: single resident C solve over the whole population ----
  mi <- list(predOnly = sm$predOnly, thetaGrad = sm$thetaGrad,
             eventTheta = sm$eventTheta)
  sc <- babelmixr2:::.nlmerSolveControl(ctl, np)
  suppressMessages(nlmixr2est::.nlmSetupEnv(start, ui, ret$dataSav, mi, sc))
  newSolve <- function() nlmixr2est::nlmerSolveGrad(tm)
  mbNew <- microbenchmark(new = newSolve(), times = 25L)
  nlmixr2est::.nlmFreeEnv()

  oldMs <- median(mbOld$time) / 1e6
  newMs <- median(mbNew$time) / 1e6
  data.frame(subjects = nsub,
             new_ms = round(newMs, 2),
             old_ms = round(oldMs, 2),
             speedup = round(oldMs / newMs, 1))
}

cat(sprintf("rxode2 threads: %d\n\n", rxode2::getRxThreads()))
out <- do.call(rbind, lapply(c(1L, 4L, 16L), benchOne))
print(out, row.names = FALSE)
