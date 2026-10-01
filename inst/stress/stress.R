# babelmixr2 NONMEM/Monolix stress test cases and runner
#
# This file is shared by the testthat stress tests (translation only)
# and by `run-stress.R` (which can also run NONMEM/Monolix end to end).
# See README.md in this directory for how to run it; from an R session
# (like RStudio) source this file and call stressCheck() and
# stressKit().
#
# Every case is a model, a data set and, for each engine, what should
# happen:
#
#  - "ok": the model and data are translated (and fit in "run" mode)
#  - "any": translated, or refused with a known error
#  - a regular expression: the model is refused with an error matching it
#
# Cases can also list regular expressions that must be found in the
# NONMEM control stream (`nonmem`) or the Monolix model file
# (`monolix`), which check that the edge case was translated the way
# it should be.

.stressEnv <- new.env(parent = emptyenv())

# ---------------------------------------------------------------------
# data sets
# ---------------------------------------------------------------------

.stressTimes <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12, 24)

#' Simulate a data set from a model and an event table
#'
#' @param model model function
#' @param ev event table (or data frame) with named observation
#'   compartments for multiple endpoint models
#' @param seed random seed
#' @return NONMEM-like data frame with `ID`, `TIME`, `EVID`, `AMT`,
#'   `CMT`, `DV` and any other event table columns
#' @noRd
.stressSim <- function(model, ev, seed=42) {
  .d <- as.data.frame(ev)
  .d <- .d[order(.d$id, .d$time, -.d$evid), ]
  .s <- rxode2::rxWithSeed(seed, suppressMessages(
    rxode2::rxSolve(model, .d, addDosing=TRUE, returnType="data.frame")))
  .obs <- which(.d$evid == 0)
  .sobs <- which(.s$evid == 0)
  stopifnot(length(.obs) == length(.sobs),
            isTRUE(all.equal(.d$time[.obs], .s$time[.sobs])))
  .d$dv <- NA_real_
  .d$dv[.obs] <- .s$sim[.sobs]
  names(.d) <- toupper(names(.d))
  .keep <- vapply(names(.d), function(n) !all(is.na(.d[[n]])), logical(1))
  .keep[c("ID", "TIME", "EVID", "AMT", "CMT", "DV")] <- TRUE
  .d <- .d[, .keep, drop=FALSE]
  rownames(.d) <- NULL
  .d
}

.stressEvOral <- function(n=24, amt=320, times=.stressTimes, cmtObs=NULL) {
  .ev <- rxode2::et(amt=amt, cmt="depot")
  if (is.null(cmtObs)) {
    .ev <- rxode2::et(.ev, times)
  } else {
    for (.c in cmtObs) .ev <- rxode2::et(.ev, times, cmt=.c)
  }
  rxode2::et(.ev, id=seq_len(n))
}

.stressEvIv <- function(n=24, amt=500, times=.stressTimes, rate=NULL, dur=NULL) {
  .ev <- if (!is.null(rate)) {
    rxode2::et(amt=amt, cmt="central", rate=rate)
  } else if (!is.null(dur)) {
    rxode2::et(amt=amt, cmt="central", dur=dur)
  } else {
    rxode2::et(amt=amt, cmt="central")
  }
  .ev <- rxode2::et(.ev, times)
  rxode2::et(.ev, id=seq_len(n))
}

.stressTheo <- function() {
  .d <- nlmixr2data::theo_sd
  # theo_sd uses evid=101 (dose into compartment 1)
  .d$EVID <- ifelse(.d$EVID == 0, 0L, 1L)
  .d
}

# ---------------------------------------------------------------------
# models
# ---------------------------------------------------------------------

.stressModels <- list(
  lin1oral=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  },
  lin1oralLagF=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; tlag <- log(0.2); tfdepot <- fix(0.8)
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      prop.sd <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      lagD <- exp(tlag)
      alag(depot) <- lagD
      f(depot) <- tfdepot
      cp <- linCmt()
      cp ~ prop(prop.sd)
    })
  },
  lin1iv=function() {
    ini({
      tcl <- log(4); tv <- log(40)
      eta.cl ~ 0.1; eta.v ~ 0.1
      prop.sd <- 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      cp <- linCmt()
      cp ~ prop(prop.sd)
    })
  },
  lin1ivK=function() {
    ini({
      tk <- log(0.1); tv <- log(40)
      eta.k ~ 0.1; eta.v ~ 0.1
      prop.sd <- 0.1
    })
    model({
      k <- exp(tk + eta.k)
      v <- exp(tv + eta.v)
      cp <- linCmt()
      cp ~ prop(prop.sd)
    })
  },
  lin2oral=function() {
    ini({
      tka <- log(1); tcl <- log(4); tv <- log(40); tq <- log(8); tvp <- log(80)
      eta.ka ~ 0.1; eta.cl ~ 0.1; eta.v ~ 0.1; eta.q ~ 0.1; eta.vp ~ 0.1
      add.sd <- 0.05; prop.sd <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      q <- exp(tq + eta.q)
      vp <- exp(tvp + eta.vp)
      cp <- linCmt()
      cp ~ add(add.sd) + prop(prop.sd)
    })
  },
  lin2iv=function() {
    ini({
      tcl <- log(4); tv <- log(40); tq <- log(8); tvp <- log(80)
      eta.cl ~ 0.1; eta.v ~ 0.1; eta.q ~ 0.1; eta.vp ~ 0.1
      prop.sd <- 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      q <- exp(tq + eta.q)
      vp <- exp(tvp + eta.vp)
      cp <- linCmt()
      cp ~ prop(prop.sd)
    })
  },
  lin2ivMicro=function() {
    ini({
      tk <- log(0.1); tk12 <- log(0.2); tk21 <- log(0.1); tv <- log(40)
      eta.k ~ 0.1; eta.k12 ~ 0.1; eta.k21 ~ 0.1; eta.v ~ 0.1
      prop.sd <- 0.1
    })
    model({
      k <- exp(tk + eta.k)
      k12 <- exp(tk12 + eta.k12)
      k21 <- exp(tk21 + eta.k21)
      v <- exp(tv + eta.v)
      cp <- linCmt()
      cp ~ prop(prop.sd)
    })
  },
  lin3oral=function() {
    ini({
      tka <- log(1); tcl <- log(4); tv <- log(40); tq <- log(8); tvp <- log(80)
      tq2 <- log(2); tvp2 <- log(200)
      eta.ka ~ 0.1; eta.cl ~ 0.1; eta.v ~ 0.1; eta.q ~ 0.1; eta.vp ~ 0.1
      eta.q2 ~ 0.1; eta.vp2 ~ 0.1
      prop.sd <- 0.1
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      q <- exp(tq + eta.q)
      vp <- exp(tvp + eta.vp)
      q2 <- exp(tq2 + eta.q2)
      vp2 <- exp(tvp2 + eta.vp2)
      cp <- linCmt()
      cp ~ prop(prop.sd)
    })
  },
  lin3iv=function() {
    ini({
      tcl <- log(4); tv <- log(40); tq <- log(8); tvp <- log(80)
      tq2 <- log(2); tvp2 <- log(200)
      eta.cl ~ 0.1; eta.v ~ 0.1; eta.q ~ 0.1; eta.vp ~ 0.1
      eta.q2 ~ 0.1; eta.vp2 ~ 0.1
      prop.sd <- 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      q <- exp(tq + eta.q)
      vp <- exp(tvp + eta.vp)
      q2 <- exp(tq2 + eta.q2)
      vp2 <- exp(tvp2 + eta.vp2)
      cp <- linCmt()
      cp ~ prop(prop.sd)
    })
  },
  lin1ivDur=function() {
    ini({
      tcl <- log(4); tv <- log(40); tdur <- log(2)
      eta.cl ~ 0.1; eta.v ~ 0.1; eta.dur ~ 0.1
      prop.sd <- 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d1 <- exp(tdur + eta.dur)
      dur(central) <- d1
      cp <- linCmt()
      cp ~ prop(prop.sd)
    })
  },
  lin1ivRate=function() {
    ini({
      tcl <- log(4); tv <- log(40); trate <- log(250)
      eta.cl ~ 0.1; eta.v ~ 0.1; eta.rate ~ 0.1
      prop.sd <- 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      r1 <- exp(trate + eta.rate)
      rate(central) <- r1
      cp <- linCmt()
      cp ~ prop(prop.sd)
    })
  },
  lin1oralWt=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; cl.wt <- 0.75
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + cl.wt * log(WT / 70))
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  },
  lin1oralEffect=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; tke0 <- log(0.5)
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; eta.ke0 ~ 0.1
      add.sd <- 0.7; add.e <- 0.3
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      ke0 <- exp(tke0 + eta.ke0)
      cp <- linCmt()
      d/dt(ce) <- ke0 * (cp - ce)
      cp ~ add(add.sd)
      ce ~ add(add.e)
    })
  },
  lin1oralOdeFirst=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; tke0 <- log(0.5)
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1; eta.ke0 ~ 0.1
      add.sd <- 0.7; add.e <- 0.3
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      ke0 <- exp(tke0 + eta.ke0)
      d/dt(ce) <- ke0 * (central / v - ce)
      cp <- linCmt()
      cp ~ add(add.sd)
      ce ~ add(add.e)
    })
  },
  lin1oralTime=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; cl.t <- 0.01
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl) * (1 + cl.t * t)
      v <- exp(tv + eta.v)
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  },
  lin1oralAmount=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      cp <- linCmt()
      cpa <- central / v
      cpa ~ add(add.sd)
    })
  },
  ode1=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
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
)

#' The ODE model with another residual error
#'
#' @param err residual error, like `"prop(prop.sd)"`
#' @param ini ini() lines for the residual error parameters
#' @return model function
#' @noRd
.stressResidual <- function(err, ini) {
  .txt <- paste0("function() {\n",
                 "  ini({\n",
                 "    tka <- 0.45; tcl <- 1; tv <- 3.45\n",
                 "    eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1\n",
                 "    ", ini, "\n",
                 "  })\n",
                 "  model({\n",
                 "    ka <- exp(tka + eta.ka)\n",
                 "    cl <- exp(tcl + eta.cl)\n",
                 "    v <- exp(tv + eta.v)\n",
                 "    d/dt(depot) <- -ka * depot\n",
                 "    d/dt(central) <- ka * depot - cl / v * central\n",
                 "    cp <- central / v\n",
                 "    cp ~ ", err, "\n",
                 "  })\n",
                 "}")
  eval(str2lang(.txt), envir=globalenv())
}

# models with an edge case in the model code (all based on ode1)
.stressCode <- list(
  ifElse=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; cl.sex <- 0.2
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      if (SEX == 1) {
        cl <- cl * (1 + cl.sex)
      }
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  },
  ifElseIf=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; cl.sex <- 0.2
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      if (SEX == 1) {
        cl <- cl * (1 + cl.sex)
      } else if (SEX == 2) {
        cl <- cl * (1 - cl.sex)
      } else {
        cl <- cl
      }
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  },
  nestedIf=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; cl.sex <- 0.2
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      if (SEX >= 1) {
        if (SEX == 1) {
          fcl <- 1 + cl.sex
        } else {
          fcl <- 1 - cl.sex
        }
      } else {
        fcl <- 1
      }
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - fcl * cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  },
  ifelseFun=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; cl.sex <- 0.2
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      fcl <- ifelse(SEX == 1, 1 + cl.sex, 1)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - fcl * cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  },
  reserved=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      s1 <- v
      ipred <- 1
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / s1 * central
      cp <- central / s1 * ipred
      cp ~ add(add.sd)
    })
  },
  longNames=function() {
    ini({
      t.absorption.rate <- 0.45; t.clearance.value <- 1; t.volume.central <- 3.45
      eta.absorption.rate ~ 0.6; eta.clearance.value ~ 0.3; eta.volume.central ~ 0.1
      add.sd.concentration <- 0.7
    })
    model({
      absorption.rate.constant <- exp(t.absorption.rate + eta.absorption.rate)
      clearance.of.drug <- exp(t.clearance.value + eta.clearance.value)
      volume.of.central <- exp(t.volume.central + eta.volume.central)
      d/dt(depot) <- -absorption.rate.constant * depot
      d/dt(central) <- absorption.rate.constant * depot - clearance.of.drug / volume.of.central * central
      concentration.in.plasma <- central / volume.of.central
      concentration.in.plasma ~ add(add.sd.concentration)
    })
  },
  fixedAndBlock=function() {
    ini({
      tka <- fix(0.45); tcl <- 1; tv <- 3.45
      eta.ka ~ fix(0.6)
      eta.cl + eta.v ~ c(0.3, 0.01, 0.1)
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
  },
  noEta=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  },
  initialCondition=function() {
    ini({
      tcl <- log(4); tv <- log(40); tbase <- log(1)
      eta.cl ~ 0.1; eta.v ~ 0.1
      prop.sd <- 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      base <- exp(tbase)
      central(0) <- base * v
      d/dt(central) <- -cl / v * central
      cp <- central / v
      cp ~ prop(prop.sd)
    })
  },
  timeVaryingCov=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; cl.crcl <- 0.5
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl) * (CRCL / 100)^cl.crcl
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  },
  probitInv=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45; tf <- 1
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      fd <- probitInv(tf)
      f(depot) <- fd
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  },
  fExpr = function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      lfdepot <- log(0.8)
      llag <- log(0.2)
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      f(depot) <- exp(lfdepot)
      alag(depot) <- exp(llag)
      d / dt(depot) <- -ka * depot
      d / dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  },
  iov=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
      iov.cl ~ 0.1 | occ
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + iov.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  },
  notMuRef=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl) + eta.cl
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  },
  tDist=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7; nu <- 3
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd) + dt(nu)
    })
  },
  fixedResidual=function() {
    ini({
      tka <- 0.45; tcl <- 1; tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- fix(0.7)
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
  },
  odeIv = function() {
    ini({
      tcl <- log(4)
      tv <- log(40)
      eta.cl ~ 0.1
      eta.v ~ 0.1
      prop.sd <- 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d / dt(central) <- -cl / v * central
      cp <- central / v
      cp ~ prop(prop.sd)
    })
  },
  odeIvDur = function() {
    ini({
      tcl <- log(4)
      tv <- log(40)
      tdur <- log(2)
      eta.cl ~ 0.1
      eta.v ~ 0.1
      eta.dur ~ 0.1
      prop.sd <- 0.1
    })
    model({
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      dur(central) <- exp(tdur + eta.dur)
      d / dt(central) <- -cl / v * central
      cp <- central / v
      cp ~ prop(prop.sd)
    })
  },
  odePkpd = function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      te0 <- log(100)
      timax <- fix(0.8)
      tic50 <- log(2)
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      eta.e0 ~ 0.05
      prop.sd <- 0.1
      add.eff <- 2
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      e0 <- exp(te0 + eta.e0)
      ic50 <- exp(tic50)
      d / dt(depot) <- -ka * depot
      d / dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      eff <- e0 * (1 - timax * cp / (ic50 + cp))
      cp ~ prop(prop.sd)
      eff ~ add(add.eff)
    })
  },
  bounded = function() {
    ini({
      tka <- c(-2, 0.45, 2)
      tcl <- c(0, 1, 3)
      tv <- c(1, 3.45, 6)
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- c(0, 0.7, 5)
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d / dt(depot) <- -ka * depot
      d / dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  },
  multiCov = function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      cl.wt <- 0.75
      cl.sex <- 0.2
      v.wt <- 1
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + cl.wt * log(WT / 70) + cl.sex * SEX2)
      v <- exp(tv + eta.v + v.wt * log(WT / 70))
      d / dt(depot) <- -ka * depot
      d / dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  },
  noEtaAtAll = function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }
)

# ---------------------------------------------------------------------
# cases
# ---------------------------------------------------------------------

#' A stress case
#'
#' @param name case name
#' @param model model function
#' @param data data set
#' @param nonmem,monolix "ok", "any" (translate or be refused with a
#'   known error) or a regular expression the refusal must match
#' @param checkNonmem,checkMonolix regular expressions that must be in
#'   the control stream/model files
#' @param run fit the case in run mode (otherwise translation only)
#' @param description description
#' @param controlNonmem,controlMonolix extra `nonmemControl()` or
#'   `monolixControl()` arguments (like `list(est="imp")`)
#' @param engines engines the case applies to
#' @param rerun in run mode, fit a second time and check that the saved
#'   output is read again (same objective function)
#' @return case list
#' @noRd
.stressCase <- function(
  name,
  model,
  data,
  nonmem = "ok",
  monolix = "ok",
  checkNonmem = NULL,
  checkMonolix = NULL,
  run = TRUE,
  description = "",
  controlNonmem = list(),
  controlMonolix = list(),
  engines = c("nonmem", "monolix"),
  rerun = FALSE
) {
  list(
    name = name,
    model = model,
    data = data,
    expect = list(nonmem = nonmem, monolix = monolix),
    check = list(nonmem = checkNonmem, monolix = checkMonolix),
    control = list(nonmem = controlNonmem, monolix = controlMonolix),
    engines = engines,
    rerun = rerun,
    run = run,
    description = description
  )
}

#' All the stress test cases
#'
#' @return list of cases (see `.stressCase()`)
#' @noRd
stressCases <- function() {
  if (!is.null(.stressEnv$cases)) return(.stressEnv$cases)
  .m <- .stressModels
  .theo <- .stressTheo()
  .theoPos <- .theo[.theo$EVID != 0 | .theo$DV > 0, ]
  .theoWt <- .theo
  .theoWt$WT <- 70 + 5 * (.theoWt$ID %% 5 - 2)
  .theoSex <- .theo
  .theoSex$SEX <- .theoSex$ID %% 3
  .theoCrcl <- .theo
  .theoCrcl$CRCL <- 90 + 2 * .theoCrcl$TIME + .theoCrcl$ID
  .theoOcc <- .theo
  .theoOcc$occ <- ifelse(.theoOcc$TIME > 8, 2, 1)
  # a second dose on day 2 is a second occasion
  .ivBolus <- .stressSim(.m$lin1iv, .stressEvIv())
  .ivInf <- .stressSim(.m$lin1iv, .stressEvIv(rate=250))
  .ivDurData <- .stressSim(.m$lin1iv, .stressEvIv(dur=2))
  .iv2 <- .stressSim(.m$lin2iv, .stressEvIv(times=c(.stressTimes, 48, 72)))
  .iv3 <- .stressSim(.m$lin3iv, .stressEvIv(times=c(.stressTimes, 48, 72, 96)))
  .oral2 <- .stressSim(.m$lin2oral, .stressEvOral(times=c(.stressTimes, 48, 72)))
  .oral3 <- .stressSim(.m$lin3oral, .stressEvOral(times=c(.stressTimes, 48, 72, 96)))
  .effect <- .stressSim(.m$lin1oralEffect, .stressEvOral(cmtObs=c("cp", "ce")))
  .effectFirst <- .effect
  # modeled duration/rate: the data say the model has the duration (-2)
  # or the rate (-1)
  .ivModelDur <- .ivBolus
  .ivModelDur$RATE <- ifelse(.ivModelDur$EVID == 1, -2, 0)
  .ivModelRate <- .ivBolus
  .ivModelRate$RATE <- ifelse(.ivModelRate$EVID == 1, -1, 0)
  # steady state and additional doses
  .ss <- .stressSim(.m$lin1oral,
                    rxode2::et(amt=320, cmt="depot", ii=24, ss=1) |>
                      rxode2::et(.stressTimes) |>
                      rxode2::et(id=seq_len(24)))
  .addl <- .stressSim(.m$lin1oral,
                      rxode2::et(amt=320, cmt="depot", ii=12, addl=3) |>
                        rxode2::et(c(.stressTimes, 36, 37, 38, 40, 44, 48)) |>
                        rxode2::et(id=seq_len(24)))
  # below the limit of quantification
  .cens <- .theoPos
  .loq <- 1
  .cens$CENS <- ifelse(.cens$EVID == 0 & .cens$DV < .loq, 1, 0)
  .cens$DV <- ifelse(.cens$CENS == 1, .loq, .cens$DV)
  .censLimit <- .cens
  .censLimit$LIMIT <- 0
  .limitOnly <- .theoPos
  .limitOnly$LIMIT <- 0
  # a reset and dose (evid=4) half way
  .reset <- .theo
  .r <- .reset[.reset$EVID != 0, ]
  .r$TIME <- 12.5
  .r$EVID <- 4
  .reset <- .reset[order(.reset$ID, .reset$TIME), ]
  .reset <- rbind(.reset, .r)
  .reset <- .reset[order(.reset$ID, .reset$TIME, -.reset$EVID), ]

  # data edge cases (all on theo_sd, fit with the ODE model)
  .theoNa <- .theo
  .w <- which(.theoNa$EVID == 0)
  .theoNa$DV[.w[seq(3, length(.w), by = 7)]] <- NA_real_
  .theoEvid2 <- .theo[.theo$EVID == 1, ]
  .theoEvid2$EVID <- 2L
  .theoEvid2$AMT <- 0
  .theoEvid2$DV <- NA_real_
  .theoEvid2$TIME <- 30
  .theoEvid2 <- rbind(.theo, .theoEvid2)
  .theoEvid2 <- .theoEvid2[
    order(.theoEvid2$ID, .theoEvid2$TIME, -.theoEvid2$EVID),
  ]
  .theoChrId <- .theo
  .theoChrId$ID <- paste0("S-", 100 + 7 * .theoChrId$ID)
  .theoShift <- .theo
  .theoShift$TIME <- .theoShift$TIME + 100
  .theoExtra <- .theo
  .theoExtra$NOTE <- ifelse(.theoExtra$EVID == 1, "dose", "sample")
  .theoExtra$STUDY <- 101
  .theoExtra$SITE <- .theoExtra$ID %% 4
  .theoMulti <- .theo
  .theoMulti$WT <- 70 + 5 * (.theoMulti$ID %% 5 - 2)
  .theoMulti$SEX2 <- .theoMulti$ID %% 2
  .theoOne <- .theo[.theo$ID == 1, ]
  # doses into the depot and into the central compartment
  .mixedRoutes <- .stressSim(
    .m$ode1,
    rxode2::et(amt = 320, cmt = "depot") |>
      rxode2::et(time = 24, amt = 100, cmt = "central") |>
      rxode2::et(c(.stressTimes, 24.25, 24.5, 25, 26, 28, 32, 36, 48)) |>
      rxode2::et(id = seq_len(24))
  )
  .odeInf <- .stressSim(.stressCode$odeIv, .stressEvIv(rate = 250))
  .odeInfSs <- .stressSim(
    .stressCode$odeIv,
    rxode2::et(amt = 500, cmt = "central", rate = 250, ii = 24, ss = 1) |>
      rxode2::et(.stressTimes) |>
      rxode2::et(id = seq_len(24))
  )
  .odeModelDur <- .stressSim(.stressCode$odeIv, .stressEvIv())
  .odeModelDur$RATE <- ifelse(.odeModelDur$EVID == 1, -2, 0)
  .pkpd <- .stressSim(
    .stressCode$odePkpd,
    .stressEvOral(cmtObs = c("cp", "eff"))
  )
  .res <- function(err, ini) .stressResidual(err, ini)

  .cases <- list(
    # linCmt() closed form ------------------------------------------
    .stressCase("linCmt 1-cmt oral", .m$lin1oral, .theo,
                checkNonmem=c("ADVAN2 TRANS1", "K=CL/V"),
                checkMonolix="pkmodel\\(V=rx_v, k=rx_k, ka=rx_ka\\)",
                description="linCmt() ~ endpoint; NONMEM ADVAN2, Monolix pkmodel()"),
    .stressCase("linCmt 1-cmt oral lag and bioavailability", .m$lin1oralLagF, .theo,
                checkNonmem=c("ADVAN2 TRANS1", "ALAG1=", "F1="),
                checkMonolix=c("Tlag=rx_tlag", "p=rx_p")),
    .stressCase("linCmt 1-cmt iv bolus", .m$lin1iv, .ivBolus,
                checkNonmem=c("ADVAN1 TRANS1", "RXE_CP=A\\(1\\)/"),
                checkMonolix="pkmodel\\(V=rx_v, k=rx_k\\)"),
    .stressCase("linCmt 1-cmt iv infusion (rate)", .m$lin1iv, .ivInf,
                checkNonmem=c("ADVAN1 TRANS1", "\\$INPUT.* RATE"),
                checkMonolix="pkmodel\\("),
    .stressCase("linCmt 1-cmt iv infusion (duration)", .m$lin1iv, .ivDurData,
                nonmem="duration/tinf",
                checkMonolix="pkmodel\\("),
    .stressCase("linCmt 1-cmt iv k/v parameterization", .m$lin1ivK, .ivBolus,
                checkNonmem="ADVAN1 TRANS1",
                checkMonolix="rx_k = k"),
    .stressCase("linCmt 2-cmt oral", .m$lin2oral, .oral2,
                checkNonmem=c("ADVAN4 TRANS1", "K23=", "K32=", "RXE_CP=A\\(2\\)/"),
                checkMonolix=c("k12=rx_k12", "ka=rx_ka")),
    .stressCase("linCmt 2-cmt iv", .m$lin2iv, .iv2,
                checkNonmem=c("ADVAN3 TRANS1", "K12=", "K21="),
                checkMonolix="k21=rx_k21"),
    .stressCase("linCmt 2-cmt iv micro-constants", .m$lin2ivMicro, .iv2,
                checkNonmem="ADVAN3 TRANS1",
                checkMonolix="pkmodel\\("),
    .stressCase("linCmt 3-cmt oral", .m$lin3oral, .oral3,
                checkNonmem=c("ADVAN12 TRANS1", "K24=", "K42="),
                checkMonolix=c("k13=rx_k13", "k31=rx_k31")),
    .stressCase("linCmt 3-cmt iv", .m$lin3iv, .iv3,
                checkNonmem=c("ADVAN11 TRANS1", "K13=", "K31="),
                checkMonolix="k13=rx_k13"),
    .stressCase("linCmt modeled duration", .m$lin1ivDur, .ivModelDur,
                checkNonmem=c("ADVAN1 TRANS1", "D1=", "\\$INPUT.* RATE"),
                checkMonolix = c("ddt_central", "Tk0=rx_dur_central"),
                description="Monolix pkmodel() cannot model the duration: ODEs"),
    .stressCase("linCmt modeled rate", .m$lin1ivRate, .ivModelRate,
                checkNonmem=c("ADVAN1 TRANS1", "R1=", "\\$INPUT.* RATE"),
                checkMonolix = c("ddt_central", "Tk0=amtDose/rx_rate_central"),
                description="Monolix pkmodel() cannot model the rate: ODEs"),
    .stressCase("linCmt weight covariate", .m$lin1oralWt, .theoWt,
                checkNonmem=c("ADVAN2 TRANS1", "\\$INPUT.* NLMIXRMUDERCOV1",
                              "\\+NLMIXRMUDERCOV1\\*THETA"),
                checkMonolix="pkmodel\\("),
    .stressCase("linCmt with effect compartment ODE", .m$lin1oralEffect, .effect,
                checkNonmem=c("ADVAN13", "DADT\\(3\\)"),
                checkMonolix=c("ddt_ce", "ddt_central"),
                description="mixed linCmt()/ODE models are solved as ODEs"),
    .stressCase("linCmt with ODE defined first", .m$lin1oralOdeFirst, .effectFirst,
                checkNonmem=c("ADVAN13", "DADT\\(3\\) = KE0", "DADT\\(1\\) = - \\(KA\\)\\*A\\(1\\)"),
                checkMonolix="ddt_ce",
                description="compartment numbers stay depot=1, central=2, ce=3"),
    # nlmixr2est's mu2 covariate hook reads `t` as a data column (and
    # finds base::t()); when that is fixed upstream this case reports
    # "not refused" and the expectation should become "ok"
    .stressCase("linCmt time-dependent clearance", .m$lin1oralTime, .theo,
                nonmem="replicate an object of type 'closure'",
                monolix="replicate an object of type 'closure'",
                checkNonmem="ADVAN13",
                checkMonolix="ddt_central",
                description="parameters changing with time need ODEs (known nlmixr2est mu2 issue with t)"),
    .stressCase("linCmt amount in the endpoint", .m$lin1oralAmount, .theo,
                checkNonmem=c("ADVAN2 TRANS1", "A\\(2\\)"),
                checkMonolix="ddt_central",
                description="Monolix pkmodel() has no amounts: ODEs"),
    .stressCase("linCmt steady state", .m$lin1oral, .ss,
                checkNonmem=c("ADVAN2 TRANS1", "\\$INPUT.* II .* SS"),
                checkMonolix="pkmodel\\("),
    .stressCase("linCmt steady state with lag time", .m$lin1oralLagF, .ss,
                checkNonmem=c("ADVAN2 TRANS1", "\\$INPUT.* II .* SS", "ALAG1="),
                checkMonolix="Tlag=rx_tlag",
                description="rxode2 splits a lagged steady state dose in two; NONMEM/Monolix get one SS dose"),
    .stressCase("linCmt additional doses", .m$lin1oral, .addl,
                checkNonmem="ADVAN2 TRANS1",
                checkMonolix="pkmodel\\("),
    # residual errors --------------------------------------------------
    .stressCase("residual prop", .res("prop(prop.sd)", "prop.sd <- 0.1"), .theoPos),
    .stressCase("residual combined1",
                .res("add(add.sd) + prop(prop.sd) + combined1()",
                     "add.sd <- 0.2; prop.sd <- 0.1"),
                .theoPos),
    .stressCase("residual combined2",
                .res("add(add.sd) + prop(prop.sd) + combined2()",
                     "add.sd <- 0.2; prop.sd <- 0.1"),
                .theoPos),
    .stressCase("residual pow", .res("pow(pow.sd, pw)", "pow.sd <- 0.1; pw <- 1"),
                .theoPos, monolix="pow"),
    .stressCase("residual lnorm", .res("lnorm(lnorm.sd)", "lnorm.sd <- 0.1"),
                .theoPos, checkMonolix="logNormal|lognormal|exponential"),
    .stressCase("residual logit",
                .res("logitNorm(logit.sd, 0, 20)", "logit.sd <- 0.1"),
                .theoPos, checkMonolix="logitNormal, min=0, max=20"),
    .stressCase("residual boxCox",
                .res("add(add.sd) + boxCox(lambda)", "add.sd <- 0.7; lambda <- 0.5"),
                .theoPos, monolix="boxCox|yeoJohnson|not supported|transform"),
    .stressCase("residual yeoJohnson",
                .res("add(add.sd) + yeoJohnson(lambda)", "add.sd <- 0.7; lambda <- 0.5"),
                .theoPos, monolix="boxCox|yeoJohnson|not supported|transform"),
    .stressCase("residual t distribution", .stressCode$tDist, .theoPos,
                nonmem="t|not supported|normal", monolix="t|not supported|normal"),
    .stressCase("fixed residual error", .stressCode$fixedResidual, .theoPos,
                checkNonmem="FIX"),
    # censoring --------------------------------------------------------
    .stressCase("censoring (CENS)", .m$lin1oral, .cens,
                checkNonmem="CENS"),
    .stressCase("censoring (CENS and LIMIT)", .m$lin1oral, .censLimit,
                checkNonmem="LIMIT"),
    .stressCase("censoring LIMIT only", .m$lin1oral, .limitOnly,
                checkNonmem = "LIMIT .GT.",
                description = "M2: a finite LIMIT without CENS"),
    # model code ---------------------------------------------------------
    .stressCase("if/else", .stressCode$ifElse, .theoSex, checkNonmem="IF \\("),
    # NONMEM (and Monolix for ifelse()) prune the if/else branches (#11)
    .stressCase("else if", .stressCode$ifElseIf, .theoSex, checkNonmem="RXL"),
    .stressCase("nested if/else", .stressCode$nestedIf, .theoSex,
                checkNonmem="RXL"),
    .stressCase("ifelse()", .stressCode$ifelseFun, .theoSex,
                checkNonmem="RXL", checkMonolix="rx_l"),
    .stressCase("NONMEM reserved names", .stressCode$reserved, .theo,
                checkNonmem="RXR1"),
    .stressCase("long dotted names", .stressCode$longNames, .theo),
    .stressCase("fixed and block random effects", .stressCode$fixedAndBlock, .theo,
                checkNonmem=c("BLOCK\\(2\\)", "FIX")),
    .stressCase("parameters without random effects", .stressCode$noEta, .theo),
    .stressCase("initial condition", .stressCode$initialCondition, .ivBolus,
                checkNonmem="A_0\\(1\\)"),
    .stressCase("time-varying covariate", .stressCode$timeVaryingCov, .theoCrcl,
                checkMonolix="CRCL"),
    .stressCase("probitInv", .stressCode$probitInv, .theo),
    .stressCase(
      "f()/alag() expressions", .stressCode$fExpr, .theo,
      checkNonmem = c("F1=", "ALAG1="),
      checkMonolix = c("Tlag=rx_lag_depot, p=rx_f_depot", "rx_f_depot = exp"),
      description = "issue #115"
    ),
    .stressCase("between-occasion variability", .stressCode$iov, .theoOcc,
                nonmem="id|occasion|level|random", monolix="id|occasion|level|random",
                run=FALSE),
    .stressCase("not mu-referenced", .stressCode$notMuRef, .theo,
                monolix="mu"),
    # data -------------------------------------------------------------
    .stressCase("reset and dose (evid=4)", .m$ode1, .reset),
    .stressCase("missing observations (DV=NA)", .m$ode1, .theoNa),
    .stressCase("other events (evid=2)", .m$ode1, .theoEvid2),
    .stressCase(
      "character IDs",
      .m$ode1,
      .theoChrId,
      description = "IDs like S-107 are renumbered for NONMEM/Monolix"
    ),
    .stressCase("time not starting at zero", .m$ode1, .theoShift),
    .stressCase(
      "extra unused data columns",
      .m$ode1,
      .theoExtra,
      description = "character, constant and unused numeric columns"
    ),
    .stressCase(
      "single subject",
      .stressCode$noEta,
      .theoOne,
      nonmem = "any",
      monolix = "any",
      description = "one subject with one random effect"
    ),
    .stressCase(
      "oral and iv doses",
      .m$ode1,
      .mixedRoutes,
      description = "doses into the depot and into the central compartment"
    ),
    .stressCase(
      "ODE infusion (rate)",
      .stressCode$odeIv,
      .odeInf,
      checkNonmem = "\\$INPUT.* RATE"
    ),
    .stressCase(
      "ODE steady state infusion",
      .stressCode$odeIv,
      .odeInfSs,
      checkNonmem = "\\$INPUT.* SS"
    ),
    .stressCase(
      "ODE modeled duration",
      .stressCode$odeIvDur,
      .odeModelDur,
      checkNonmem = c("D1=", "\\$INPUT.* RATE"),
      checkMonolix = "Tk0="
    ),
    .stressCase(
      "ODE PK/PD two endpoints",
      .stressCode$odePkpd,
      .pkpd,
      description = "prop error on cp, add error on eff (CMT is the endpoint)"
    ),
    .stressCase(
      "bounded thetas",
      .stressCode$bounded,
      .theo,
      checkNonmem = "\\(-2,"
    ),
    .stressCase("several covariates", .stressCode$multiCov, .theoMulti),
    # NONMEM gets a control stream without $OMEGA; run mode shows
    # whether NONMEM fits it
    .stressCase(
      "no random effects",
      .stressCode$noEtaAtAll,
      .theo,
      nonmem = "any",
      monolix = "mixed effect",
      description = "fixed effects only; Monolix needs random effects"
    ),
    # estimation options -----------------------------------------------
    .stressCase(
      "NONMEM est=imp",
      .m$ode1,
      .theo,
      engines = "nonmem",
      controlNonmem = list(est = "imp"),
      checkNonmem = "METHOD=IMP"
    ),
    .stressCase(
      "NONMEM est=its",
      .m$ode1,
      .theo,
      engines = "nonmem",
      controlNonmem = list(est = "its"),
      checkNonmem = "METHOD=ITS"
    ),
    .stressCase(
      "NONMEM est=posthoc",
      .m$ode1,
      .theo,
      engines = "nonmem",
      controlNonmem = list(est = "posthoc"),
      checkNonmem = "MAXEVALS=0"
    ),
    .stressCase(
      "NONMEM no covariance step",
      .m$ode1,
      .theo,
      engines = "nonmem",
      controlNonmem = list(cov = "")
    ),
    .stressCase(
      "NONMEM ADVAN6",
      .m$ode1,
      .theo,
      engines = "nonmem",
      controlNonmem = list(advanOde = "advan6"),
      checkNonmem = "ADVAN6"
    ),
    .stressCase(
      "NONMEM linCmt as ODEs",
      .m$lin2oral,
      .oral2,
      engines = "nonmem",
      controlNonmem = list(linCmt = "ode"),
      checkNonmem = "ADVAN13"
    ),
    .stressCase(
      "Monolix linearization",
      .m$ode1,
      .theo,
      engines = "monolix",
      controlMonolix = list(useLinearization = TRUE),
      checkMonolix = "method = Linearization"
    ),
    .stressCase(
      "Monolix linCmt as ODEs",
      .m$lin2oral,
      .oral2,
      engines = "monolix",
      controlMonolix = list(linCmt = "ode"),
      checkMonolix = "ddt_central"
    ),
    .stressCase(
      "Monolix stiff ODEs",
      .m$ode1,
      .theo,
      engines = "monolix",
      controlMonolix = list(stiff = TRUE),
      checkMonolix = "odeType = stiff"
    ),
    .stressCase(
      "Monolix decreasing variability",
      .m$ode1,
      .theo,
      engines = "monolix",
      controlMonolix = list(variability = "decreasing"),
      checkMonolix = "variability = decreasing"
    ),
    # a second fit reads the saved output instead of running again
    .stressCase("rerun reads saved output", .m$lin1oral, .theo, rerun = TRUE)
  )
  names(.cases) <- vapply(.cases, function(x) x$name, character(1))
  .stressEnv$cases <- .cases
  .cases
}

# ---------------------------------------------------------------------
# nlmixr2lib sweep
# ---------------------------------------------------------------------

#' Models from nlmixr2lib to stress the translation
#'
#' @param all use every model (otherwise every linCmt() model and a
#'   sample of the others)
#' @param n number of models per category in the sample
#' @return data frame of the nlmixr2lib models to use
#' @noRd
stressLibModels <- function(all=FALSE, n=3L) {
  .db <- nlmixr2lib::modeldb
  if (all) return(.db)
  .lin <- .db[.db$linCmt, ]
  .other <- .db[!.db$linCmt, ]
  .cat <- ifelse(is.na(.other$category), "", .other$category)
  .other <- do.call(rbind, lapply(split(.other, .cat), function(d) {
    d[seq_len(min(nrow(d), n)), ]
  }))
  .ret <- rbind(.lin, .other)
  rownames(.ret) <- NULL
  .ret
}

#' Simulated data for an nlmixr2lib model
#'
#' Doses go into the model's dosing compartments; every endpoint is
#' observed; covariates are set to typical values.
#'
#' @param ui rxode2 ui of the model
#' @param dosing dosing compartments (from nlmixr2lib's modeldb)
#' @return data frame for nlmixr2
#' @noRd
stressLibData <- function(ui, dosing) {
  .state <- rxode2::rxState(ui)
  .dosing <- trimws(strsplit(paste(dosing), ",")[[1]])
  .dosing <- .dosing[.dosing %in% .state]
  if (length(.dosing) == 0L) .dosing <- .state[1]
  # models without compartments have no doses
  .ev <- if (length(.state) == 0L) rxode2::et() else rxode2::et(amt=100, cmt=.dosing[1])
  .ends <- ui$predDf
  if (is.null(.ends) || nrow(.ends) == 0L) stop("no endpoints", call.=FALSE)
  .multi <- nrow(.ends) > 1L
  for (.i in seq_len(nrow(.ends))) {
    if (.multi) {
      .ev <- rxode2::et(.ev, .stressTimes, cmt=paste(.ends$cond[.i]))
    } else if (.i == 1L) {
      .ev <- rxode2::et(.ev, .stressTimes)
    }
  }
  .ev <- as.data.frame(rxode2::et(.ev, id=1:6))
  # everything the model uses that is not estimated comes from the data
  .mv <- rxode2::rxModelVars(ui)
  .covs <- setdiff(.mv$params, c(ui$iniDf$name, "t", "time", "tad", "tafd",
                                 "tlast", "tfirst", "podo", "dosenum"))
  for (.c in unique(c(ui$allCovs, .covs))) {
    .ev[[.c]] <- .stressCovValue(.c)
  }
  .ev$dv <- ifelse(.ev$evid == 0, 1, NA_real_)
  # the event items are upper case; covariates keep the model's names
  .items <- c("id", "time", "evid", "amt", "cmt", "dv", "ii", "ss", "rate", "dur", "addl")
  .w <- names(.ev) %in% .items
  names(.ev)[.w] <- toupper(names(.ev)[.w])
  .ev
}

.stressCovValue <- function(name) {
  .n <- toupper(name)
  if (grepl("^(WT|BW|WEIGHT|BWT|TBW|FFM|LBM)", .n)) return(70)
  if (grepl("^(HT|HEIGHT)", .n)) return(170)
  if (grepl("AGE", .n)) return(40)
  if (grepl("^(CRCL|EGFR|GFR|CLCR)", .n)) return(100)
  if (grepl("^ALB", .n)) return(4)
  if (grepl("BMI", .n)) return(24)
  if (grepl("BSA", .n)) return(1.8)
  1
}

#' Stress cases from nlmixr2lib
#'
#' @inheritParams stressLibModels
#' @return list of cases; the expected outcome is "any", that is a
#'   successful translation or a documented refusal
#' @noRd
stressLibCases <- function(all=FALSE, n=3L) {
  .db <- stressLibModels(all=all, n=n)
  .ret <- lapply(seq_len(nrow(.db)), function(i) {
    .name <- .db$name[i]
    .model <- try(suppressMessages(rxode2::rxode2(nlmixr2lib::readModelDb(.name))), silent=TRUE)
    .data <- try(stressLibData(.model, .db$dosing[i]), silent=TRUE)
    .stressCase(paste0("nlmixr2lib::", .name), .model, .data,
                nonmem="any", monolix="any", run=FALSE,
                description=.db$description[i])
  })
  names(.ret) <- vapply(.ret, function(x) x$name, character(1))
  .ret
}

# errors that are documented refusals (the translation is not supported);
# anything else from an nlmixr2lib model is reported as an error
.stressKnownRefusals <- c(
  # up-front assertions
  "needs to be a completely mu-referenced model",
  "needs to be a mixed effect model",
  "needs to be a \\(transformably\\) normal model",
  "can only have random effects on ID",
  "residual parameters cannot depend on the model calculated parameters",
  "cannot use the residual (error|transformation)",
  "cannot fix residual error parameters",
  # translation limits
  "is not supported by babelmixr2",
  "is not supported in (NONMEM|monolix)<->nlmixr",
  "will not allow `else if` or `else` statements",
  "will not handle nested if/else",
  "unknown rxode2 assignment type",
  "does not support the parameter transformation",
  "cannot be translated to (NONMEM|monolix)",
  "strings in nlmixr<->",
  "ylo/yup not implemented",
  "probit not supported in nonmem",
  "residual type '[^']+' is not supported",
  # data limits
  "does not support a duration/tinf data item",
  "events are not supported in",
  "does not support a fixed duration",
  "events are not supported",
  "steady state infusions are not supported",
  "complex steady state \\(ss=2\\) are not supported",
  "all transformations need to be fixed",
  # priors
  "cannot translate the prior on", "PRIOR NWPRI gives priors",
  "cannot be written for Monolix", "the model puts a prior on",
  "the prior on the residual error parameter",
  # linCmt()
  "could not translate this linCmt\\(\\) model"
)

# errors from other nlmixr2 packages (not babelmixr2); they are reported
# with the status "upstream" and do not fail the tests
.stressKnownUpstream <- c(
  # nlmixr2est's mu2 covariate hook reads `t` as a data column
  "replicate an object of type 'closure'",
  "were in the ini block but not in the model block",
  "rxode2 syntax error"
)

# ---------------------------------------------------------------------
# checks of the files that are written
# ---------------------------------------------------------------------

.stressBalanced <- function(line) {
  .open <- nchar(gsub("[^(]", "", line))
  .close <- nchar(gsub("[^)]", "", line))
  .open == .close
}

#' Problems in a NONMEM control stream
#'
#' @param lines control stream lines
#' @return character vector of problems (empty when none)
#' @noRd
stressLintNonmem <- function(lines) {
  .ret <- character(0)
  # only the abbreviated code records are checked for syntax
  .rec <- cumsum(grepl("^\\$", lines))
  .recName <- sub("^\\$([A-Z]+).*$", "\\1", lines[grepl("^\\$", lines)])
  .isCode <- .rec > 0 & .recName[pmax(.rec, 1)] %in% c("PK", "PRED", "DES", "ERROR")
  .code <- sub(";.*$", "", lines)
  .code[!.isCode] <- ""
  .bad <- which(!vapply(.code, .stressBalanced, logical(1)))
  if (length(.bad) > 0L) .ret <- c(.ret, paste0("unbalanced parentheses: ", lines[.bad]))
  .bad <- grep("<-|~|linCmt|\\bNA\\b|\\bNaN\\b|\\bInf\\b|rxLinCmt[^O]", .code)
  if (length(.bad) > 0L) .ret <- c(.ret, paste0("R syntax left in code: ", lines[.bad]))
  .bad <- grep("=\\s*$", .code)
  if (length(.bad) > 0L) .ret <- c(.ret, paste0("assignment without value: ", lines[.bad]))
  # NONMEM cannot use a logical expression as a number (only in IF ())
  .if <- paste0("^\\s*(ELSE )?IF\\s*\\(.*\\)( THEN)?",
                "(\\s+[A-Z_0-9]+\\s*=\\s*[0-9.]+([ED][-+]?[0-9]+)?)?\\s*$")
  .bad <- grep("\\.(EQ|NE|GT|GE|LT|LE|AND|OR|NOT)\\.", sub(.if, "", .code))
  if (length(.bad) > 0L) {
    .ret <- c(.ret, paste0("logical expression used as a number: ", lines[.bad]))
  }
  if (!any(grepl("^\\$(PK|PRED)", lines))) .ret <- c(.ret, "no $PK/$PRED record")
  if (any(grepl("^\\$SUBROUTINES ADVAN(1|2|3|4|11|12) ", lines)) &&
        any(grepl("^\\$DES", lines))) {
    .ret <- c(.ret, "closed-form ADVAN with $DES")
  }
  .ret
}

#' Problems in a Monolix model file
#'
#' @param lines model file lines
#' @return character vector of problems (empty when none)
#' @noRd
stressLintMonolix <- function(lines) {
  .ret <- character(0)
  .code <- sub(";.*$", "", lines)
  .bad <- which(!vapply(.code, .stressBalanced, logical(1)))
  if (length(.bad) > 0L) .ret <- c(.ret, paste0("unbalanced parentheses: ", lines[.bad]))
  .bad <- grep("<-|~|linCmt\\(|\\bNA\\b", .code)
  if (length(.bad) > 0L) .ret <- c(.ret, paste0("R syntax left in code: ", lines[.bad]))
  # macro arguments run together, like empty(adm=1target=...)
  .bad <- grep("=[^,(){}= ]+[A-Za-z]+=", .code)
  if (length(.bad) > 0L) .ret <- c(.ret, paste0("macro arguments without a comma: ", lines[.bad]))
  .bad <- grep("=\\s*$", .code)
  if (length(.bad) > 0L) .ret <- c(.ret, paste0("assignment without value: ", lines[.bad]))
  if (any(grepl("pkmodel\\(", .code)) && any(grepl("^\\s*ddt_", .code))) {
    .ret <- c(.ret, "pkmodel() with ODEs")
  }
  .ret
}

.stressLintMlxtran <- function(lines) {
  .ret <- character(0)
  .bad <- grep("min=[^,]+max=", lines)
  if (length(.bad) > 0L) .ret <- c(.ret, paste0("distribution arguments without a comma: ", lines[.bad]))
  .ret
}

#' Problems in the data written for NONMEM/Monolix
#'
#' @param data data frame read from the written csv
#' @return character vector of problems (empty when none)
#' @noRd
.stressLintData <- function(data) {
  .ret <- character(0)
  .n <- toupper(names(data))
  names(data) <- .n
  if (any(is.na(data$ID)) || any(is.na(data$TIME))) {
    .ret <- c(.ret, "missing ID or TIME in the data")
  }
  if (all(c("EVID", "AMT") %in% .n)) {
    .dose <- data[!is.na(data$EVID) & data$EVID %in% c(1, 4), ]
    .key <- intersect(c("ID", "TIME", "EVID", "AMT", "CMT", "ADM"), .n)
    if (nrow(.dose) > 0L && anyDuplicated(.dose[, .key]) > 0L) {
      .ret <- c(.ret, "duplicated dose records in the data")
    }
  }
  .ret
}

# ---------------------------------------------------------------------
# running a case
# ---------------------------------------------------------------------

.stressControl <- function(
  engine,
  mode,
  modelName,
  runCommand,
  extra = list()
) {
  .cmd <- if (mode == "translate") NA else runCommand
  if (engine == "nonmem") {
    if (is.null(.cmd)) {
      .cmd <- getOption("babelmixr2.nonmem", "")
    }
    do.call(
      babelmixr2::nonmemControl,
      c(list(runCommand = .cmd, modelName = modelName), extra)
    )
  } else {
    if (is.null(.cmd)) {
      .cmd <- getOption("babelmixr2.monolix", "")
    }
    do.call(
      babelmixr2::monolixControl,
      c(list(runCommand = .cmd, modelName = modelName), extra)
    )
  }
}

.stressFit <- function(case, engine, mode, dir, modelName, runCommand) {
  withr::with_dir(dir, {
    tryCatch(
      suppressWarnings(suppressMessages(
        nlmixr2est::nlmixr2(
          case$model,
          case$data,
          est = engine,
          control = .stressControl(
            engine,
            mode,
            modelName,
            runCommand,
            case$control[[engine]]
          )
        )
      )),
      error = function(e) e
    )
  })
}

#' How well rxode2 reproduces the NONMEM/Monolix predictions
#'
#' @param fit babelmixr2 fit
#' @param engine "nonmem" or "monolix"
#' @return median relative differences (%) of IPRED and PRED
#' @noRd
.stressPredDiff <- function(fit, engine) {
  .f <- if (engine == "nonmem") {
    ".nonmemMergePredsAndCalcRelativeErr"
  } else {
    ".monolixMergePredsAndCalcRelativeErr"
  }
  .d <- try(
    suppressWarnings(
      utils::getFromNamespace(.f, "babelmixr2")(fit)
    ),
    silent = TRUE
  )
  if (inherits(.d, "try-error") || is.null(.d)) {
    return(c(NA_real_, NA_real_))
  }
  # the quantiles are 0, lower ci, median, upper ci, 1
  c(unname(.d$individualRel[3]), unname(.d$popRel[3]))
}

.stressFiles <- function(engine, modelName) {
  if (engine == "nonmem") {
    list(model=file.path(paste0(modelName, "-nonmem"), paste0(modelName, ".nmctl")),
         data=file.path(paste0(modelName, "-nonmem"), paste0(modelName, ".csv")))
  } else {
    list(model=paste0(modelName, "-monolix.txt"),
         mlxtran=paste0(modelName, "-monolix.mlxtran"),
         data=paste0(modelName, "-monolix.csv"))
  }
}

.stressModelName <- function(name) {
  .n <- gsub("[^A-Za-z0-9]+", "_", name)
  .n <- gsub("^_|_$", "", .n)
  substr(.n, 1, 60)
}

#' Run one stress case for one engine
#'
#' @param case case (see `.stressCase()`)
#' @param engine "nonmem" or "monolix"
#' @param mode "translate" (write the files only) or "run" (fit with
#'   the external program)
#' @param dir directory to write the files into
#' @param runCommand run command for the engine (`NULL` uses the
#'   babelmixr2 options)
#' @param reference fit the model with nlmixr2 as well and compare
#'   the population estimates (run mode only)
#' @param predTol in run mode, the largest median relative difference
#'   (in %) allowed between the rxode2 and NONMEM/Monolix IPRED
#' @return one row data frame with the result
#' @noRd
stressRunCase <- function(
  case,
  engine,
  mode = "translate",
  dir = tempfile("stress"),
  runCommand = NULL,
  reference = FALSE,
  predTol = 5
) {
  .expect <- case$expect[[engine]]
  .modelName <- .stressModelName(case$name)
  .ret <- data.frame(
    case = case$name,
    engine = engine,
    mode = mode,
    expect = .expect,
    status = NA_character_,
    message = "",
    problems = "",
    seconds = NA_real_,
    objf = NA_real_,
    ipredRelDiff = NA_real_,
    predRelDiff = NA_real_,
    rerunSeconds = NA_real_,
    maxRelDiffTheta = NA_real_,
    stringsAsFactors = FALSE
  )
  .setup <- Filter(function(x) inherits(x, "try-error"), list(case$model, case$data))
  if (length(.setup) > 0L) {
    .ret$status <- "setup"
    .ret$message <- paste(vapply(.setup, paste, character(1)), collapse = " ")
    return(.ret)
  }
  if (mode == "run" && !isTRUE(case$run)) {
    .ret$status <- "skipped"
    .ret$message <- "case is translation only"
    return(.ret)
  }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  .time <- proc.time()
  .fit <- .stressFit(case, engine, mode, dir, .modelName, runCommand)
  .ret$seconds <- (proc.time() - .time)[["elapsed"]]
  if (inherits(.fit, "error")) {
    .ret$message <- gsub("\n", " ", conditionMessage(.fit))
    .ret$status <- .stressErrorStatus(.expect, conditionMessage(.fit))
    return(.ret)
  }
  if (!(.expect %in% c("ok", "any"))) {
    .ret$status <- "not refused"
    .ret$message <- paste0("expected an error matching '", .expect, "'")
    return(.ret)
  }
  .problems <- .stressCheckFiles(case, engine, dir, .modelName)
  if (mode == "run") {
    .run <- .stressCheckFit(
      .fit, case, engine, dir, .modelName, runCommand, reference, predTol
    )
    .ret[names(.run$values)] <- .run$values
    .problems <- c(.problems, .run$problems)
  }
  .ret$problems <- paste(.problems, collapse = " | ")
  .ret$status <- if (length(.problems) > 0L) "problem" else "ok"
  .ret
}

#' Status of a case that gave an error
#'
#' @param expect what the case expects ("ok", "any" or a regular
#'   expression the refusal must match)
#' @param msg error message
#' @return "refused", "upstream" or "error"
#' @noRd
.stressErrorStatus <- function(expect, msg) {
  if (identical(expect, "ok")) {
    return("error")
  }
  if (!identical(expect, "any")) {
    return(if (grepl(expect, msg)) "refused" else "error")
  }
  if (any(vapply(.stressKnownRefusals, grepl, logical(1), x = msg))) {
    return("refused")
  }
  .upstream <- vapply(.stressKnownUpstream, grepl, logical(1), x = msg, fixed = TRUE)
  if (any(.upstream)) "upstream" else "error"
}

#' Check the files written for NONMEM/Monolix
#'
#' @inheritParams stressRunCase
#' @param modelName model name (file names)
#' @return character vector of problems
#' @noRd
.stressCheckFiles <- function(case, engine, dir, modelName) {
  .files <- .stressFiles(engine, modelName)
  .problems <- character(0)
  .modelFile <- file.path(dir, .files$model)
  if (!file.exists(.modelFile)) {
    .problems <- paste0("missing ", .files$model)
  } else {
    .lines <- readLines(.modelFile)
    .problems <- if (engine == "nonmem") {
      stressLintNonmem(.lines)
    } else {
      stressLintMonolix(.lines)
    }
    .mlxtran <- file.path(dir, .files$mlxtran)
    if (engine == "monolix" && file.exists(.mlxtran)) {
      .lines <- c(.lines, readLines(.mlxtran))
      .problems <- c(.problems, .stressLintMlxtran(readLines(.mlxtran)))
    }
    .found <- vapply(case$check[[engine]], function(re) any(grepl(re, .lines)), logical(1))
    .problems <- c(.problems, sprintf("missing '%s'", case$check[[engine]][!.found]))
  }
  .dataFile <- file.path(dir, .files$data)
  if (!file.exists(.dataFile)) {
    return(c(.problems, paste0("missing ", .files$data)))
  }
  c(.problems, .stressLintData(utils::read.csv(.dataFile, na.strings = ".")))
}

#' Check a NONMEM/Monolix fit (run mode)
#'
#' @param fit the fit
#' @inheritParams stressRunCase
#' @param modelName model name (file names)
#' @return list with `values` (result columns) and `problems`
#' @noRd
.stressCheckFit <- function(
  fit,
  case,
  engine,
  dir,
  modelName,
  runCommand,
  reference,
  predTol
) {
  if (!inherits(fit, "nlmixr2FitData")) {
    return(list(values = list(), problems = "the fit did not return an nlmixr2 fit"))
  }
  .problems <- character(0)
  if (!is.finite(fit$objf)) {
    .problems <- "objective function is not finite"
  }
  .bad <- names(fit$theta)[!is.finite(fit$theta)]
  if (length(.bad) > 0L) {
    .problems <- c(.problems, paste0("non-finite estimates: ", paste(.bad, collapse = ", ")))
  }
  .pd <- .stressPredDiff(fit, engine)
  if (is.na(.pd[1])) {
    .problems <- c(.problems, "cannot compare the predictions with rxode2")
  } else if (.pd[1] > predTol) {
    .problems <- c(
      .problems,
      sprintf(
        "rxode2 IPRED differs from %s by %.2f%% (median; limit %g%%)",
        engine,
        .pd[1],
        predTol
      )
    )
  }
  .values <- list(objf = fit$objf, ipredRelDiff = .pd[1], predRelDiff = .pd[2])
  if (isTRUE(case$rerun)) {
    .rerun <- .stressCheckRerun(fit, case, engine, dir, modelName, runCommand)
    .values$rerunSeconds <- .rerun$seconds
    .problems <- c(.problems, .rerun$problems)
  }
  if (reference) {
    .values$maxRelDiffTheta <- .stressCompare(fit, case, engine)
  }
  list(values = .values, problems = .problems)
}

#' Fit a case again: the saved output should be read again
#'
#' @inheritParams .stressCheckFit
#' @return list with `seconds` and `problems`
#' @noRd
.stressCheckRerun <- function(fit, case, engine, dir, modelName, runCommand) {
  .time <- proc.time()
  .fit2 <- .stressFit(case, engine, "run", dir, modelName, runCommand)
  .seconds <- (proc.time() - .time)[["elapsed"]]
  .problems <- if (inherits(.fit2, "error")) {
    paste0("rerun failed: ", conditionMessage(.fit2))
  } else if (!isTRUE(all.equal(fit$objf, .fit2$objf))) {
    "rerun gave a different objective function"
  }
  list(seconds = .seconds, problems = .problems)
}

#' Compare an external fit to the same model fit in nlmixr2
#'
#' @param fit babelmixr2 fit
#' @param case stress case
#' @param engine engine used for the fit
#' @return maximum relative difference of the population parameters
#' @noRd
.stressCompare <- function(fit, case, engine) {
  .est <- if (engine == "nonmem") "focei" else "saem"
  .ref <- try(suppressWarnings(suppressMessages(
    nlmixr2est::nlmixr2(case$model, case$data, est=.est,
                        control=list(print=0)))), silent=TRUE)
  if (inherits(.ref, "try-error")) return(NA_real_)
  .t1 <- fit$theta
  .t2 <- .ref$theta[names(.t1)]
  max(abs(.t1 - .t2) / pmax(abs(.t2), 1e-3), na.rm=TRUE)
}

#' Run the stress cases
#'
#' @param cases list of cases
#' @param engines engines to use
#' @inheritParams stressRunCase
#' @param progress print progress
#' @return data frame of results (one row per case and engine)
#' @noRd
stressRun <- function(cases=stressCases(), engines=c("nonmem", "monolix"),
                      mode="translate", dir=tempfile("stress"),
                      runCommand=list(nonmem=NULL, monolix=NULL),
                      reference = FALSE, predTol = 5,
                      progress = interactive()) {
  .ret <- list()
  for (.case in cases) {
    for (.engine in intersect(engines, .case$engines)) {
      .dir <- file.path(dir, .engine, .stressModelName(.case$name))
      .r <- stressRunCase(.case, .engine, mode=mode, dir=.dir,
                          runCommand = runCommand[[.engine]],
                          reference = reference, predTol = predTol)
      if (progress) {
        message(sprintf("%-8s %-12s %s %s", .engine, .r$status, .case$name,
                        ifelse(.r$status %in% c("ok", "skipped"), "",
                               paste(.r$message, .r$problems))))
      }
      .ret[[length(.ret) + 1L]] <- .r
    }
  }
  do.call(rbind, .ret)
}

#' Is a stress result a failure?
#'
#' @param res data frame from `stressRun()`
#' @return logical vector
#' @noRd
stressFailed <- function(res) {
  res$status %in% c("error", "problem", "not refused", "setup") &
    !(res$status == "setup" & res$expect == "any")
}

# ---------------------------------------------------------------------
# finding NONMEM/Monolix and the report
# ---------------------------------------------------------------------

#' The command that runs NONMEM
#'
#' @return the babelmixr2.nonmem option, an nmfe7* on the PATH or in
#'   the usual install directories, or "" when none is found
#' @noRd
stressFindNonmem <- function() {
  .o <- getOption("babelmixr2.nonmem", "")
  if (is.character(.o) && nzchar(.o)) {
    return(.o)
  }
  .w <- Sys.which(paste0("nmfe7", 9:0))
  .w <- .w[.w != ""]
  if (length(.w) > 0L) {
    return(unname(.w[1]))
  }
  .globs <- c(
    "/opt/NONMEM/*/run/nmfe7*",
    "/opt/nm*/run/nmfe7*",
    "/usr/local/NONMEM/*/run/nmfe7*",
    "/usr/local/nm*/run/nmfe7*",
    "~/nm*/run/nmfe7*",
    "~/NONMEM/*/run/nmfe7*",
    "C:/nm*/run/nmfe7*.bat",
    "C:/NONMEM/*/run/nmfe7*.bat"
  )
  .f <- Sys.glob(path.expand(.globs))
  .f <- .f[!grepl("\\.(f90|o|obj)$", .f) & file.exists(.f)]
  if (length(.f) == 0L) {
    return("")
  }
  # newest NONMEM first
  sort(.f, decreasing = TRUE)[1]
}

#' How Monolix is run
#'
#' @return description of the Monolix run command or lixoftConnectors,
#'   or "" when Monolix is not found
#' @noRd
stressMonolixStatus <- function() {
  .o <- getOption("babelmixr2.monolix", "")
  if (is.character(.o) && nzchar(.o)) {
    return(paste0("command ", .o))
  }
  if (!requireNamespace("lixoftConnectors", quietly = TRUE)) {
    return("")
  }
  .x <- try(
    suppressMessages(
      lixoftConnectors::initializeLixoftConnectors(
        software = "monolix",
        force = TRUE
      )
    ),
    silent = TRUE
  )
  if (inherits(.x, "try-error") || isFALSE(.x)) {
    return("")
  }
  paste0("lixoftConnectors ", utils::packageVersion("lixoftConnectors"))
}

#' Package versions for the report
#'
#' @return markdown list lines
#' @noRd
stressVersions <- function() {
  .p <- c(
    "babelmixr2",
    "rxode2",
    "nlmixr2est",
    "lotri",
    "nonmem2rx",
    "monolix2rx",
    "nlmixr2lib",
    "lixoftConnectors"
  )
  .v <- vapply(
    .p,
    function(p) {
      if (requireNamespace(p, quietly = TRUE)) {
        as.character(utils::packageVersion(p))
      } else {
        "-"
      }
    },
    character(1)
  )
  .sha <- utils::packageDescription("babelmixr2")$RemoteSha
  c(
    paste0(
      "- ",
      .p,
      " ",
      .v,
      ifelse(
        .p == "babelmixr2" & !is.null(.sha),
        paste0(" (", substr(.sha, 1, 7), ")"),
        ""
      )
    ),
    paste0("- R ", getRversion(), " on ", R.version$platform)
  )
}


#' Markdown summary of a stress test
#'
#' @param res results (from `stressRun()`, with a `failed` column)
#' @param modes modes that were run
#' @param engines engines that were used
#' @param nonmem,monolix how NONMEM/Monolix were run
#' @return markdown lines
#' @noRd
stressSummary <- function(res, modes, engines, nonmem, monolix) {
  .md <- c(
    "# babelmixr2 NONMEM/Monolix stress test",
    "",
    paste0("- date: ", format(Sys.time())),
    paste0("- mode: ", paste(modes, collapse = ", ")),
    stressVersions(),
    if ("run" %in% modes && "nonmem" %in% engines) {
      paste0("- NONMEM: ", nonmem)
    },
    if ("run" %in% modes && "monolix" %in% engines) {
      paste0("- Monolix: ", monolix)
    },
    "",
    "## Summary",
    ""
  )
  .col <- paste(res$mode, res$engine)
  .tab <- table(res$status, .col)
  .md <- c(
    .md,
    paste0("| status | ", paste(colnames(.tab), collapse = " | "), " |"),
    paste0("|---|", paste(rep("---", ncol(.tab)), collapse = "|"), "|"),
    vapply(
      rownames(.tab),
      function(r) {
        paste0("| ", r, " | ", paste(.tab[r, ], collapse = " | "), " |")
      },
      character(1)
    ),
    "",
    "## Failures",
    ""
  )
  .fail <- res[res$failed, ]
  if (nrow(.fail) == 0L) {
    .md <- c(.md, "None.")
  } else {
    .md <- c(
      .md,
      "| case | engine | mode | status | message |",
      "|---|---|---|---|---|",
      sprintf(
        "| %s | %s | %s | %s | %s |",
        .fail$case,
        .fail$engine,
        .fail$mode,
        .fail$status,
        gsub("\\|", "/", substr(paste(.fail$message, .fail$problems), 1, 300))
      )
    )
  }
  if ("run" %in% modes) {
    .ok <- res[res$mode == "run" & res$status %in% c("ok", "problem"), ]
    .num <- function(x, fmt) ifelse(is.na(x), "", sprintf(fmt, x))
    .md <- c(
      .md,
      "",
      "## Fits",
      "",
      paste(
        "| case | engine | status | seconds | objective | IPRED diff % |",
        "PRED diff % | rerun seconds | max rel. diff vs nlmixr2 |"
      ),
      "|---|---|---|---|---|---|---|---|---|",
      sprintf(
        "| %s | %s | %s | %.1f | %s | %s | %s | %s | %s |",
        .ok$case,
        .ok$engine,
        .ok$status,
        .ok$seconds,
        .num(.ok$objf, "%.3f"),
        .num(.ok$ipredRelDiff, "%.3f"),
        .num(.ok$predRelDiff, "%.3f"),
        .num(.ok$rerunSeconds, "%.1f"),
        .num(.ok$maxRelDiffTheta, "%.3f")
      )
    )
  }
  .md
}

#' Zip the output directory to send back
#'
#' @param out output directory
#' @return the zip file (invisibly), or NULL when it cannot be made
#' @noRd
stressBundle <- function(out) {
  .zip <- paste0(out, ".zip")
  .ok <- withr::with_dir(dirname(out), {
    try(utils::zip(.zip, basename(out), flags = "-r9Xq"), silent = TRUE)
  })
  if (inherits(.ok, "try-error") || !file.exists(.zip)) {
    message("could not zip the output; send the directory ", out, " instead")
    return(invisible(NULL))
  }
  message("send this file back: ", .zip)
  invisible(.zip)
}

# ---------------------------------------------------------------------
# the kit (run from an R session or from run-stress.R)
# ---------------------------------------------------------------------

#' Report the versions and whether NONMEM and Monolix are found
#'
#' @param nonmem command that runs NONMEM (like "nmfe75" or its full
#'   path); `NULL` looks for it (see `stressFindNonmem()`)
#' @param monolix command that runs Monolix; `NULL` uses
#'   lixoftConnectors (or `options(babelmixr2.monolix=)`)
#' @return list with `nonmem` and `monolix` (`""` when not found),
#'   invisibly
#' @noRd
stressCheck <- function(nonmem = NULL, monolix = NULL) {
  .found <- .stressEngines(nonmem, monolix)
  message(paste(stressVersions(), collapse = "\n"))
  message(
    "- NONMEM: ",
    if (nzchar(.found$nonmem)) {
      .found$nonmem
    } else {
      "not found (give it with nonmem=)"
    }
  )
  message(
    "- Monolix: ",
    if (nzchar(.found$monolix)) {
      .found$monolix
    } else {
      "not found (load lixoftConnectors or give it with monolix=)"
    }
  )
  if (!("linCmtMicro" %in% getNamespaceExports("rxode2"))) {
    message(
      "- rxode2 has no linCmtMicro(): the closed-form linCmt() checks will fail"
    )
  }
  invisible(.found)
}

.stressEngines <- function(nonmem = NULL, monolix = NULL) {
  list(
    nonmem = if (is.null(nonmem)) stressFindNonmem() else nonmem,
    monolix = if (is.null(monolix)) {
      stressMonolixStatus()
    } else {
      paste0("command ", monolix)
    },
    monolixCommand = monolix
  )
}

#' The cases to use
#'
#' @param cases regular expression of the case names (`NULL` is all)
#' @param nlmixr2lib "none", "sample" or "all" nlmixr2lib models
#' @return list of cases
#' @noRd
.stressSelect <- function(cases, nlmixr2lib) {
  .ret <- stressCases()
  if (nlmixr2lib != "none") {
    if (!requireNamespace("nlmixr2lib", quietly = TRUE)) {
      stop("nlmixr2lib= needs the nlmixr2lib package", call. = FALSE)
    }
    .ret <- c(.ret, stressLibCases(all = (nlmixr2lib == "all")))
  }
  if (!is.null(cases)) {
    .ret <- .ret[grepl(cases, names(.ret))]
  }
  .ret
}

#' List the stress cases
#'
#' @inheritParams .stressSelect
#' @return data frame with what each engine should do, invisibly
#' @noRd
stressList <- function(cases = NULL, nlmixr2lib = "none") {
  .c <- .stressSelect(cases, match.arg(nlmixr2lib, c("none", "sample", "all")))
  .e <- function(case, engine) {
    if (!(engine %in% case$engines)) {
      return("-")
    }
    .x <- case$expect[[engine]]
    if (.x %in% c("ok", "any")) .x else "refuse"
  }
  .ret <- data.frame(
    case = names(.c),
    nonmem = vapply(.c, .e, character(1), engine = "nonmem"),
    monolix = vapply(.c, .e, character(1), engine = "monolix"),
    description = vapply(.c, function(x) x$description, character(1)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  print(.ret, right = FALSE)
  invisible(.ret)
}

#' Run the stress test
#'
#' The kit: translate every case (and nlmixr2lib models), fit every case
#' with the engines that are found, compare with nlmixr2 and zip the
#' output to send back.  Works the same from an R session (like
#' RStudio) and from `run-stress.R`.
#'
#' @param engines engines to use; `NULL` uses the ones that are found
#'   (in run mode)
#' @param modes "translate" and/or "run"
#' @param cases regular expression of the case names (`NULL` is all)
#' @param nlmixr2lib "none", "sample" or "all" nlmixr2lib models
#'   (translation only); `NULL` is "sample" when nlmixr2lib is installed
#' @param out output directory
#' @param reference compare the fits with nlmixr2 (focei for NONMEM, saem
#'   for Monolix)
#' @param predTol largest median relative difference (%) between the
#'   rxode2 and NONMEM/Monolix IPRED
#' @param bundle zip the output directory
#' @inheritParams stressCheck
#' @return data frame of the results (invisibly), with the attributes
#'   `out` (output directory) and `zip` (the zip file)
#' @noRd
stressKit <- function(
  nonmem = NULL,
  monolix = NULL,
  engines = NULL,
  modes = c("translate", "run"),
  cases = NULL,
  nlmixr2lib = NULL,
  out = paste0("babelmixr2-stress-", format(Sys.time(), "%Y%m%d-%H%M%S")),
  reference = TRUE,
  predTol = 5,
  bundle = TRUE
) {
  modes <- match.arg(modes, c("translate", "run"), several.ok = TRUE)
  if (is.null(nlmixr2lib)) {
    nlmixr2lib <- if (requireNamespace("nlmixr2lib", quietly = TRUE)) {
      "sample"
    } else {
      "none"
    }
  }
  nlmixr2lib <- match.arg(nlmixr2lib, c("none", "sample", "all"))
  .found <- .stressEngines(nonmem, monolix)
  engines <- .stressKitEngines(engines, modes, .found)
  .cases <- .stressSelect(cases, "none")
  .libCases <- if (nlmixr2lib == "none") {
    list()
  } else {
    .stressSelect(cases, nlmixr2lib)
  }
  .libCases <- .libCases[setdiff(names(.libCases), names(.cases))]

  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  out <- normalizePath(out)
  writeLines(
    utils::capture.output(utils::sessionInfo()),
    file.path(out, "sessionInfo.txt")
  )
  message(paste(stressVersions(), collapse = "\n"))
  message(
    length(.cases) + length(.libCases),
    " cases; engines: ",
    paste(engines, collapse = ", "),
    "; mode: ",
    paste(modes, collapse = ", "),
    "; output: ",
    out
  )
  .runCommand <- list(
    nonmem = if (nzchar(.found$nonmem)) .found$nonmem,
    monolix = .found$monolixCommand
  )
  .res <- lapply(modes, function(mode) {
    # the nlmixr2lib models are translation only
    stressRun(
      if (mode == "translate") c(.cases, .libCases) else .cases,
      engines = engines,
      mode = mode,
      dir = if (length(modes) > 1L) file.path(out, mode) else out,
      runCommand = .runCommand,
      reference = reference,
      predTol = predTol,
      progress = TRUE
    )
  })
  .res <- do.call(rbind, .res)
  rownames(.res) <- NULL
  .res$failed <- stressFailed(.res)
  utils::write.csv(.res, file.path(out, "results.csv"), row.names = FALSE)
  .md <- stressSummary(.res, modes, engines, .found$nonmem, .found$monolix)
  writeLines(.md, file.path(out, "summary.md"))
  message("\n", paste(.md, collapse = "\n"))
  message("\nresults: ", file.path(out, "results.csv"))
  attr(.res, "out") <- out
  attr(.res, "zip") <- if (bundle) stressBundle(out)
  invisible(.res)
}

#' Which engines the kit uses
#'
#' @param engines engines asked for (`NULL`: the ones found in run mode,
#'   both in translate mode)
#' @param modes modes
#' @param found from `.stressEngines()`
#' @return engines
#' @noRd
.stressKitEngines <- function(engines, modes, found) {
  if (!("run" %in% modes)) {
    return(
      if (is.null(engines)) {
        c("nonmem", "monolix")
      } else {
        match.arg(engines, c("nonmem", "monolix"), several.ok = TRUE)
      }
    )
  }
  .found <- c(nonmem = nzchar(found$nonmem), monolix = nzchar(found$monolix))
  if (is.null(engines)) {
    engines <- names(.found)[.found]
    if (length(engines) == 0L) {
      stop("found neither NONMEM nor Monolix; see stressCheck()", call. = FALSE)
    }
  }
  engines <- match.arg(engines, c("nonmem", "monolix"), several.ok = TRUE)
  .how <- c(
    nonmem = "give it with nonmem= (like nonmem = \"nmfe75\")",
    monolix = "load lixoftConnectors or give it with monolix="
  )
  .missing <- engines[!.found[engines]]
  if (length(.missing) > 0L) {
    stop(
      paste0(
        c(nonmem = "NONMEM", monolix = "Monolix")[.missing[1]],
        " is not found; ",
        .how[.missing[1]]
      ),
      call. = FALSE
    )
  }
  engines
}
