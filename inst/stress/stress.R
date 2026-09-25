# babelmixr2 NONMEM/Monolix stress test cases and runner
#
# This file is shared by the testthat stress tests (translation only)
# and by `run-stress.R` (which can also run NONMEM/Monolix end to end).
# See README.md in this directory for how to run it.
#
# Every case is a model, a data set and, for each engine, what should
# happen:
#
#  - "ok": the model and data are translated (and fit in "run" mode)
#  - a regular expression: the model is refused with an error matching it
#
# Cases can also list regular expressions that must be found in the
# NONMEM control stream (`nonmem`) or the Monolix model file
# (`monolix`), which check that the edge case was translated the way
# it should be.

.stressEnv <- new.env(parent=emptyenv())

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
  }
)

# ---------------------------------------------------------------------
# cases
# ---------------------------------------------------------------------

.stressCase <- function(name, model, data, nonmem="ok", monolix="ok",
                        checkNonmem=NULL, checkMonolix=NULL, run=TRUE,
                        description="") {
  list(name=name, model=model, data=data,
       expect=list(nonmem=nonmem, monolix=monolix),
       check=list(nonmem=checkNonmem, monolix=checkMonolix),
       run=run, description=description)
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

  .res <- function(err, ini) .stressResidual(err, ini)

  .cases <- list(
    # linCmt() closed form ------------------------------------------
    .stressCase("linCmt 1-cmt oral", .m$lin1oral, .theo,
                checkNonmem=c("ADVAN2 TRANS1", "K=CL/V"),
                checkMonolix=c("pkmodel\\(V=rx_v, k=rx_k, ka=rx_ka\\)"),
                description="linCmt() ~ endpoint; NONMEM ADVAN2, Monolix pkmodel()"),
    .stressCase("linCmt 1-cmt oral lag and bioavailability", .m$lin1oralLagF, .theo,
                checkNonmem=c("ADVAN2 TRANS1", "ALAG1=", "F1="),
                checkMonolix=c("Tlag=rx_tlag", "p=rx_p")),
    .stressCase("linCmt 1-cmt iv bolus", .m$lin1iv, .ivBolus,
                checkNonmem=c("ADVAN1 TRANS1", "RXE_CP=A\\(1\\)/"),
                checkMonolix=c("pkmodel\\(V=rx_v, k=rx_k\\)")),
    .stressCase("linCmt 1-cmt iv infusion (rate)", .m$lin1iv, .ivInf,
                checkNonmem=c("ADVAN1 TRANS1", "\\$INPUT.* RATE"),
                checkMonolix=c("pkmodel\\(")),
    .stressCase("linCmt 1-cmt iv infusion (duration)", .m$lin1iv, .ivDurData,
                nonmem="duration/tinf",
                checkMonolix=c("pkmodel\\(")),
    .stressCase("linCmt 1-cmt iv k/v parameterization", .m$lin1ivK, .ivBolus,
                checkNonmem=c("ADVAN1 TRANS1"),
                checkMonolix=c("rx_k = k")),
    .stressCase("linCmt 2-cmt oral", .m$lin2oral, .oral2,
                checkNonmem=c("ADVAN4 TRANS1", "K23=", "K32=", "RXE_CP=A\\(2\\)/"),
                checkMonolix=c("k12=rx_k12", "ka=rx_ka")),
    .stressCase("linCmt 2-cmt iv", .m$lin2iv, .iv2,
                checkNonmem=c("ADVAN3 TRANS1", "K12=", "K21="),
                checkMonolix=c("k21=rx_k21")),
    .stressCase("linCmt 2-cmt iv micro-constants", .m$lin2ivMicro, .iv2,
                checkNonmem=c("ADVAN3 TRANS1"),
                checkMonolix=c("pkmodel\\(")),
    .stressCase("linCmt 3-cmt oral", .m$lin3oral, .oral3,
                checkNonmem=c("ADVAN12 TRANS1", "K24=", "K42="),
                checkMonolix=c("k13=rx_k13", "k31=rx_k31")),
    .stressCase("linCmt 3-cmt iv", .m$lin3iv, .iv3,
                checkNonmem=c("ADVAN11 TRANS1", "K13=", "K31="),
                checkMonolix=c("k13=rx_k13")),
    .stressCase("linCmt modeled duration", .m$lin1ivDur, .ivModelDur,
                checkNonmem=c("ADVAN1 TRANS1", "D1=", "\\$INPUT.* RATE"),
                checkMonolix=c("ddt_central", "Tk0="),
                description="Monolix pkmodel() cannot model the duration: ODEs"),
    .stressCase("linCmt modeled rate", .m$lin1ivRate, .ivModelRate,
                checkNonmem=c("ADVAN1 TRANS1", "R1=", "\\$INPUT.* RATE"),
                checkMonolix=c("ddt_central", "Tk0=amtDose/"),
                description="Monolix pkmodel() cannot model the rate: ODEs"),
    .stressCase("linCmt weight covariate", .m$lin1oralWt, .theoWt,
                checkNonmem=c("ADVAN2 TRANS1", "\\$INPUT.* NLMIXRMUDERCOV1",
                              "\\+NLMIXRMUDERCOV1\\*THETA"),
                checkMonolix=c("pkmodel\\(")),
    .stressCase("linCmt with effect compartment ODE", .m$lin1oralEffect, .effect,
                checkNonmem=c("ADVAN13", "DADT\\(3\\)"),
                checkMonolix=c("ddt_ce", "ddt_central"),
                description="mixed linCmt()/ODE models are solved as ODEs"),
    .stressCase("linCmt with ODE defined first", .m$lin1oralOdeFirst, .effectFirst,
                checkNonmem=c("ADVAN13", "DADT\\(3\\) = KE0", "DADT\\(1\\) = - \\(KA\\)\\*A\\(1\\)"),
                checkMonolix=c("ddt_ce"),
                description="compartment numbers stay depot=1, central=2, ce=3"),
    # nlmixr2est's mu2 covariate hook reads `t` as a data column (and
    # finds base::t()); when that is fixed upstream this case reports
    # "not refused" and the expectation should become "ok"
    .stressCase("linCmt time-dependent clearance", .m$lin1oralTime, .theo,
                nonmem="replicate an object of type 'closure'",
                monolix="replicate an object of type 'closure'",
                checkNonmem=c("ADVAN13"),
                checkMonolix=c("ddt_central"),
                description="parameters changing with time need ODEs (known nlmixr2est mu2 issue with t)"),
    .stressCase("linCmt amount in the endpoint", .m$lin1oralAmount, .theo,
                checkNonmem=c("ADVAN2 TRANS1", "A\\(2\\)"),
                checkMonolix=c("ddt_central"),
                description="Monolix pkmodel() has no amounts: ODEs"),
    .stressCase("linCmt steady state", .m$lin1oral, .ss,
                checkNonmem=c("ADVAN2 TRANS1", "\\$INPUT.* II .* SS"),
                checkMonolix=c("pkmodel\\(")),
    .stressCase("linCmt steady state with lag time", .m$lin1oralLagF, .ss,
                checkNonmem=c("ADVAN2 TRANS1", "\\$INPUT.* II .* SS", "ALAG1="),
                checkMonolix=c("Tlag=rx_tlag"),
                description="rxode2 splits a lagged steady state dose in two; NONMEM/Monolix get one SS dose"),
    .stressCase("linCmt additional doses", .m$lin1oral, .addl,
                checkNonmem=c("ADVAN2 TRANS1"),
                checkMonolix=c("pkmodel\\(")),
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
                .theoPos, checkMonolix=c("logNormal|lognormal|exponential")),
    .stressCase("residual logit",
                .res("logitNorm(logit.sd, 0, 20)", "logit.sd <- 0.1"),
                .theoPos, checkMonolix=c("logitNormal, min=0, max=20")),
    .stressCase("residual boxCox",
                .res("add(add.sd) + boxCox(lambda)", "add.sd <- 0.7; lambda <- 0.5"),
                .theoPos, monolix="boxCox|yeoJohnson|not supported|transform"),
    .stressCase("residual yeoJohnson",
                .res("add(add.sd) + yeoJohnson(lambda)", "add.sd <- 0.7; lambda <- 0.5"),
                .theoPos, monolix="boxCox|yeoJohnson|not supported|transform"),
    .stressCase("residual t distribution", .stressCode$tDist, .theoPos,
                nonmem="t|not supported|normal", monolix="t|not supported|normal"),
    .stressCase("fixed residual error", .stressCode$fixedResidual, .theoPos,
                checkNonmem=c("FIX")),
    # censoring --------------------------------------------------------
    .stressCase("censoring (CENS)", .m$lin1oral, .cens,
                checkNonmem=c("CENS")),
    .stressCase("censoring (CENS and LIMIT)", .m$lin1oral, .censLimit,
                checkNonmem=c("LIMIT")),
    .stressCase("censoring LIMIT only", .m$lin1oral, .limitOnly,
                nonmem="ylo|yup|limit|laplacian"),
    # model code ---------------------------------------------------------
    .stressCase("if/else", .stressCode$ifElse, .theoSex, checkNonmem=c("IF \\(")),
    .stressCase("else if", .stressCode$ifElseIf, .theoSex, nonmem="else|if"),
    .stressCase("NONMEM reserved names", .stressCode$reserved, .theo,
                checkNonmem=c("RXR1")),
    .stressCase("long dotted names", .stressCode$longNames, .theo),
    .stressCase("fixed and block random effects", .stressCode$fixedAndBlock, .theo,
                checkNonmem=c("BLOCK\\(2\\)", "FIX")),
    .stressCase("parameters without random effects", .stressCode$noEta, .theo),
    .stressCase("initial condition", .stressCode$initialCondition, .ivBolus,
                checkNonmem=c("A_0\\(1\\)")),
    .stressCase("time-varying covariate", .stressCode$timeVaryingCov, .theoCrcl,
                checkMonolix=c("CRCL")),
    .stressCase("probitInv", .stressCode$probitInv, .theo),
    .stressCase("between-occasion variability", .stressCode$iov, .theoOcc,
                nonmem="id|occasion|level|random", monolix="id|occasion|level|random",
                run=FALSE),
    .stressCase("not mu-referenced", .stressCode$notMuRef, .theo,
                monolix="mu"),
    # data -------------------------------------------------------------
    .stressCase("reset and dose (evid=4)", .m$ode1, .reset)
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

.stressControl <- function(engine, mode, modelName, runCommand) {
  .cmd <- if (mode == "translate") NA else runCommand
  if (engine == "nonmem") {
    if (is.null(.cmd)) .cmd <- getOption("babelmixr2.nonmem", "")
    babelmixr2::nonmemControl(runCommand=.cmd, modelName=modelName)
  } else {
    if (is.null(.cmd)) .cmd <- getOption("babelmixr2.monolix", "")
    babelmixr2::monolixControl(runCommand=.cmd, modelName=modelName)
  }
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
#' @return one row data frame with the result
#' @noRd
stressRunCase <- function(case, engine, mode="translate", dir=tempfile("stress"),
                          runCommand=NULL, reference=FALSE) {
  .expect <- case$expect[[engine]]
  .modelName <- .stressModelName(case$name)
  .ret <- data.frame(case=case$name, engine=engine, mode=mode,
                     expect=.expect, status=NA_character_, message="",
                     problems="", seconds=NA_real_, objf=NA_real_,
                     maxRelDiffTheta=NA_real_, stringsAsFactors=FALSE)
  if (inherits(case$model, "try-error") || inherits(case$data, "try-error")) {
    .ret$status <- "setup"
    .ret$message <- paste(c(if (inherits(case$model, "try-error")) paste(case$model),
                            if (inherits(case$data, "try-error")) paste(case$data)),
                          collapse=" ")
    return(.ret)
  }
  if (mode == "run" && !isTRUE(case$run)) {
    .ret$status <- "skipped"
    .ret$message <- "case is translation only"
    return(.ret)
  }
  dir.create(dir, showWarnings=FALSE, recursive=TRUE)
  .time <- proc.time()
  .fit <- withr::with_dir(dir, {
    tryCatch(suppressWarnings(suppressMessages(
      nlmixr2est::nlmixr2(case$model, case$data, est=engine,
                          control=.stressControl(engine, mode, .modelName, runCommand)))),
      error=function(e) e)
  })
  .ret$seconds <- (proc.time() - .time)[["elapsed"]]
  if (inherits(.fit, "error")) {
    .msg <- conditionMessage(.fit)
    .ret$message <- gsub("\n", " ", .msg)
    if (identical(.expect, "ok")) {
      .ret$status <- "error"
    } else if (identical(.expect, "any")) {
      .known <- any(vapply(.stressKnownRefusals, function(r) grepl(r, .msg), logical(1)))
      .upstream <- any(vapply(.stressKnownUpstream, function(r) grepl(r, .msg, fixed=TRUE), logical(1)))
      .ret$status <- if (.known) "refused" else if (.upstream) "upstream" else "error"
    } else {
      .ret$status <- if (grepl(.expect, .msg)) "refused" else "error"
    }
    return(.ret)
  }
  if (!(.expect %in% c("ok", "any"))) {
    .ret$status <- "not refused"
    .ret$message <- paste0("expected an error matching '", .expect, "'")
    return(.ret)
  }
  .files <- .stressFiles(engine, .modelName)
  .problems <- character(0)
  .modelFile <- file.path(dir, .files$model)
  if (!file.exists(.modelFile)) {
    .problems <- c(.problems, paste0("missing ", .files$model))
  } else {
    .lines <- readLines(.modelFile)
    .problems <- c(.problems,
                   if (engine == "nonmem") stressLintNonmem(.lines) else stressLintMonolix(.lines))
    if (engine == "monolix") {
      .mlxtran <- file.path(dir, .files$mlxtran)
      if (file.exists(.mlxtran)) {
        .lines <- c(.lines, readLines(.mlxtran))
        .problems <- c(.problems, .stressLintMlxtran(readLines(.mlxtran)))
      }
    }
    for (.re in case$check[[engine]]) {
      if (!any(grepl(.re, .lines))) {
        .problems <- c(.problems, paste0("missing '", .re, "'"))
      }
    }
  }
  .dataFile <- file.path(dir, .files$data)
  if (!file.exists(.dataFile)) {
    .problems <- c(.problems, paste0("missing ", .files$data))
  } else {
    .problems <- c(.problems, .stressLintData(utils::read.csv(.dataFile, na.strings=".")))
  }
  .ret$problems <- paste(.problems, collapse=" | ")
  .ret$status <- if (length(.problems) > 0L) "problem" else "ok"
  if (mode == "run" && inherits(.fit, "nlmixr2FitData")) {
    .ret$objf <- .fit$objf
    if (reference) {
      .ret$maxRelDiffTheta <- .stressCompare(.fit, case, engine)
    }
  }
  .ret
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
                      reference=FALSE, progress=interactive()) {
  .ret <- list()
  for (.case in cases) {
    for (.engine in engines) {
      .dir <- file.path(dir, .engine, .stressModelName(.case$name))
      .r <- stressRunCase(.case, .engine, mode=mode, dir=.dir,
                          runCommand=runCommand[[.engine]], reference=reference)
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
