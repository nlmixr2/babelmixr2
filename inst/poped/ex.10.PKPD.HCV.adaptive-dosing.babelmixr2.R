## HCV example (Nyberg et al., Br. J. Clin. Pharm., 2014) where the whole
## regimen -- dose amount, infusion duration *and* dosing interval -- is
## registered from inside the model with rxode2's adaptive dosing
## functions, and all three are PopED design ("a") variables.
## See https://github.com/nlmixr2/babelmixr2/issues/131
##
## Compare with ex.10.PKPD.HCV.dose-and-tinf.babelmixr2.R, which gets the
## amount and the duration out of a dose record with `f()`/`dur()`.  That
## approach is simpler and needs nothing from the solver, but the dosing
## *interval* stays baked into the event table.  Here the design dataset
## carries no dose records at all: `infuseDur()` pushes them, so `ii` is
## just another covariate and can be optimized too.
##
## Requires rxode2 > 5.1.7 (the pushed-dose timing fix, rxode2#1214).

library(babelmixr2)
library(PopED)

f <- function() {
  ini({
    tp <- fix(100)
    td <- fix(0.001)
    te <- fix(1e-7)
    ts <- fix(20000)

    tKA <- log(0.8)
    tKE <- log(0.15)
    tVD <- log(100)
    tEC50 <- log(0.12)
    tn <- log(2)
    tdelta <- log(0.2)
    tc <- log(7)

    eta.KA ~ 0.25
    eta.KE ~ 0.25
    eta.VD ~ 0.25
    eta.EC50 ~ 0.25
    eta.n ~ 0.25
    eta.delta ~ 0.25
    eta.c ~ 0.25

    add.sd.pk <- sqrt(0.04) # nlmixr2 uses sd
    add.sd.pd <- sqrt(0.04)
  })
  model({
    p <- tp
    d <- td
    e <- te
    s <- ts
    KA <- exp(tKA + eta.KA)
    KE <- exp(tKE + eta.KE)
    VD <- exp(tVD + eta.VD)
    EC50 <- exp(tEC50 + eta.EC50)
    n <- exp(tn + eta.n)
    delta <- exp(tdelta + eta.delta)
    c <- exp(tc + eta.c)

    ## The design variables have to be copied into model variables before
    ## they are handed to infuseDur(); a covariate that appears *only* as
    ## an adaptive dosing argument is not registered as a parameter and
    ## the model will not compile (rxode2#1231).
    amtI <- DOSE
    durI <- TINF
    iiI  <- TAU

    d/dt(depot)   <- -KA*depot
    d/dt(central) <-  KA*depot - KE*central
    d/dt(TC) <- s - TC*(e*VP + d)              # target cells (TC)
    d/dt(IC) <- e*TC*VP - delta*IC             # productively infected cells
    d/dt(VP) <- p*(1 - (pow(central/VD, n)/(pow(central/VD, n) +
                                              pow(EC50, n))))*IC - c*VP

    TC(0) <- c*delta/(p*e)
    IC(0) <- (s*e*p - d*c*delta)/(p*delta*e)
    VP(0) <- (s*e*p - d*c*delta)/(c*delta*e)

    ## four infusions of DOSE over TINF, every TAU
    if (t <= 0) {
      infuseDur(amtI, durI, cmt=depot, ii=iiI, addl=3)
    }

    conc <- central/VD
    eff  <- log10(VP)
    conc ~ add(add.sd.pk)
    eff  ~ add(add.sd.pd)
  })
}

## observation-only design: there are no dose records to write down
tms <- c(0, 0.25, 0.5, 1, 2, 3, 4, 7, 10, 14, 21, 28)
e1 <- as.data.frame(et(tms))
e1$dvid <- 1
e2 <- e1
e2$dvid <- 2
e <- rbind(e1, e2)

babel.db <- nlmixr2(f, e, "poped",
                    popedControl(groupsize=30,
                                 a=list(c(DOSE=180, TINF=1, TAU=7)),
                                 mina=c(DOSE=20,  TINF=0.1, TAU=1),
                                 maxa=c(DOSE=600, TINF=12,  TAU=14)))

plot_model_prediction(babel.db, facet_scales="free")

evaluate_design(babel.db)$ofv
#> [1] 88.27007
## ...which is the same design as the `f()`/`dur()` version of this
## example (88.26974); the two parameterizations agree.

## the interval is now a design variable like any other
for (ta in c(1, 2, 3, 5, 7, 10, 14)) {
  .db <- babel.db
  .db$design$a[1, "TAU"] <- ta
  message(sprintf("TAU = %5.1f   ofv = %.4f", ta, evaluate_design(.db)$ofv))
}
#> TAU =   1.0   ofv = 91.0976
#> TAU =   2.0   ofv = 90.6507
#> TAU =   3.0   ofv = 90.2115
#> TAU =   5.0   ofv = 89.4562
#> TAU =   7.0   ofv = 88.2701
#> TAU =  10.0   ofv = 89.1156
#> TAU =  14.0   ofv = 88.8393

r <- poped_optim(babel.db, opt_a=TRUE, opt_xt=FALSE, parallel=FALSE)
#> Optimized Covariates:
#> Group 1: 1 : 600 : 7.1429 : 1
#>
#> OFV = 100.785
#
# (a random search, so the exact optimum moves between runs)
