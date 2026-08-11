## HCV example (Nyberg et al., Br. J. Clin. Pharm., 2014) with the dose
## amount *and* the infusion duration treated as PopED design ("a")
## variables.  See https://github.com/nlmixr2/babelmixr2/issues/131
##
## The trick is to keep the dose *record* in the design dataset, but to
## make both of its interesting quantities model quantities:
##
##  - `amt=1`   in the data + `f(depot)   <- DOSE` in the model
##      => the administered amount is the covariate `DOSE`
##  - `rate=-2` in the data + `dur(depot) <- TINF` in the model
##      => the infusion duration is the covariate `TINF`
##
## Because `DOSE` and `TINF` are ordinary covariates as far as
## rxode2/babelmixr2 are concerned, PopED sees them as elements of `a`
## and can optimize over them with `opt_a=TRUE`.
##
## The dosing *interval* is still fixed by the event table here.  See
## ex.10.PKPD.HCV.adaptive-dosing.babelmixr2.R for the variant that
## pushes the whole regimen from inside the model with `infuseDur()`,
## which makes `ii` a design variable too.

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

    d/dt(depot)   <- -KA*depot
    d/dt(central) <-  KA*depot - KE*central
    d/dt(TC) <- s - TC*(e*VP + d)              # target cells (TC)
    d/dt(IC) <- e*TC*VP - delta*IC             # productively infected cells
    d/dt(VP) <- p*(1 - (pow(central/VD, n)/(pow(central/VD, n) +
                                              pow(EC50, n))))*IC - c*VP

    ## the two design variables enter here
    f(depot)   <- DOSE  # the dataset carries amt=1
    dur(depot) <- TINF  # the dataset carries rate=-2

    TC(0) <- c*delta/(p*e)
    IC(0) <- (s*e*p - d*c*delta)/(p*delta*e)
    VP(0) <- (s*e*p - d*c*delta)/(c*delta*e)

    conc <- central/VD
    eff  <- log10(VP)
    conc ~ add(add.sd.pk)
    eff  ~ add(add.sd.pd)
  })
}

TAU <- 7
tms <- c(0, 0.25, 0.5, 1, 2, 3, 4, 7, 10, 14, 21, 28)

## note `amt=1` and `rate=-2`; the actual amount/duration come from the model
e1 <- et(amt=1, rate=-2, ii=TAU, addl=3, cmt="depot") |>
  et(tms) |>
  as.data.frame()
e1$dvid <- 1

e2 <- e1[e1$evid == 0, ]
e2$dvid <- 2

e <- rbind(e1, e2)

babel.db <- nlmixr2(f, e, "poped",
                    popedControl(groupsize=30,
                                 a=list(c(DOSE=180, TINF=1)),
                                 mina=c(DOSE=20,  TINF=0.1),
                                 maxa=c(DOSE=600, TINF=12)))

plot_model_prediction(babel.db, facet_scales="free")

evaluate_design(babel.db)
#> $ofv
#> [1] 88.26974

## the criterion really does depend on the infusion duration
for (ti in c(0.1, 0.5, 1, 2, 4, 8, 12)) {
  .db <- babel.db
  .db$design$a[1, "TINF"] <- ti
  message(sprintf("TINF = %5.2f   ofv = %.4f", ti, evaluate_design(.db)$ofv))
}

## ...so it can be optimized over
r <- poped_optim(babel.db, opt_a=TRUE, opt_xt=FALSE, parallel=FALSE)
#> Optimized Covariates:
#> Group 1: 1 : 446.115 : 0.456544
#>
#> OFV = 95.5746
#
# (the default algorithm is a random search, so the exact DOSE/TINF the
# optimizer lands on varies from run to run; the point is that both are
# now optimized instead of being fixed features of the event table)

## Note: `cmt` in the design dataset must be given by *name*
## (`cmt="depot"`); a bare compartment number is not translated in the
## PopED path and silently produces a design with no dosing.
