## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation
##   for population pharmacokinetics-pharmacodynamics studies",
##   Br. J. Clin. Pharm., 2014.

## Optimization using an additive + proportional reidual error to
##   avoid sample times at very low concentrations (time 0 or very late samoples).

## Model described with an ODE
library(PopED)
library(babelmixr2)

f <- function() {
  ini({
    tCl <- 0.15
    tV <- 8
    tKA <- 1.0
    tFavail <- fix(1)
    eta.cl ~ 0.07
    eta.v ~ 0.02
    eta.ka ~ 0.6
    prop.sd <- sqrt(0.01) # nlmixr2 uses sd
    add.sd <- sqrt(0.25)
  })
  model({
    CL <- tCl*exp(eta.cl)
    V <- tV*exp(eta.v)
    KA <- tKA*exp(eta.ka)
    Favail <- tFavail
    d/dt(depot) <- -KA*depot
    d/dt(central) <- KA*depot - (CL/V)*central
    depot(0) <- Favail*DOSE
    y <- central/V
    y ~ prop(prop.sd) + add(add.sd)
  })
}

## -- Define initial design  and design space
e <-  et(c(0.5, 1,2,6,24,36,72,120)) %>%
  as.data.frame()

babel.db <- nlmixr2(f, e, "poped",
                    popedControl(groupsize=32,
                                 minxt=0,
                                 maxxt=120,
                                 a=70,
                                 mina=0,
                                 maxa=100))

##  create plot of model without variability
plot_model_prediction(babel.db)

##  create plot of model with variability
plot_model_prediction(babel.db,IPRED=T,DV=T)

## evaluate initial design (much faster than pure R solution)
tic(); design_ode_compiled <- evaluate_design(babel.db); toc()

## making optimization times more reasonable
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output <- poped_optim(babel.db, opt_xt =TRUE, parallel=TRUE, method = "LS")
