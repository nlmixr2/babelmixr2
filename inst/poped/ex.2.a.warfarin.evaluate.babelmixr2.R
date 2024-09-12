## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation
##   for population pharmacokinetics-pharmacodynamics studies",
##   Br. J. Clin. Pharm., 2014.

library(babelmixr2)
library(PopED)

##-- Model: One comp first order absorption

f <- function() {
  ini({
    tCl <- 0.15
    tV <- 8
    tKA <- 1.0
    tFavail <- fix(1)
    eta.cl ~ 0.07
    eta.v ~ 0.02
    eta.ka ~ 0.6

    prop.sd <- sqrt(0.01)
  })
  model({
    CL <- tCl*exp(eta.cl)
    V <- tV*exp(eta.v)
    KA <- tKA*exp(eta.ka)
    Favail <- tFavail
    y <- (DOSE*Favail*KA/(V*(KA-CL/V)))*(exp(-CL/V*time)-exp(-KA*time))
    y ~ prop(prop.sd)

  })
}

e <-  et(c(0.5, 1,2,6,24,36,72,120)) %>%
  as.data.frame()

## -- Define initial design  and design space
babel.db <- nlmixr2(f, e, "poped",
                    control=popedControl(
                      groupsize=32,
                      minxt=0,
                      maxxt=120,
                      a=70))

##  create plot of model without variability
plot_model_prediction(babel.db)

##  create plot of model with variability
plot_model_prediction(babel.db,IPRED=T,DV=T)

## get predictions from model
model_prediction(babel.db)

## evaluate initial design
evaluate_design(babel.db)
shrinkage(babel.db)

## Evaluate with full FIM
evaluate_design(babel.db, fim.calc.type=0)

# Examine efficiency of sampling windows
plot_efficiency_of_windows(babel.db,xt_windows=0.5)
