## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation
##   for population pharmacokinetics-pharmacodynamics studies",
##   Br. J. Clin. Pharm., 2014.

## Optimization using an additive + proportional residual error to
##   avoid sample times at very low concentrations (time 0 or very late samples).

library(babelmixr2)
library(PopED)

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
    y <- (DOSE*Favail*KA/(V*(KA-CL/V)))*(exp(-CL/V*time)-exp(-KA*time))
    y ~ prop(prop.sd) + add(add.sd)
  })
}
# First define standard controler from nlmixr2
e <-  et(c(0.5, 1,2,6,24,36,72,120)) %>%
  as.data.frame()

babel.db <- nlmixr2(f, e, "poped",
                    popedControl(groupsize=32,
                                 minxt=0,
                                 maxxt=120,
                                 a=70,
                                 mina=0,
                                 maxa=100,
                                 # selecting important/unimportant
                                 # parameters assumes Ds optimal design.
                                 important=c("tCl", "tV", "tKa")))

##  create plot of model without variability
plot_model_prediction(babel.db)

##  create plot of model with variability
plot_model_prediction(babel.db,IPRED=T,DV=T)

# Original RSEs
## $rse
##       CL          V         KA       d_CL        d_V       d_KA SIGMA[1,1] SIGMA[2,2]
## 5.096246   3.031164  14.260384  29.761226  36.681388  26.748640  32.011719  25.637971

## evaluate initial design
evaluate_design(babel.db)

# RS+SG+LS optimization of sample times
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output <- poped_optim(babel.db, opt_xt=T, parallel=T)

summary(output)

plot_model_prediction(output$poped.db)
