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

#########################################
## NOTE All PopED output for residuals
## (add or prop) are VARIANCES instead of
## standard deviations!
#########################################

## get predictions from model
## Original:
## > model_prediction(poped.db)
## >    Time      PRED Group Model a_i
## 1   0.5 3.4254357     1     1  70
## 2   1.0 5.4711041     1     1  70
## 3   2.0 7.3821834     1     1  70
## 4   6.0 7.9462805     1     1  70
## 5  24.0 5.6858561     1     1  70
## 6  36.0 4.5402483     1     1  70
## 7  72.0 2.3116966     1     1  70
## 8 120.0 0.9398657     1     1  70
model_prediction(babel.db)

## evaluate initial design
# Original:
## $rse
## CL          V         KA       d_CL        d_V       d_KA SIGMA[1,1]
## 4.738266   2.756206  13.925829  25.627205  30.344316  25.777327  11.170784
evaluate_design(babel.db)

shrinkage(babel.db)

## Evaluate with full FIM
evaluate_design(babel.db, fim.calc.type=0)

# Examine efficiency of sampling windows
plot_efficiency_of_windows(babel.db,xt_windows=0.5)
