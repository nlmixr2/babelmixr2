library(PopED)

library(babelmixr2)

f <- function() {
  ini({
    tEmax <- 100
    tED50 <- 20
    tGamma <- 4.5
    tBase <- 1
    eta.emax ~ 0.0625
    eta.ed50 ~ 0.0625
    eta.base ~ 0.0625
    prop.sd <- sqrt(0.01)
    add.sd <- sqrt(0.1)
  })
  model({
    EMAX <- tEmax*exp(eta.emax)
    ED50 <- tED50*exp(eta.ed50)
    GAMMA <- tGamma
    BASE <- tBase+eta.base
    y <- time
    DOSE <- time
    y <- BASE + EMAX*DOSE^(GAMMA)/(ED50^(GAMMA) + DOSE^(GAMMA))
    y ~ add(add.sd) + prop(prop.sd)
  })
}

e <- et(seq(0,50,length.out=8))

babel.db <- nlmixr2(f, e, "poped",
                    popedControl(groupsize=100,
                                 minxt=0,
                                 maxxt=50,
                                 ourzero=0))



library(ggplot2)
plot1 <- plot_model_prediction(babel.db,IPRED=T,DV=T)
plot1 + xlab("Dose")

## evaluate initial design

## $rse
##     EMAX       ED50      GAMMA       BASE     d_EMAX     d_ED50     d_BASE SIGMA[1,1]
## 2.588804   2.554668   1.112444   3.922208  14.825482  14.375072  33.257104   7.072694
## SIGMA[2,2]
## 18.487445
evaluate_design(babel.db)

# Optimization of doses
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output <- poped_optim(babel.db, opt_xt = T, parallel = T)

summary(output)
get_rse(output$FIM,output$poped.db)
plot_model_prediction(output$poped.db) + xlab("Dose")
