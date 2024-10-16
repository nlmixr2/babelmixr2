## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

library(PopED)
# library(deSolve)
library(babelmixr2)

##-- Model: One comp first order absorption
##-- Model: One comp first order absorption, analytic solution
f_without <- function() {
      ini({
            tV <- 8
            tKa <- 1
            tCl <- 0.15
            tF <- fix(1)
            
            eta.v ~ 0.02
            eta.ka ~ 0.6
            eta.cl ~0.07
            
            prop.sd <- sqrt(0.01)
            
      })
      model({
            V<-tV*exp(eta.v)
            KA<-tKa*exp(eta.ka)
            CL<-tCl*exp(eta.cl)
            
            KE=CL/V
            xt=t # set xt as time
            Favail <- tF
            
            y <- (DOSE*Favail*KA/(V*(KA-KE)))*(exp(-KE*xt)-exp(-KA*xt))
            
            y ~ prop(prop.sd)
      })
}

f_with <- function() {
      ini({
            tV <- 8
            tKa <- 1
            tCl <- 0.15
            tF <- fix(1)
            
            ## For correlated parameters, you specify the names of each
            ## correlated parameter separated by a addition operator `+`
            ## and the left handed side specifies the lower triangular
            ## matrix initial of the covariance matrix.
            eta.cl + eta.v + eta.ka ~ c(0.07,
                                       0.03, 0.02,
                                       0.1,  0.09, 0.6)
            
            prop.sd <- sqrt(0.01)
            
      })
      model({
            V<-tV*exp(eta.v)
            KA<-tKa*exp(eta.ka)
            CL<-tCl*exp(eta.cl)
            
            KE=CL/V
            xt=t # set xt as time
            Favail <- tF
            
            y <- (DOSE*Favail*KA/(V*(KA-KE)))*(exp(-KE*xt)-exp(-KA*xt))
            
            y ~ prop(prop.sd)
      })
}

# minxt, maxxt
e <- et(list(c(0, 120),
             c(0, 120),
             c(0, 120),
             c(0, 120),
             c(0, 120),
             c(0, 120),
             c(0, 120),
             c(0, 120))) %>%
      as.data.frame()

#xt
e$time <-  c(0.5,1,2,6,24,36,72,120)

babel.db.without <- nlmixr2(f_without, e, "poped",
                            popedControl(groupsize=32,
                                   bUseGrouped_xt=TRUE,
                                   a=c(DOSE=70),
                                   maxa=c(DOSE=200),
                                   mina=c(DOSE=0)))

babel.db.with <- nlmixr2(f_with, e, "poped",
                         popedControl(groupsize=32,
                                         bUseGrouped_xt=TRUE,
                                         a=c(DOSE=70),
                                         maxa=c(DOSE=200),
                                         mina=c(DOSE=0)))


# What do the covariances mean?
(IIV <- babel.db.with$parameters$param.pt.val$d)
cov2cor(IIV)

##  create plot of model with variability 
library(ggplot2)
p1 <- plot_model_prediction(babel.db.without,IPRED=T)+ylim(0,12)
p2 <- plot_model_prediction(babel.db.with,IPRED=T) +ylim(0,12)
gridExtra::grid.arrange(p1, p2, nrow = 1)

## evaluate initial design
evaluate_design(babel.db.without)
evaluate_design(babel.db.with)

## Evaluate with full FIM
evaluate_design(babel.db.without, fim.calc.type=0)
evaluate_design(babel.db.with, fim.calc.type=0)

