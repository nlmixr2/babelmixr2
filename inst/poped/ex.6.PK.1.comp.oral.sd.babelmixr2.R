library(PopED)
library(babelmixr2)

##-- Model: One comp first order absorption
##-- Model: One comp first order absorption, analytic solution
f1 <- function() {
      ini({
            tV <- 8
            tKa <- 1
            tCl <- 0.15
            tF <- fix(1)
            
            eta.v ~ 0.02
            eta.ka ~ 0.6
            eta.cl ~0.07
            
            prop.sd <- sqrt(0.01)
            add.sd <- sqrt(0.25)
            
      })
      model({
            V<-tV*exp(eta.v)
            KA<-tKa*exp(eta.ka)
            CL<-tCl*exp(eta.cl)
            
            KE=CL/V
            xt=t # set xt as time
            Favail <- tF
            
            y <- (DOSE*Favail*KA/(V*(KA-KE)))*(exp(-KE*xt)-exp(-KA*xt))
            
            y ~ prop(prop.sd) + add(add.sd)
      })
}

##-- Model: One comp first order absorption, analytic solution
## reparameterization
f2 <- function() {
      ini({
            tV <- 8
            tKa <- 1
            tKe <- 0.15/8
            tF <- fix(1)
            
            eta.v ~ 0.02
            eta.ka ~ 0.6
            eta.ke ~0.07
            
            prop.sd <- sqrt(0.01)
            add.sd <- sqrt(0.25)
            
      })
      model({
            V<-tV*exp(eta.v)
            KA<-tKa*exp(eta.ka)
            KE=tKe*exp(eta.ke)
            
            Favail <- tF
            xt=t # set xt as time
            
            y <- (DOSE*Favail*KA/(V*(KA-KE)))*(exp(-KE*xt)-exp(-KA*xt))
            
            y ~ prop(prop.sd) + add(add.sd)
      })
}

# ODE solution
f3 <- function() {
      ini({
            tV <- 8
            tKa <- 1
            tCl <- 0.15
            tF <- fix(1)
            
            eta.v ~ 0.02
            eta.ka ~ 0.6
            eta.cl ~0.07
            
            prop.sd <- sqrt(0.01)
            add.sd <- sqrt(0.25)
            
      })
      model({
            V<-tV*exp(eta.v)
            KA<-tKa*exp(eta.ka)
            CL<-tCl*exp(eta.cl)
            Favail <- tF
            d/dt(depot) <- -KA*depot
            d/dt(central) <- KA*depot - (CL/V)*central
            f(depot) <- Favail*DOSE
            y <- central/V
            y ~ prop(prop.sd) + add(add.sd)
      })
}

# minxt, maxxt
e <- et(list(c(0, 25),
             c(0, 25),
             c(0, 25),
             c(0, 120),
             c(0, 120),
             c(0, 120),
             c(0, 120),
             c(0, 120))) %>%
      # users must be careful to not add a dosing record
      # et(amt=70, cmt="depot") %>%
      as.data.frame()

#xt
e$time <-  c(1,2,3,6,24,36,72,120)

babel.db.1 <- nlmixr2(f1, e, "poped",
                      popedControl(groupsize=32,
                                   bUseGrouped_xt=TRUE,
                                   a=c(DOSE=70),
                                   maxa=c(DOSE=200),
                                   mina=c(DOSE=0)))


babel.db.2 <- nlmixr2(f2, e, "poped",
                    popedControl(groupsize=32,
                                 bUseGrouped_xt=TRUE,
                                 a=c(DOSE=70),
                                 maxa=c(DOSE=200),
                                 mina=c(DOSE=0)))


# Add a dosing event as placeholder
e2 <- et(list(c(0, 25),
             c(0, 25),
             c(0, 25),
             c(0, 120),
             c(0, 120),
             c(0, 120),
             c(0, 120),
             c(0, 120))) %>%
      et(amt=1, cmt="depot") %>%
      as.data.frame()

babel.db.3 <- nlmixr2(f3, e2, "poped",
                    popedControl(groupsize=32,
                                 bUseGrouped_xt=TRUE,
                                 a=c(DOSE=70),
                                 maxa=c(DOSE=200),
                                 mina=c(DOSE=0)))


##  create plot of models
(plot1 <- plot_model_prediction(babel.db.1,IPRED=T))
(plot2 <- plot_model_prediction(babel.db.2,IPRED=T))
(plot3 <- plot_model_prediction(babel.db.3,IPRED=T))

# different results for different parameterizations
library(ggplot2)
library(gridExtra)
plot1 <- plot1 + ggtitle("CL Parameterization")
plot2 <- plot2 + ggtitle("KE Parameterization")
grid.arrange(plot1,plot2)


## evaluate initial designs
# different results for different parameterizations
evaluate_design(babel.db.1)
evaluate_design(babel.db.2)
evaluate_design(babel.db.3)
  

# Optimization of sample times
# different results for different parameterizations
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output.1 <- poped_optim(poped.db.1,opt_xt=T,parallel=T)
output.2 <- poped_optim(poped.db.2,opt_xt=T,parallel=T)

summary(output.1)
summary(output.2)
