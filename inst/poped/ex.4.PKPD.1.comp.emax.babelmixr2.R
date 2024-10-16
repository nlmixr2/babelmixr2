library(PopED)
library(babelmixr2)
library(tidyverse)

# Define the model
f <- function() {
      
      ini({
            # Initial estimates for population parameters (fixed effects)
            tCL <- 0.5  # Clearance
            tV <- 0.2   # Volume of distribution
            tE0 <- 1    # Baseline effect
            tEMAX <- 1  # Maximum effect
            tEC50 <- 1  # Concentration at 50% of maximum effect
            
            # Variability (random effects)
            eta.CL ~ 0.09
            eta.V ~ 0.09
            eta.E0 ~ 0.04
            eta.EC50 ~ 0.09
            
            # Residual unexplained variability (RUV)
            eps.pk.prop <- sqrt(0.15)   # PK proportional error for(SD/mean)
            eps.pdadd  <- sqrt(0.015)  # PD additive error
            
      })
      
      model({
            
            # PK model
            # Using a simple analytical funtion. Notice the use of the variable "DOSE". This will be defined later in the poped controler
            V <- tV * exp(eta.V)
            CL <- tCL * exp(eta.CL)
            cp = DOSE / V * exp(-CL / V * time)
            
            cp ~ prop(eps.pk.prop)
            
            # PD model
            E0   <- tE0 * exp(eta.E0)
            EMAX <- tEMAX
            EC50 <- tEC50 * exp(eta.EC50)
            
            effect = E0 + cp * EMAX / (EC50 + cp)
            effect ~ add(eps.pdadd)
            
      })
}

# e <-  et(c(1,3,8)) %>%
#       as.data.frame()

event.table.optim.pk <- 
      # amt is a placeholder here
      et(amt = 1) %>% 
      # specify observation records with appropriate bounds - these are observation windows
      et(replicate(4, c(0, 5), simplify = FALSE)) %>% 
      # make sure that the appropriate times are kept
      mutate(time = rep(c(0, 0.33,0.66,0.9,5), 1)) %>% 
      mutate(dvid = ifelse(time == 0, "dose", "cp"))

event.table.optim.pd <- 
      # specify observation records with appropriate bounds - these are observation windows
      et(replicate(4, c(0, 5), simplify = FALSE)) %>% 
      # make sure that the appropriate times are kept
      mutate(time = rep(c(0.1, 1, 2, 5), 1), 
             dvid = 'effect')

e <- event.table.optim.pk %>% bind_rows(event.table.optim.pd) %>% arrange(id, time)
# event.table.optim2 <- event.table.optim %>% filter(id ==2)
str(e)


## -- Define initial design  and design space
babel.db <- nlmixr2(f, e, "poped",
                    control=popedControl(
                          groupsize=20, 
                          m = 3,
                          bUseGrouped_xt = T,
                          minxt=0,
                          maxxt=5, 
                          ourzero = 0,
                          a=list(c(DOSE=0),c(DOSE=1),c(DOSE=2)), 
                          maxa = c(DOSE=10),
                          mina = c(DOSE=0)) )

plot_model_prediction(babel.db,facet_scales="free")
plot_model_prediction(babel.db,IPRED=T,DV=T,facet_scales="free",separate.groups=T)

## evaluate initial design
evaluate_design(babel.db)
shrinkage(babel.db)

# Optimization of sample times and doses
output <- poped_optim(babel.db, opt_xt = T, opt_a = T, parallel = T,method = c("LS"))

get_rse(output$FIM,output$babel.db)
plot_model_prediction(output$babel.db,facet_scales="free")

