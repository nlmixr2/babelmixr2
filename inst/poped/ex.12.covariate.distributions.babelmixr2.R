library(PopED)
library(babelmixr2)
## Introduction
# Lets assume that we have a model with a covariate included in the model description.  
# Here we define a one-compartment PK model that has weight on both clearance and volume of distribution. 

f <- function() {
      ini({
            tCl <- 3.8
            tV <- 20
            WT_CL <- 0.75
            WT_V <- 1
            
            eta.cl ~ 0.05
            eta.v ~ 0.05

            # Residual unexplained variability (RUV)
            add.sd  <- fix(sqrt(0.0015))  # Additive error
            prop.sd <- sqrt(0.015)        # Proportional error
            
      })
      model({
            
            CL <- tCl * (WT/70)^(WT_CL) * exp(eta.cl)
            V <-  tV  * (WT/70)^(WT_V) * exp(eta.v)
            
            DOSE = 1000*(WT/70)

            y = DOSE/V*exp(-CL/V*time) 
            y ~ prop(prop.sd) + add(add.sd)
            
      })
}

e <-  et(c( 1,2,4,6,8,24)) %>%
      as.data.frame()

## -- Define initial design  and design space
babel.db <- nlmixr2(f, e, "poped",
                    control=popedControl(
                          groupsize=50,
                          minxt=0,
                          maxxt=24,
                          bUseGrouped_xt = T,
                          a=c(WT=70)))

# We can create a plot of the model typical predictions:
plot_model_prediction(babel.db)

# And evaluate the initial design
evaluate_design(babel.db)

# We see that the covariate parameters can not be estimated 
# according to this design calculation (RSE of bpop[3]=0 and bpop[4]=0).  
# Why is that? Well, the calculation being done is assuming that every 
# individual in the group has the same covariate (to speed 
# up the calculation).  This is clearly a poor prediction in this case!

# distribution of covariates: We can improve the computation by assuming a
#distribution of the covariate (WT) in the individuals in the study. We set
#`groupsize=1`, the number of groups to be 50 (`m=50`) and assume that WT is
#sampled from a normal distribution with mean=70 and sd=10
#(`a=as.list(rnorm(50,mean = 70,sd=10)`).
babel.db.2 <- nlmixr2(f, e, "poped",
                      control=popedControl(
                            groupsize=1,
                            m=50,
                            minxt=0,
                            maxxt=24,
                            bUseGrouped_xt = T,
                            a = lapply(rnorm(50, mean = 70, sd = 10), function(x) c(WT=x)) ))


evaluate_design(babel.db.2)

#Here we see that, given this distribution of weights, the covariate effect
#parameters (bpop[3] and bpop[4]) would be well estimated.

#However, we are only looking at one sample of 50 individuals.  Maybe a better
#approach is to look at the distribution of RSEs over a number of experiments
#given the expected weight distribution.
nsim <- 10
rse_list <- c()
for(i in 1:nsim){
      poped_db_tmp <- nlmixr2(f, e, "poped",
                            control=popedControl(
                                  groupsize=1,
                                  m=50,
                                  minxt=0,
                                  maxxt=24,
                                  bUseGrouped_xt = T,
                                  a = lapply(rnorm(50, mean = 70, sd = 10), function(x) c(WT=x)) ))
      
      rse_tmp <- evaluate_design(poped_db_tmp)$rse
      rse_list <- rbind(rse_list,rse_tmp)
}
apply(rse_list,2,quantile)

