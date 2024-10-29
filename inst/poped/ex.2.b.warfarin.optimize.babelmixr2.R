## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation
##   for population pharmacokinetics-pharmacodynamics studies",
##   Br. J. Clin. Pharm., 2014.

## Optimization using an additive + proportional reidual error to
##   avoid sample times at very low concentrations (time 0 or very late samoples).
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
plot_model_prediction(babel.db,PI=TRUE)

#########################################
## NOTE All PopED output for residuals
## (add or prop) are VARIANCES instead of
## standard deviations!
#########################################

## evaluate initial design
evaluate_design(babel.db)

##############
# Optimization
##############

# below are a number of ways to optimize the problem

# Optimization of sample times
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output <- poped_optim(babel.db, opt_xt = T, parallel = T)
get_rse(output$FIM,output$poped.db)
plot_model_prediction(output$poped.db)

# Examine efficiency of sampling windows
plot_efficiency_of_windows(output$poped.db,xt_windows=0.5)


# Optimization of DOSE and sampling times
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output_D_T <- poped_optim(babel.db, opt_xt = T, opt_a = T)

summary(output_D_T)

plot_model_prediction(output_D_T$poped.db)

# Discrete optimization with only integer times allowed
# and Dose in units of 10

# This requires the original design space to only include integers, so
# this needs to be updated:


## -- Define discrete design space
e <-  et(c(1,2,6,24,36,72,120)) %>%
  as.data.frame()

# Also as a note, since babelmixr2 speeds up solving by preloading
# design points, the PopED design database from babelmixr2 should be
# updated to avoid any errors/issues (instead of creating using
# PopED::create.poped.database).  If the original design points are
# the same, this step does not need to be performed

babel.db.discrete <- nlmixr2(f, e, "poped",
                             popedControl(groupsize=32,
                                          minxt=0,
                                          maxxt=120,
                                          a=70,
                                          mina=0,
                                          maxa=100,
                                          discrete_xt=list(0:120),
                                          discrete_a=list(seq(10, 100, 10))))

# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output_discrete <- poped_optim(babel.db.discrete, opt_xt = T, opt_a = T, parallel = T)

get_rse(output_discrete$FIM,output_discrete$poped.db)

plot_model_prediction(output_discrete$poped.db)


# Optimization using a genetic algorithm
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output_ga <- poped_optim(babel.db, opt_xt = T, parallel = T, method = "GA")

summary(output_ga)

plot_model_prediction(output_ga$poped.db)
