library(babelmixr2)
library(PopED)


## define the ODE
f <- function() {
  ini({
    tV <- 72.8
    tKa <- 0.25
    tCl <- 3.75
    tF <- fix(0.9)

    eta.v ~ 0.09
    eta.ka ~ 0.09
    eta.cl ~0.25^2

    prop.sd <- fix(sqrt(0.04))
    add.sd <- fix(sqrt(5e-6))

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
e <- et(list(c(0, 10),
             c(0, 10),
             c(0, 10),
             c(240, 248),
             c(240, 248))) %>%
  et(amt=1/0.9, ii=24, until=248,cmt="depot") %>%
  as.data.frame()

#xt
e$time <-  c(0, 1,2,8,240,245)


babel.db <- nlmixr2(f, e, "poped",
                    popedControl(groupsize=20,
                                 bUseGrouped_xt=TRUE,
                                 a=list(c(DOSE=20,TAU=24),
                                        c(DOSE=40, TAU=24)),
                                 maxa=c(DOSE=200,TAU=24),
                                 mina=c(DOSE=0,TAU=24)))

##  create plot of model without variability
plot_model_prediction(babel.db, model_num_points = 300)

##  create plot of model with variability
plot_model_prediction(babel.db, IPRED=T, DV=T, separate.groups=T, model_num_points = 300)

## evaluate initial design
evaluate_design(babel.db)

shrinkage(babel.db)

# Optimization of sample times
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output <- poped_optim(babel.db, opt_xt =TRUE, parallel=TRUE, method = c("LS"))

# Evaluate optimization results
summary(output)

get_rse(output$FIM,output$poped.db)

plot_model_prediction(output$poped.db)


# Optimization of sample times and doses
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output_2 <- poped_optim(output$poped.db, opt_xt =TRUE, opt_a = TRUE, parallel=TRUE, method = c("LS"))

summary(output_2)

get_rse(output_2$FIM,output_2$poped.db)

plot_model_prediction(output_2$poped.db)


# Optimization of sample times with only integer time points in design space
# faster than continuous optimization in this case
babel.db.discrete <- create.poped.database(babel.db,discrete_xt = list(0:248))

# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output_discrete <- poped_optim(babel.db.discrete, opt_xt=T, parallel=TRUE, method = c("LS"))


summary(output_discrete)

get_rse(output_discrete$FIM,output_discrete$poped.db)

plot_model_prediction(output_discrete$poped.db)


# Efficiency of sampling windows
plot_efficiency_of_windows(output_discrete$poped.db, xt_windows=1)
