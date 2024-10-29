## using libary models and reparameterizing the problen to KA, KE and V
## optimization of dose and dose interval

library(babelmixr2)

library(PopED)

f <- function() {
  ini({
    tV <- 72.8
    tKa <- 0.25
    tKe <- 3.75/72.8
    tFavail <- fix(0.9)

    eta.v ~ 0.09
    eta.ka ~ 0.09
    eta.ke ~ 0.25^2

    prop.sd <- fix(sqrt(0.04))
    add.sd <- fix(sqrt(5e-6))
  })
  model({
    V <- tV*exp(eta.v)
    KA <- tKa*exp(eta.ka)
    KE <- tKe*exp(eta.ke)
    Favail <- tFavail
    N <- floor(time/TAU)+1
    y <- (DOSE*Favail/V)*(KA/(KA - KE)) *
      (exp(-KE * (time - (N - 1) * TAU)) * (1 - exp(-N * KE * TAU))/(1 - exp(-KE * TAU)) -
         exp(-KA * (time - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))

    y ~ prop(prop.sd) + add(add.sd)
  })
}

# minxt, maxxt
e <- et(list(c(0, 10),
             c(0, 10),
             c(0, 10),
             c(240, 248),
             c(240, 248))) %>%
  as.data.frame()

#xt
e$time <-  c(1,2,8,240,245)


babel.db <- nlmixr2(f, e, "poped",
                    popedControl(groupsize=20,
                                 bUseGrouped_xt=TRUE,
                                 a=list(c(DOSE=20,TAU=24),
                                        c(DOSE=40, TAU=24)),
                                 maxa=c(DOSE=200,TAU=24),
                                 mina=c(DOSE=0,TAU=24)))


##  create plot of model without variability
plot_model_prediction(babel.db)

##  create plot of model with variability
plot_model_prediction(babel.db,IPRED=T,DV=T,separate.groups=T)

## evaluate initial design
evaluate_design(babel.db)

shrinkage(babel.db)

# Optimization of sample times
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output <- poped_optim(babel.db, opt_xt =TRUE, parallel=TRUE)

# Evaluate optimization results
summary(output)

# from original
# V        KA        KE       d_V      d_KA      d_KE 
# 6.303480  7.899813  5.750915 29.156362 43.794939 33.388714
get_rse(output$FIM,output$poped.db)

plot_model_prediction(output$poped.db)

# Optimization of sample times, doses and dose intervals
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output_2 <- poped_optim(output$poped.db, opt_xt =TRUE, opt_a = TRUE, parallel=TRUE)

summary(output_2)
get_rse(output_2$FIM,output_2$poped.db)
plot_model_prediction(output_2$poped.db)

# Optimization of sample times with only integer time points in design space
# faster than continuous optimization in this case
babel.db.discrete <- create.poped.database(babel.db,discrete_xt = list(0:248))

# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output_discrete <- poped_optim(babel.db.discrete, opt_xt=T, parallel=TRUE)

summary(output_discrete)

# from original
# V        KA        KE       d_V      d_KA      d_KE 
# 6.367801  8.230615  5.835965 29.067957 45.786136 33.418356 
get_rse(output_discrete$FIM,output_discrete$poped.db)

plot_model_prediction(output_discrete$poped.db)


# Efficiency of sampling windows
plot_efficiency_of_windows(output_discrete$poped.db, xt_windows=1)
