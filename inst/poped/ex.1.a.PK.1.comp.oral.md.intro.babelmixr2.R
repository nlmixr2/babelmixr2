library(babelmixr2)

library(PopED)

##-- Model: One comp first order absorption
## -- Analytic solution for both mutiple and single dosing

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
    N <-  floor(time/TAU)+1
    y <- (DOSE*Favail/V)*(KA/(KA - CL/V)) *
      (exp(-CL/V * (time - (N - 1) * TAU)) *
         (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) -
         exp(-KA * (time - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))

    y ~ prop(prop.sd) + add(add.sd)
  })
}

e <- et(list(c(0, 10),
             c(0, 10),
             c(0, 10),
             c(240, 248),
             c(240, 248))) %>%
  as.data.frame()

# PopED xt equivalent
e$time <-  c(1,2,8,240,245)


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
## From original
## $rse
##        V        KA        CL       d_V      d_KA      d_CL
## 8.215338 10.090955  4.400304 39.833230 60.089601 27.391518

evaluate_design(babel.db)

## original: > shrinkage(poped.db)
## # A tibble: 9 × 5
##     d_V  d_KA   d_CL type       group
##   <dbl> <dbl>  <dbl> <chr>      <chr>
## 1 0.364 0.578 0.184  shrink_var all_groups
## 2 0.364 0.579 0.184  shrink_var grp_1
## 3 0.363 0.577 0.183  shrink_var grp_2
## 4 0.202 0.350 0.0965 shrink_sd  all_groups
## 5 0.202 0.351 0.0967 shrink_sd  grp_1
## 6 0.202 0.350 0.0963 shrink_sd  grp_2
## 7 0.181 0.228 0.107  se         all_groups
## 8 0.181 0.228 0.107  se         grp_1
## 9 0.181 0.228 0.107  se         grp_2
shrinkage(babel.db)

# Optimization of sample times
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output <- poped_optim(babel.db, opt_xt =TRUE, parallel=TRUE)

# Evaluate optimization results
summary(output)

## From original
# V        KA        CL       d_V      d_KA      d_CL 
# 6.281944  7.726279  4.295908 32.416232 49.062880 26.363021 
get_rse(output$FIM,output$poped.db)

plot_model_prediction(output$poped.db)

# Optimization of sample times and doses
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output_2 <- poped_optim(output$poped.db, opt_xt =TRUE, opt_a = TRUE, parallel = TRUE)

summary(output_2)

# From original
# V        KA        CL       d_V      d_KA      d_CL 
# 6.252332  7.547072  4.240929 32.205996 47.014629 25.684326 
get_rse(output_2$FIM,output_2$poped.db)

plot_model_prediction(output_2$poped.db)

# Optimization of sample times with only integer time points in design space
# faster than continuous optimization in this case
babel.db.discrete <- create.poped.database(babel.db,discrete_xt = list(0:248))

# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output_discrete <- poped_optim(babel.db.discrete, opt_xt=T, parallel = TRUE)

summary(output_discrete)

# V        KA        CL       d_V      d_KA      d_CL 
# 6.331614  8.009220  4.297905 32.351741 51.795028 26.386514 
get_rse(output_discrete$FIM,output_discrete$poped.db)

plot_model_prediction(output_discrete$poped.db)

# Efficiency of sampling windows
plot_efficiency_of_windows(output_discrete$poped.db, xt_windows=1)
