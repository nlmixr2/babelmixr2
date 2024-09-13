
library(babelmixr2)

library(PopED)

##-- Model: One comp first order absorption + inhibitory imax
## -- works for both mutiple and single dosing
f <- function() {
  ini({
    tV <- 72.8
    tKa <- 0.25
    tCl <- 3.75
    tFavail <- fix(0.9)
    tE0 <- 1120
    tImax <- 0.87
    tIC50 <- 0.0993

    eta.v ~ 0.09
    eta.ka ~ 0.09
    eta.cl ~ 0.25^2
    eta.e0 ~ 0.09

    conc.prop.sd <- fix(sqrt(0.04))
    conc.sd <- fix(sqrt(5e-6))

    eff.prop.sd <- fix(sqrt(0.09))
    eff.sd <- fix(sqrt(100))
  })
  model({
    V<- tV*exp(eta.v)
    KA <- tKa*exp(eta.ka)
    CL <- tCl*exp(eta.cl)
    Favail <- tFavail
    E0 <- tE0*exp(eta.e0)
    IMAX <- tImax
    IC50 <- tIC50
    # PK
    N <- floor(time/TAU)+1
    CONC <- (DOSE*Favail/V)*(KA/(KA - CL/V)) *
      (exp(-CL/V * (time - (N - 1) * TAU)) *
         (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) -
         exp(-KA * (time - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))
    # PD model
    EFF <- E0*(1 - CONC*IMAX/(IC50 + CONC))

    CONC ~ add(conc.sd) + prop(conc.prop.sd)
    EFF ~ add(eff.sd) + prop(eff.prop.sd)

  })
}

# Note that design point 240 is repeated
e1 <- et(c( 1,2,8,240, 240)) %>%
  as.data.frame() %>%
  dplyr::mutate(dvid=1)

e1$low <- c(0,0,0,240, 240)
e1$high <- c(10,10,10,248, 248)
# Since the design point is repeated, there needs to be a grouping
# variable which is defined in the dataset as G_xt since it is defined
# in PopED as G_xt
e1$G_xt <- seq_along(e1$low)

e2 <- e1
e2$dvid <- 2
e <- rbind(e1, e2)

babel.db <- nlmixr2(f, e, "poped",
                    popedControl(
                      groupsize=20,
                      discrete_xt = list(0:248),
                      bUseGrouped_xt=TRUE,
                      a=list(c(DOSE=20,TAU=24),
                             c(DOSE=40, TAU=24),
                             c(DOSE=0, TAU=24)),
                      maxa=c(DOSE=200,TAU=40),
                      mina=c(DOSE=0,TAU=2),
                      ourzero=0
                    ))


##  create plot of model and design
plot_model_prediction(babel.db,facet_scales="free")

plot_model_prediction(babel.db,IPRED=T,DV=T,facet_scales="free",separate.groups=T)

## evaluate initial design
evaluate_design(babel.db)
shrinkage(babel.db)

# Optimization
output <- poped_optim(babel.db, opt_xt = T, parallel = T)

summary(output)

get_rse(output$FIM,output$poped.db)

plot_model_prediction(output$poped.db,facet_scales="free")
