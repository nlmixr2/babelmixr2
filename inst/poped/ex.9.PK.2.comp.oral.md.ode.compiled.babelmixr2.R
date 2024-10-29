#' An implementation of a two compartment model with oral absorption using ODEs.
library(PopED)

library(babelmixr2)

f <- function() {
  ini({
    #c(CL=10,V1=100,KA=1,Q= 3.0, V2= 40.0, Favail=1),
    tCl <- 10
    tV1 <- 100
    tKa <- 1
    tQ <- 3
    tV2 <- 40
    tFavail <- fix(1)
    eta.cl ~ 0.15^2
    eta.ka ~ 0.25^2
    prop.sd <- 0.1
    add.sd <- 0.05
  })
  model({
    CL <- tCl*exp(eta.cl)
    V1 <- tV1
    KA <- tKa*exp(eta.ka)
    Q <- tQ
    V2 <- tV2
    Favail <- tFavail
    f(depot) <- DOSE
    d/dt(depot) <- -KA*depot
    d/dt(central) <- KA*depot + periph * Q/V2 - central*(CL/V1+Q/V1)
    d/dt(periph) <- central* Q/V1 - periph* Q/V2
    cp <- central/(V1/Favail)
    cp ~ add(add.sd) + prop(prop.sd)
  })
}

# not quite the same since tau can be optimized over in original
# example.
e <- et(c( 48,50,55,65,70,85,90,120)) %>%
  et(amt=1, ii=24, until=120)

babel.db <- nlmixr2(f, e, "poped",
                    popedControl(groupsize=20,
                                 minxt=0,
                                 maxxt=144,
                                 discrete_xt = list(0:144),
                                 a=c(DOSE=100),
                                 maxa=c(DOSE=1000),
                                 mina=c(DOSE=0),
                                 discrete_a = list(DOSE=seq(0,1000,by=100))))


#' plot intial design just PRED
plot_model_prediction(babel.db,model_num_points = 500)

#' plot intial design with BSV and RUV in model
plot_model_prediction(babel.db,IPRED=T,DV=T)

#' how long does one evaluation of the FIM take?
tic(); (eval <- evaluate_design(babel.db)); toc()

##  create plot of model without variability
plot_model_prediction(babel.db,model_num_points = 500)

##  create plot of model with variability
plot_model_prediction(babel.db,IPRED=T,DV=T,model_num_points = 500)

#' how long does one evaluation of the FIM take?
tic(); (eval_compiled <- evaluate_design(babel.db)); toc()

#' very small differences in computation value
#' but a large difference in computation time (8 times faster with the compiled code)
(eval_compiled$ofv-eval$ofv)/eval$ofv

#' making optimization times more reasonable
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
#output <- poped_optim(babel.db,opt_xt=T, opt_a=T, parallel=T)
