library(babelmixr2)
library(PopED)
library(ggplot2)

# Does not work properly for now. PopED has a special way to handle IOV. This feature will be added in the future

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
            
            # IOV
            theta.iov.cl <- sqrt(c(0.09))
            
            iov.cl1 <-  fix(1)
            # iov.cl2 <- fix(1)
            # iov.cl3 <- fix (1) #  as many as occasions as needed
            
            prop.sd <- fix(sqrt(0.04))
            add.sd <- fix(sqrt(5e-6))
            
      })
      model({
            # No of doses
            N = floor(t/TAU) + 1
            OCC1 = ifelse(N > 6, 1, 0)
            iov <- theta.iov.cl * (OCC1 * iov.cl1) 
            
            V<-tV*exp(eta.v)
            KA<-tKa*exp(eta.ka)
            CL = tCl *exp(eta.cl + iov) 
            Favail <- tF

            d/dt(depot) <- -KA*depot
            d/dt(central) <- KA*depot - (CL/V)*central
            f(depot) <- Favail*DOSE
            y <- central/V
            y ~ prop(prop.sd) + add(add.sd)
      })
}

# careful with the et. the amt is scaled to the F
# minxt, maxxt
e <- et(list(c(0, 10),
             c(0, 10),
             c(0, 10),
             c(240, 248),
             c(240, 248))) %>%
      et(amt=1, ii=24, until=248,cmt="depot") %>%
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
plot_model_prediction(babel.db, model_num_points = 300) + coord_cartesian(ylim = c(0,0.5))

##  Visulaize the IOV 
set.seed(12345679)
plot_model_prediction(babel.db, PRED=F,IPRED=F, 
                      separate.groups=T, model_num_points = 500, 
                      groupsize_sim = 1,
                      IPRED.lines = T, alpha.IPRED.lines=0.6,
                      sample.times = F
) + geom_vline(xintercept = 24*6,color="red")

## evaluate initial design
evaluate_design(babel.db)

