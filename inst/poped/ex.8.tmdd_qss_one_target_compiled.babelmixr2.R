
library(babelmixr2)

library(PopED)

f <- function() {
  ini({
    tCl     <- 0.3
    tVc     <- 3
    tQ      <- 0.2
    tVp     <- 3
    tFavail <- 0.7
    tKa     <- 0.5
    tVmax   <- fix(0)
    tKmss   <- fix(0)
    tR0     <- 0.1
    tKsss   <- 0.015
    tKdeg   <- 10
    tKint   <- 0.05

    eta.cl     ~ 0.09
    eta.vc     ~ 0.09
    eta.q      ~ 0.04
    eta.vp     ~ 0.04
    eta.favail ~ 0.04
    eta.ka     ~ 0.16
    eta.vmax   ~ fix(0)
    eta.kmss   ~ fix(0)
    eta.r0     ~ 0.09
    eta.ksss   ~ 0.09
    eta.kdeg   ~ 0.04
    eta.kint   ~ 0.04

    rtot.sd <- sqrt(0.04)
    free.sd <- sqrt(0.0225)
  })
  model({
    Cl     <- tCl     * exp(eta.cl)
    Vc     <- tVc     * exp(eta.vc)
    Q      <- tQ      * exp(eta.q)
    Vp     <- tVp     * exp(eta.vp)
    Favail <- tFavail * exp(eta.favail)
    Ka     <- tKa     * exp(eta.ka)
    Vmax   <- tVmax   * exp(eta.vmax)
    Kmss   <- tKmss   * exp(eta.kmss)
    R0     <- tR0     * exp(eta.r0)
    Ksss   <- tKsss   * exp(eta.ksss)
    Kdeg   <- tKdeg   * exp(eta.kdeg)
    Kint   <- tKint   * exp(eta.kint)

    Ctot          <- central/Vc
    Cfree         <- 0.5*((Ctot-Rtot-Ksss)+sqrt((Ctot-Rtot-Ksss)^2+4*Ksss*Ctot))
    depot(0)      <- DOSE*SC_FLAG
    d/dt(depot)   <- -Ka*depot
    central(0)    <- DOSE*(1-SC_FLAG)
    d/dt(central) <- Ka*depot*Favail + (Q/Vp)*periph - (Cl/Vc+Q/Vc)*Cfree*Vc -
      Rtot*Kint*Cfree*Vc/(Ksss+Cfree)
    d/dt(periph)  <- (Q/Vc)*Cfree*Vc - (Q/Vp)*periph
    Rtot(0)       <- R0
    d/dt(Rtot)    <- R0*Kdeg - Kdeg*Rtot - (Kint-Kdeg)*(Rtot*Cfree/(Ksss+Cfree))

    Rtot         ~  lnorm(rtot.sd)
    Cfree        ~  lnorm(free.sd)
  })
}

f <- f()

e1 <- et(c(0.0417,0.25,0.5,1,3,7,14,21,28,35,42,49,56)) %>%
  as.data.frame() %>% dplyr::mutate(dvid=1)

e2 <- e1
e2$dvid <- 2

e0 <- rbind(e1, e2) %>%
  dplyr::mutate(ID=1)

e1 <- et(c(0.0417,1,1,7,14,21,28,56,63,70,77,84,91,98,105)) %>%
  as.data.frame() %>% dplyr::mutate(dvid=1)

e2 <- e1
e2$dvid <- 2

e <- rbind(e1, e2) %>%
  dplyr::mutate(ID=2)

e <- rbind(e0, e)






#################################################
# for study 1 in gibiansky,JPKPD,2012 table 2 
#################################################
db.1 <- nlmixr2(f, e, "poped",
              control=popedControl(
                    groupsize=6,
                    a=list(c(ID=1, DOSE=100, SC_FLAG=0),
                           c(ID=1, DOSE=300, SC_FLAG=0),
                           c(ID=1, DOSE=600, SC_FLAG=0),
                           c(ID=1, DOSE=1000, SC_FLAG=1)),
                    discrete_a = list(DOSE=seq(100,1000,by=100),
                                      SC_FLAG=c(0,1)),
              ))

plot_model_prediction(db.1,facet_scales="free")

tic(); eval <- evaluate_design(db.1); toc()


# This is a longer running example, so the choice to check times for
# changes each time may not be ideal.

# For this first example, we will not check time differences, which
# speeds up solving a bit

db <- babel.poped.database(db, optTime=FALSE)

tic();e1 <- evaluate_design(db);toc()

# Original example:
# $rse
# CL         V1          Q         V2     FAVAIL         KA         R0       KSSS       KDEG       KINT       d_CL 
# 7.125675   6.867018   7.408125  10.435998  11.471400  19.592490   8.408247  11.028650   8.175673   7.343506  32.695140 
# d_V1        d_Q       d_V2   d_FAVAIL       d_KA       d_R0     d_KSSS     d_KDEG     d_KINT SIGMA[1,1] SIGMA[2,2] 
# 33.548134  74.864908  84.903024  98.434131  78.750729  36.359732  49.607416  66.478853  49.674898   8.935624   9.676856 

#################################################
# for study 1 + 2 in gibiansky,JPKPD,2012 table 2 
#################################################

db <- nlmixr2(f, e, "poped",
              control=popedControl(
                groupsize=rbind(6,6,6,6,100,100),
                a=list(c(ID=1, DOSE=100, SC_FLAG=0),
                       c(ID=1, DOSE=300, SC_FLAG=0),
                       c(ID=1, DOSE=600, SC_FLAG=0),
                       c(ID=1, DOSE=1000, SC_FLAG=1),
                       c(ID=2, DOSE=600, SC_FLAG=0),
                       c(ID=2, DOSE=1000, SC_FLAG=1)),
                discrete_a = list(DOSE=seq(100,1000,by=100),
                                  SC_FLAG=c(0,1)),
                ))

# This is a longer running example, so the choice to check times for
# changes each time may not be ideal.

# For this first example, we will not check time differences, which
# speeds up solving a bit

db <- babel.poped.database(db, optTime=FALSE)

tic();e1 <- evaluate_design(db);toc()

# You can see that without allowing times to be changed, you close the
# same value as the original ex.8.tmdd from PopED (with a little rounding error)

# Original example:
## $rse
##   CL         V1          Q         V2     FAVAIL         KA         R0       KSSS
## 2.418464   2.493863   2.390597   2.757785   3.023260   4.838732   2.800056   3.122798
## KDEG       KINT       d_CL       d_V1        d_Q       d_V2   d_FAVAIL       d_KA
## 2.711192   2.431927  10.830828  12.269644  22.271256  20.032780  24.240763  19.409710
##      d_R0     d_KSSS     d_KDEG     d_KINT SIGMA[1,1] SIGMA[2,2]
## 11.890735  13.446074  20.396815  18.125143   2.731176   2.906398

# If you check for time changes, this is slower and gives a slightly
# different answer
db <- babel.poped.database(db, optTime=TRUE)

tic();e2 <- evaluate_design(db);toc()

# now optimization in parallel for unix/mac
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output <- poped_optim(db,opt_xt = F, opt_a = T, parallel=T, method = c("LS")) 

# optimization for windows
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
# output <- poped_optim(poped.db.2,opt_xt = F, opt_a = T, parallel=T, method = c("LS"), dlls = c('tmdd_qss_one_target'))
plot_model_prediction(output$poped.db,facet_scales="free")

# You can see the differences between the two with microbenchmark.
m1 <- microbenchmark::microbenchmark(compareTime=evaluate_design(db), times=25)

# Models created by babelmixr2 check for time and model_switch
# changes.  If you are not optimizing over times, some can be lost
db <- babel.poped.database(db, optTime=FALSE)

m2 <- microbenchmark::microbenchmark(assumeTimesDontChange=evaluate_design(db), times=25)

# But this may not matter in the problem (like this one).  In most
# cases it is OK to leave the time comparison option enabled to allow
# optimizing over time

db <- babel.poped.database(db, optTime=TRUE)

library(ggplot2)

m3 <- rbind(m1, m2)

autoplot(m3)

# Going back to the model where time differences are not checked:
db <- babel.poped.database(db, optTime=FALSE)

plot_model_prediction(db, model_num_points=300, PI=TRUE, facet_scales="free")
