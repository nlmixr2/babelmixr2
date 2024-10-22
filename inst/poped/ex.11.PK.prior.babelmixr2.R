library(babelmixr2)

library(PopED)
# This example shows how to include a prior FIM into the design evaluation.
# We look at PK assessment in pediatrics, where we are mainly interested
# in assessing if there is a substantial difference of more than 20% in
# clearance (CL) between children and adults.

# First we setup the general model, then the designs for adults and pediatrics,
# and finally we can evaluate the separate designs and the pooled data.

##-- Model: One comp first order absorption
## -- Analytic solution for both mutiple and single dosing

# Define the model
f <- function() {
  ini({
    tV <- 72.8
    tKa <- 0.25
    tCl <- 3.75
    tF <- fix(0.9)
    pedCL <- 0.8

    eta.v ~ 0.09
    eta.ka ~ 0.09
    eta.cl ~0.25^2

    prop.sd <- fix(sqrt(0.04))
    add.sd <- fix(sqrt(5e-6))

  })
  model({
    V<-tV*exp(eta.v)
    KA<-tKa*exp(eta.ka)
    CL<-tCl*exp(eta.cl)* (pedCL^(isPediatric)) # add covariate for pediatrics
    Favail <- tF

    N <-  floor(t/TAU)+1
    y <- (DOSE*Favail/V)*(KA/(KA - CL/V)) *
      (exp(-CL/V * (t - (N - 1) * TAU)) *
         (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) -
         exp(-KA * (t - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))

    y ~ prop(prop.sd) + add(add.sd)
  })
}

e <- et(c( 1,8,10,240,245))

babel.db <- nlmixr2(f, e, "poped",
                    popedControl(m = 2,
                                 groupsize=20,
                                 bUseGrouped_xt=TRUE,
                                 a=list(c(DOSE=20,TAU=24,isPediatric = 0),
                                        c(DOSE=40, TAU=24,isPediatric = 0))))


# Note, to be able to use the adults FIM to combine with the pediatrics,
# both have to have the parameter "pedCL" defined and set notfixed_bpop to 1.

## Define pediatric model/design (isPediatric = 1)
## One arm, 4 time points only

e.ped <- et(c( 1,2,6,240))

babel.db.ped <- nlmixr2(f, e.ped, "poped",
                        popedControl(m = 1,
                                     groupsize=6,
                                     bUseGrouped_xt=TRUE,
                                     a=list(c(DOSE=40,TAU=24,isPediatric = 1))))


##  Create plot of model of adult data without variability
plot_model_prediction(babel.db, model_num_points = 300)

##  Create plot of model of pediatric (single dose level) data without variability
plot_model_prediction(babel.db.ped, model_num_points = 300)

## To store FIM from adult design - need FIM from evaluate_design

## $rse
## V        KA        CL     pedCL       d_V      d_KA      d_CL
## 6.634931  8.587203  4.354792        NA 33.243601 55.689432 27.133255

## Warning message:
##           The following parameters are not estimable:
##                                              pedCL
## Is the design adequate to estimate all parameters?

(outAdult <- evaluate_design(babel.db))
# It is obvious that we cannot estimate the pediatric covariate from adult
# data only - therefore the message from the calculation.
# You can also note the zeros in the 4th column and 4th row of the FIM.

# We can evaluate the adult design without warning, by setting the pedCLs
# parameter to be fixed (i.e., not estimated)c
f2 <- f %>% ini(pedCL = fix(0.8))

babel.db.2 <- nlmixr2(f2, e, "poped",
                      popedControl(m = 2,
                                   groupsize=20,
                                   bUseGrouped_xt=TRUE,
                                   a=list(c(DOSE=20,TAU=24,isPediatric = 0),
                                          c(DOSE=40, TAU=24,isPediatric = 0))))

## PopED
##
## $rse
##        V        KA        CL       d_V      d_KA      d_CL
## 6.634931  8.587203  4.354792 33.243601 55.689432 27.133255
evaluate_design(babel.db.2)
# One obtains good estimates for all parameters for adults (<60% RSE for all).

## evaluate design of pediatrics only - insufficient
# Similarly as before with only pediatrics we cannot estimate the covariate effect,
# so we fix it.
babel.db.ped.2 <- nlmixr2(f2, e.ped, "poped",
                          popedControl(m = 1,
                                       groupsize=6,
                                       bUseGrouped_xt=TRUE,
                                       a=list(c(DOSE=40,TAU=24,isPediatric = 1)),
                                       literalFix=FALSE))


## PopED
## $rse
##        V        KA        CL       d_V      d_KA      d_CL
## 24.72088  30.84953  11.94767 116.23095 181.19778  77.29188

evaluate_design(babel.db.ped.2)
# Due to having less subjects, less samples per subject, and only one dose level
# the variability in pediatrics cannot be estimated well.

## Add adult prior
# Now we combined the two studies, where we assume that we only assess
# a difference between adults and pediatrics in CL.  We can set the
# prior FIM to the adult one: Note to setup babelmixr2 models
# appropriately we need to use `babel.poped.database` instead of
# `create.poped.database`
babel.db.all <- babel.poped.database(
  babel.db.ped,
  prior_fim = outAdult$fim
)


## evaluate design using prior FIM from adults

## PopED
## $rse
##        V        KA        CL     pedCL       d_V      d_KA      d_CL
## 6.381388  8.222819  4.354761 12.591940 31.808871 52.858399 25.601551

(out.all <- evaluate_design(babel.db.all))
# Obviously, the pooled data leads to much higher precision in parameter estimates
# compared to the pediatrics only.

# One can also obtain the power for estimating the covariate to be different from 1.

## PopED
##
## $rse
##        V        KA        CL     pedCL       d_V      d_KA      d_CL
## 6.381388  8.222819  4.354761 12.591940 31.808871 52.858399 25.601551

## $power
##       Value      RSE power_pred power_want need_rse min_N_tot
## pedCL   0.8 12.59194   51.01851         80 8.923519        14
evaluate_power(babel.db.all,
               bpop_idx=babelBpopIdx(babel.db.all, "pedCL"),
               h0=1,out=out.all)
