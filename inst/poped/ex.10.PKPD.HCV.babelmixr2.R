## HCV example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation 
##   for population pharmacokinetics-pharmacodynamics studies", 
##   Br. J. Clin. Pharm., 2014. 

library(babelmixr2)
library(PopED)

f <- function() {
      ini({
            
            tp     <- fix(100)
            td     <- fix(0.001)
            te     <- fix(1e-7)
            ts     <- fix(20000)
            
            tKA <- log(0.8)
            tKE <- log(0.15)
            tVD <- log(100) #VD=log(100000),
            tEC50 <- log(0.12) #EC50=log(0.00012),
            tn <- log(2)
            tdelta <- log(0.2)
            tc <- log(7)
            
            eta.KA ~ 0.25
            eta.KE ~ 0.25
            eta.VD ~ 0.25
            eta.EC50 ~0.25
            eta.n ~ 0.25
            eta.delta ~ 0.25
            eta.c ~ 0.25
            
            add.sd.pk <- sqrt(0.04) # nlmixr2 uses sd
            add.sd.pd <- sqrt(0.04)
      })
      model({
            
            p <- tp
            d <- td
            e <- te
            s <- ts
            KA <- exp(tKA + eta.KA)
            KE <- exp(tKE + eta.KE)
            VD <- exp(tVD + eta.VD)
            EC50 <- exp(tEC50 + eta.EC50)
            n <- exp(tn + eta.n)
            delta <- exp(tdelta + eta.delta)
            c <- exp(tc + eta.c)
            
            # # Apply the conditional logic -
            # # This does not work right now - could be addressed in the future
            # if ( (time - (floor(time / TAU) * TAU)) < TINF) {
            #       In <- DOSE / TINF  # Infusion rate
            # } else {
            #       In <- 0  # No infusion
            # }
            
            d/dt(depot)   = -KA*depot #  + In
            d/dt(central) =  KA*depot - KE*central
            d/dt(TC) = s - TC*(e*VP+d)    # target cells (TC)
            d/dt(IC) = e*TC*VP-delta*IC # productively infected cells (IC)
            d/dt(VP) = p*(1-(pow(central/VD,n)/(pow(central/VD,n)+pow(EC50,n))))*IC-c*VP # viral particles (VP)
            
            # Initial conditions
            dur(depot) = 1
            # depot(0) <- DOSE
            
            TC(0) =c*delta/(p*e)
            IC(0) =(s*e*p-d*c*delta)/(p*delta*e)
            VP(0) = (s*e*p-d*c*delta)/(c*delta*e)
            
            # Concentration and effect are calculated
            conc = central/VD
            eff  = log10(VP)
            ## And both are assumed to follow additive error
            conc ~ add(add.sd.pk)
            eff ~ add(add.sd.pd)
      })
}

rxClean()
f <- f()

TAU = 7
TINF = 1
e1 <- et(c( 0,0.25,0.5,1,2,3,4,7,10,14,21,28)) %>%
      add.dosing(dose=180, nbr.doses=4, dosing.interval=TAU, rate = 180/TINF) %>%
      as.data.frame() %>%
      dplyr::mutate(dvid=1)

e2 <- et(c( 0,0.25,0.5,1,2,3,4,7,10,14,21,28)) %>%
      add.dosing(dose=0, nbr.doses=4, dosing.interval=TAU, rate = 180/TINF) %>%
      as.data.frame() %>%
      dplyr::mutate(dvid=2) %>%
      dplyr::filter(is.na(ii))

e <- rbind(e1, e2)

babel.db <- nlmixr2(f, e, "poped",
                    popedControl(
                          groupsize=30
                          # ,a=list(c(DOSE=180,TINF=1,TAU=7))
                    ))

##  create plot of model without variability 
plot_model_prediction(babel.db, facet_scales = "free", model_num_points = 1000)
evaluate_design(babel.db)
shrinkage(babel.db)

#' #####################################
#' The reduced FIM 
#' ####################################
 # computation time
tic(); FIM_babel <- evaluate.fim(babel.db); toc()

#' design evaluation
crit <- det(FIM_babel)^(1/length(get_unfixed_params(babel.db)[["all"]]))
crit
rse <- get_rse(FIM_babel,babel.db) # this is for the log of the fixed effect parameters
rse_norm <- sqrt(diag(inv(FIM_babel)))*100 # this is approximately the RSE for the normal scale of the fixed effects 
rse[1:7] <- rse_norm[1:7] # replace the log scale for the normal scale RSE values
rse

#' Evaluation comparison to 
#' Nyberg et al., "Methods and software tools for design evaluation 
#' for population pharmacokinetics-pharmacodynamics studies", 
#' Br. J. Clin. Pharm., 2014. 
crit_reference_reduced <- 248.8
rse_reference_reduced <- c(12.1,10.5,10.0,15.8,10.4,9.4,11.0,40.0,30.8,28.8,60.4,28.8,27.2,32.8,8.5,9.0)

# the relative differences in percent
(rse - rse_reference_reduced)/rse_reference_reduced * 100
(crit - crit_reference_reduced)/crit_reference_reduced * 100

#' #####################################
#' The full FIM. 
#' #####################################
FIM_babel_full <- evaluate.fim(babel.db,fim.calc.type = 0) 

#' design evaluation
crit_full <- det(FIM_babel_full)^(1/length(get_unfixed_params(babel.db)[["all"]]))
crit_full
rse_full <- get_rse(FIM_babel_full,babel.db) # this is for the log of the fixed effect parameters
rse_norm_full <- sqrt(diag(inv(FIM_babel_full)))*100 # this is approximately the RSE for the normal scale of the fixed effects 
rse_full[1:7] <- rse_norm_full[1:7] # replace the log scale for the normal scale RSE values
rse_full

#' Evaluation compared to 
#' Nyberg et al., "Methods and software tools for design evaluation 
#' for population pharmacokinetics-pharmacodynamics studies", 
#' Br. J. Clin. Pharm., 2014. 
crit_reference_full <- 318.2
rse_reference_full <- c(8.6,6.9,8.4,13.5,7.5,8.5,8.7,43.2,37.2,33.2,66.4,32.8,31.6,33.6,8.5,9.3)

# the relative differences in percent
(rse_full - rse_reference_full)/rse_reference_full * 100
(crit_full - crit_reference_full)/crit_reference_full * 100



