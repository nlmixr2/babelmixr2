## Warfarin example from software comparison in:
## Nyberg et al., "Methods and software tools for design evaluation
##   for population pharmacokinetics-pharmacodynamics studies",
##   Br. J. Clin. Pharm., 2014.

## Evaluating with uncertainty around parameter values in the model

library(PopED)

library(babelmixr2)

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
# First define standard controler from nlmixr2
e <-  et(c(0.5, 1,2,6,24,36,72,120)) %>%
  as.data.frame()

babel.db <- nlmixr2(f, e, "poped",
                    popedControl(groupsize=32,
                                 minxt=0,
                                 maxxt=120,
                                 a=70,
                                 mina=0,
                                 maxa=100))


# Adding 10% Uncertainty to all fixed effects (not Favail)
bpop_vals_ed <- babel.db$parameters$bpop
for (n in row.names(bpop_vals_ed)) {
  if (n %in% c("tCl", "tV", "tKA")) {
    bpop_vals_ed[n,] <- c(4, # log-normal distribution;
                          # note 1: normal distribution
                          bpop_vals_ed[n,2], # original value
                          (bpop_vals_ed[n,2]*0.1)^2 # 10% of original value
                          )
  }
}

# Now update the database to include these new values

babel.db <- babel.poped.database(babel.db, bpop=bpop_vals_ed,
                                 ED_samp_size=20)



## -- Define initial design  and design space
## ElnD: E(ln(det(FIM))) evaluate.
## result is inaccurate (run several times to see)
## increase ED_samp_size for a more accurate calculation
tic();evaluate_design(babel.db,d_switch=FALSE,ED_samp_size=20); toc()

## Note that this is a simulation procedure so each time you run it
## you may get different results

## Also note that babelmixr2 values are similar

## $rse
## CL          V         KA       d_CL        d_V       d_KA SIGMA[1,1] SIGMA[2,2]
## 4.991010   2.977982  14.014207  29.802546  36.711408  26.754059  31.477157  25.297312

## optimization with line search search
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output_ls <- poped_optim(babel.db, opt_xt=T, parallel=T, method = "LS", d_switch=F, ED_samp_size=20)

## laplace does not seem to work with babelmixr2 or example:
## ED: E(det(FIM)) using Laplace approximation
## deterministic calculation, relatively fast
## can be more stable for optimization
tic(); evaluate_design(babel.db,d_switch=FALSE,use_laplace=TRUE); toc()

## optimization with Laplace
# Note: The parallel option does not work well with Windows machines at this moment. 
# Please set parallel = FALSE if you are working on a Windows machine
output_ls <- poped_optim(babel.db, opt_xt=T, parallel=T, method = "LS", d_switch=F, use_laplace=T)
