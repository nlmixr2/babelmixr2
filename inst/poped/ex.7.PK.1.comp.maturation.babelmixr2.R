library(babelmixr2)

library(PopED)

f <- function() {
  ini({
    tCl <- 3.8
    tV <- 20
    tTM50 <- 60/52
    tHill <- 3
    eta.cl ~ 0.05
    eta.v ~ 0.05
    eta.tm50 ~ 0.05
    prop.sd <- sqrt(0.015)
    add.sd <- fix(sqrt(0.0015))
  })
  model({

    CL=tCl*exp(eta.cl)
    V= tV*exp(eta.v)
    TM50= tTM50*exp(eta.tm50)
    HILL= tHill

    # Maturation
    WTmax1 <- 2.76
    TM50wt1 <- 38.5
    HILLwt1 <- 12.9
    HILL2wt1 <- 2.74
    PMA1 <- PMA*52

    WTmax2 <- 16.4
    TM50wt2 <- 2.1
    HILLwt2 <- 2.04
    HILL2wt2 <- 1

    WTmax3 <- 40.2
    TM50wt3 <- 12.4
    HILLwt3 <- 2.87
    HILL2wt3 <- 0

    TLAGwt4 <- 12.4
    WTmax4 <- 33.6
    THALFwt4 <- 3.61

    FFEM <- 1
    if(SEX==1) FFEM <- 0.884

    wt1 <- FFEM*WTmax1/(1+(TM50wt1/PMA1)^((PMA1<TM50wt1)*HILLwt1+(PMA1>=TM50wt1)*HILL2wt1))
    wt2 <- WTmax2/(1+(TM50wt2/PMA)^((PMA<TM50wt2)*HILLwt2+(PMA>=TM50wt2)*HILL2wt2))
    wt3 <- WTmax3/(1+(TM50wt3/PMA)^((PMA<TM50wt3)*HILLwt3+(PMA>=TM50wt3)*HILL2wt3))
    wt4 <- (PMA > TLAGwt4)*(FFEM*WTmax4*(1-exp(-log(2)/THALFwt4*(PMA-TLAGwt4))))

    WT <- wt1 + wt2 + wt3 + wt4

    y=time

    CL=CL*(WT/70)^(3/4)*(PMA^HILL)/(TM50^HILL+PMA^HILL)
    V=V*(WT/70)
    DOSE=1000*(WT/70)
    y = DOSE/V*exp(-CL/V*time)

    y ~ add(add.sd) + prop(prop.sd)
  })
}

e <- et(c( 1,2,4,6,8,24)) %>% as.data.frame()

babel.db <- nlmixr2(f, e, "poped",
                    popedControl(
                      groupsize=rbind(50,20,20,20),
                      minxt=0,
                      maxxt=24,
                      bUseGrouped_xt=TRUE,
                      a=list(c(PMA=25,SEX=2),
                             c(PMA=15,SEX=2),
                             c(PMA=10,SEX=2),
                             c(PMA=5,SEX=2)),
                      maxa=c(PMA=30,SEX=2),
                      mina=c(PMA=1,SEX=2)))

plot_model_prediction(babel.db)

plot_model_prediction(babel.db,IPRED=T,DV=T,separate.groups=T)

## $rse
#       CL            V         TM50         HILL         d_CL          d_V       d_TM50
# 3.724808     2.251799  3048.507973  2117.356641    15.751133    14.997327 27968.366291
# SIGMA[1,1]
# 6.950869

## evaluate initial design
evaluate_design(babel.db)

shrinkage(babel.db)

# Optimization of sample times and WT
output <- poped_optim(babel.db,opt_xt=T,opt_a=T)

summary(output)
get_rse(output$FIM,output$poped.db)
plot_model_prediction(output$poped.db)
