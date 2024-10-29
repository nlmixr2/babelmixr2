library(PopED)
library(babelmixr2)

## following example 1.4.2 in PFIM
f <- function() {
      ini({
            tKA <- 2
            tK  <- 0.25
            tV <- 15

            eta.ka ~ 1
            eta.k ~ 0.1
            eta.v ~ 0.25
            
            # Residual unexplained variability (RUV)
            add.sd  <- fix(0.5)   # Additive error
            prop.sd <- fix(0.15)  # Proportional error
            
      })
      model({
            
            KA <-  tKA * exp(eta.ka)
            K  <-  tK  * exp(eta.k)
            V  <-  tV  * exp(eta.v)

            y = (DOSE/V*KA/(KA-K)*(exp(-K*time)-exp(-KA*time)))
            y ~ prop(prop.sd) + add(add.sd)
            
      })
}

e <-  et(c(1,3,8)) %>%
      as.data.frame()

## -- Define initial design  and design space
babel.db <- nlmixr2(f, e, "poped",
                    control=popedControl(
                          groupsize=1,
                          minxt=0,
                          maxxt=10,
                          a=100))


(shr_out <- shrinkage(babel.db))
evaluate_design(babel.db)
plot_model_prediction(babel.db)


