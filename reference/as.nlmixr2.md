# Convert an object to a nlmixr2 fit object

Convert an object to a nlmixr2 fit object

## Usage

``` r
as.nlmixr2(
  x,
  ...,
  table = nlmixr2est::tableControl(),
  rxControl = rxode2::rxControl(),
  ci = 0.95
)

as.nlmixr(
  x,
  ...,
  table = nlmixr2est::tableControl(),
  rxControl = rxode2::rxControl(),
  ci = 0.95
)
```

## Arguments

- x:

  Object to convert

- ...:

  Other arguments

- table:

  is the
  [`nlmixr2est::tableControl()`](https://nlmixr2.github.io/nlmixr2est/reference/tableControl.html)
  options

- rxControl:

  is the
  [`rxode2::rxControl()`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)
  options, which is generally needed for how `addl` doses are handled in
  the translation

- ci:

  is the confidence interval of the residual differences calculated (by
  default 0.95)

## Value

nlmixr2 fit object

## Author

Matthew L. Fidler

## Examples

``` r

# \donttest{

# First read in the model (but without residuals)
mod <- nonmem2rx(system.file("mods/cpt/runODE032.ctl", package="nonmem2rx"),
                 determineError=FALSE, lst=".res", save=FALSE)
#> ℹ getting information from  '/home/runner/work/_temp/Library/nonmem2rx/mods/cpt/runODE032.ctl'
#> ℹ reading in xml file
#> ℹ done
#> ℹ reading in ext file
#> ℹ done
#> ℹ reading in phi file
#> ℹ done
#> ℹ reading in lst file
#> ℹ abbreviated list parsing
#> ℹ done
#> ℹ reading in grd file
#> ℹ done
#> ℹ splitting control stream by records
#> ℹ done
#> ℹ Processing record $INPUT
#> ℹ Processing record $MODEL
#> ℹ Processing record $gTHETA
#> ℹ Processing record $OMEGA
#> ℹ Processing record $SIGMA
#> ℹ Processing record $PROBLEM
#> ℹ Processing record $DATA
#> ℹ Processing record $SUBROUTINES
#> ℹ Processing record $PK
#> ℹ Processing record $DES
#> ℹ Processing record $ERROR
#> ℹ Processing record $ESTIMATION
#> ℹ Ignore record $ESTIMATION
#> ℹ Processing record $COVARIANCE
#> ℹ Ignore record $COVARIANCE
#> ℹ Processing record $TABLE
#> ℹ change initial estimate of `theta1` to `1.37034036528946`
#> ℹ change initial estimate of `theta2` to `4.19814911033061`
#> ℹ change initial estimate of `theta3` to `1.38003493562413`
#> ℹ change initial estimate of `theta4` to `3.87657341967489`
#> ℹ change initial estimate of `theta5` to `0.196446108190896`
#> ℹ change initial estimate of `eta1` to `0.101251418415006`
#> ℹ change initial estimate of `eta2` to `0.0993872449483344`
#> ℹ change initial estimate of `eta3` to `0.101302674763154`
#> ℹ change initial estimate of `eta4` to `0.0730497519364148`
#> ℹ read in nonmem input data (for model validation): /home/runner/work/_temp/Library/nonmem2rx/mods/cpt/Bolus_2CPT.csv
#> ℹ ignoring lines that begin with a letter (IGNORE=@)'
#> ℹ applying names specified by $INPUT
#> ℹ subsetting accept/ignore filters code: .data[-which((.data$SD == 0)),]
#> ℹ renaming 'ytype' to 'nmytype'
#> ℹ done
#>  
#>  
#> ℹ read in nonmem IPRED data (for model validation): /home/runner/work/_temp/Library/nonmem2rx/mods/cpt/runODE032.csv
#> ℹ done
#> ℹ changing most variables to lower case
#> ℹ done
#> ℹ replace theta names
#> ℹ done
#> ℹ replace eta names
#> ℹ done (no labels)
#> ℹ renaming compartments
#> ℹ done
#>  
#>  
#> ℹ solving ipred problem
#> ℹ done
#> ℹ solving pred problem
#> ℹ done

# define the model with residuals (and change the name of the
# parameters) In this step you need to be careful to not change the
# estimates and make sure the residual estimates are correct (could
# have to change var to sd).

 mod2 <-function() {
   ini({
     lcl <- 1.37034036528946
     lvc <- 4.19814911033061
     lq <- 1.38003493562413
     lvp <- 3.87657341967489
     RSV <- c(0, 0.196446108190896, 1)
     eta.cl ~ 0.101251418415006
     eta.v ~ 0.0993872449483344
     eta.q ~ 0.101302674763154
     eta.v2 ~ 0.0730497519364148
   })
   model({
     cmt(CENTRAL)
     cmt(PERI)
     cl <- exp(lcl + eta.cl)
     v <- exp(lvc + eta.v)
     q <- exp(lq + eta.q)
     v2 <- exp(lvp + eta.v2)
     v1 <- v
     scale1 <- v
     k21 <- q/v2
     k12 <- q/v
     d/dt(CENTRAL) <- k21 * PERI - k12 * CENTRAL - cl * CENTRAL/v1
     d/dt(PERI) <- -k21 * PERI + k12 * CENTRAL
     f <- CENTRAL/scale1
     f ~ prop(RSV)
   })
 }

# now we create another nonmem2rx object that validates the model above:

new <- as.nonmem2rx(mod2, mod)
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> ℹ copy 'dfSub' to nonmem2rx model
#> ℹ copy 'thetaMat' to nonmem2rx model
#> ℹ copy 'dfObs' to nonmem2rx model
#>  
#>  
#> ℹ solving ipred problem
#> ℹ done
#> ℹ solving pred problem
#> ℹ done

# once that is done, you can translate to a full nlmixr2 fit (if you wish)

fit <- as.nlmixr2(new)
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → finding duplicate expressions in EBE model...
#> → optimizing duplicate expressions in EBE model...
#> → compiling EBE model...
#>  
#>  
#> ✔ done
#> rxode2 5.1.4 using 2 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`
#> → Calculating residuals/tables
#> ✔ done

print(fit)
#> ── nlmixr² nonmem2rx reading NONMEM ver 7.4.3 ──
#> 
#>               OBJF      AIC      BIC Log-likelihood Condition#(Cov)
#> nonmem2rx 15977.28 20185.64 20237.23      -10083.82        335.4129
#>           Condition#(Cor)
#> nonmem2rx        2.096559
#> 
#> ── Time (sec $time): ──
#> 
#>             setup postprocess table compress NONMEM as.nlmixr2
#> elapsed 0.6365934       0.016 0.336    0.002 100.95      1.557
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>       Est.       SE   %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> lcl  1.370  0.02979  2.174    3.937 (3.713, 4.173)    32.64      1.935 
#> lvc  4.198  0.02952 0.7032    66.56 (62.82, 70.53)    32.33      2.465 
#> lq   1.380  0.05471  3.965    3.975 (3.571, 4.425)    32.65      40.50 
#> lvp  3.877  0.03483 0.8986    48.26 (45.07, 51.67)    27.53      28.37 
#> RSV 0.1964 0.003153  1.605 0.1964 (0.1903, 0.2026)                     
#>  
#>   Covariance Type ($covMethod): nonmem2rx
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     
#> 
#>  WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
#> 
#>  (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
#> 
#>     
#> 0MINIMIZATION SUCCESSFUL
#>  NO. OF FUNCTION EVALUATIONS USED:      320
#>  NO. OF SIG. DIGITS IN FINAL EST.:  2.5
#> 
#>     IPRED relative difference compared to Nonmem IPRED: 0%; 95% percentile: (0%,0%); rtol=6.43e-06
#>     PRED relative difference compared to Nonmem PRED: 0%; 95% percentile: (0%,0%); rtol=6.41e-06
#>     IPRED absolute difference compared to Nonmem IPRED: 95% percentile: (2.25e-05, 0.0418); atol=0.00167
#>     PRED absolute difference compared to Nonmem PRED: 95% percentile: (1.41e-07,0.00382); atol=6.41e-06
#>     nonmem2rx model file: '/home/runner/work/_temp/Library/nonmem2rx/mods/cpt/runODE032.ctl' 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 2,280 × 25
#>   ID     TIME    DV  PRED    RES IPRED  IRES  IWRES eta.cl eta.v  eta.q eta.v2
#>   <fct> <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl>  <dbl>  <dbl> <dbl>  <dbl>  <dbl>
#> 1 1      0.25 1041. 1750. -710.  1215. -175. -0.732 -0.144 0.375 0.0650  0.241
#> 2 1      0.5  1629  1700.  -70.8 1192.  437.  1.87  -0.144 0.375 0.0650  0.241
#> 3 1      0.75  878. 1651. -774.  1169. -291. -1.27  -0.144 0.375 0.0650  0.241
#> # ℹ 2,277 more rows
#> # ℹ 13 more variables: f <dbl>, CENTRAL <dbl>, PERI <dbl>, cl <dbl>, v <dbl>,
#> #   q <dbl>, v2 <dbl>, v1 <dbl>, scale1 <dbl>, k21 <dbl>, k12 <dbl>, tad <dbl>,
#> #   dosenum <dbl>

# }
```
