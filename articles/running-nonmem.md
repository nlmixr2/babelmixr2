# Running NONMEM with nlmixr2

``` r

library(babelmixr2)
```

## Step 0: What do you need to do to have `nlmixr2` run `NONMEM` from a nlmixr2 model

To use `NONMEM` in nlmixr, you do not need to change your data or your
`nlmixr2` dataset. `babelmixr2` will do the heavy lifting here.

You do need to setup how to run `NONMEM`. For many cases this is easy;
You simply have to figure out the command to run `NONMEM` (it is often
useful to use the full command path). You can set it in
`options("babelmixr2.nonmem"="nmfe743")` or use
`nonmemControl(runCommand="nmfe743")`. I prefer the
[`options()`](https://rdrr.io/r/base/options.html) method since you only
need to set it once. This could also be a function if you prefer (but I
will not cover using the function here).

## Step 1: Run a `nlmixr2` in NONMEM

Lets take the classic warfarin example to start the comparison.

The model we use in the `nlmixr2` vignettes is:

``` r

library(babelmixr2)
pk.turnover.emax3 <- function() {
  ini({
    tktr <- log(1)
    tka <- log(1)
    tcl <- log(0.1)
    tv <- log(10)
    ##
    eta.ktr ~ 1
    eta.ka ~ 1
    eta.cl ~ 2
    eta.v ~ 1
    prop.err <- 0.1
    pkadd.err <- 0.1
    ##
    temax <- logit(0.8)
    tec50 <- log(0.5)
    tkout <- log(0.05)
    te0 <- log(100)
    ##
    eta.emax ~ .5
    eta.ec50  ~ .5
    eta.kout ~ .5
    eta.e0 ~ .5
    ##
    pdadd.err <- 10
  })
  model({
    ktr <- exp(tktr + eta.ktr)
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    emax = expit(temax+eta.emax)
    ec50 =  exp(tec50 + eta.ec50)
    kout = exp(tkout + eta.kout)
    e0 = exp(te0 + eta.e0)
    ##
    DCP = center/v
    PD=1-emax*DCP/(ec50+DCP)
    ##
    effect(0) = e0
    kin = e0*kout
    ##
    d/dt(depot) = -ktr * depot
    d/dt(gut) =  ktr * depot -ka * gut
    d/dt(center) =  ka * gut - cl / v * center
    d/dt(effect) = kin*PD -kout*effect
    ##
    cp = center / v
    cp ~ prop(prop.err) + add(pkadd.err)
    effect ~ add(pdadd.err) | pca
  })
}
```

Now you can run the `nlmixr2` model using `NONMEM` you simply can run it
directly:

``` r

try(nlmixr(pk.turnover.emax3, nlmixr2data::warfarin, "nonmem",
           nonmemControl(readRounding=FALSE, modelName="pk.turnover.emax3")),
    silent=TRUE)
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> 
#> 
#>  WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
#> 
#>  (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
#> 
#> 
#> 0MINIMIZATION TERMINATED
#>  DUE TO ROUNDING ERRORS (ERROR=134)
#>  NO. OF FUNCTION EVALUATIONS USED:     1088
#>  NO. OF SIG. DIGITS UNREPORTABLE
#> 0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
#> 
#> nonmem model: 'pk.turnover.emax3-nonmem/pk.turnover.emax3.nmctl'
#> → terminated with rounding errors, can force nlmixr2/rxode2 to read with nonmemControl(readRounding=TRUE)
```

That this is the same way you would run an ordinary `nlmixr2` model, but
it is simply a new estimation method `"nonmem"` with a new controller
([`nonmemControl()`](https://nlmixr2.github.io/babelmixr2/reference/nonmemControl.md))
to setup options for estimation.

A few options in the
[`nonmemControl()`](https://nlmixr2.github.io/babelmixr2/reference/nonmemControl.md)
here is `modelName` which helps control the output directory of `NONMEM`
(if not specified `babelmixr2` tries to guess based on the model name
based on the input).

If you try this yourself, you see that `NONMEM` fails with rounding
errors. You could do the standard approach of changing `sigdig`, `sigl`,
`tol` etc, to get a successful `NONMEM` model convergence, of course
that is supported. But with `babelmixr2` you can *do more*.

## Optional Step 2: Recover a failed NONMEM run

One of the other approaches is to **ignore** the rounding errors that
have occurred and read into `nlmixr2` anyway:

``` r

# Can still load the model to get information (possibly pipe) and create a new model
f <- nlmixr(pk.turnover.emax3, nlmixr2data::warfarin, "nonmem",
            nonmemControl(readRounding=TRUE, modelName="pk.turnover.emax3"))
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → finding duplicate expressions in EBE model...
#> → compiling EBE model...
#> ✔ done
#> rxode2 5.1.8 using 2 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`
#> → Calculating residuals/tables
#> ✔ done
#> → compress origData in nlmixr2 object, save 27800
#> → compress parHistData in nlmixr2 object, save 5048
```

You may see more work happening than you expected to need for an already
completed model. When reading in a NONMEM model, `babelmixr2` grabs:

- `NONMEM`’s objective function value
- `NONMEM`’s covariance (if available)
- `NONMEM`’s optimization history
- `NONMEM`’s final parameter estimates (including the ETAs)
- `NONMEM`’s `PRED` and `IPRED` values (for validation purposes)

These are used to solve the ODEs *as if they came from an nlmixr2*
optimization procedure.

This means that you can compare the `IPRED` and `PRED` values of
`nlmixr2`/`rxode2` and *know immediately* if your model validates.

This is similar to the procedure Kyle Baron advocates for validating a
NONMEM model against a `mrgsolve` model (see
<https://mrgsolve.org/blog/posts/2022-05-validate-translation/> and
<https://mrgsolve.org/blog/posts/2023-update-validation.html>),

The advantage of this method is that you need to simply write one model
to get a validated `roxde2`/`nlmixr2` model.

In this case you can see the validation when you print the fit object:

``` r

print(f)
#> ── nlmixr² nonmem ver 7.4.3 ──
#> 
#>                 OBJF      AIC      BIC Log-likelihood Condition#(Cov)
#> nonmem focei 1326.91 2252.605 2332.025      -1107.302              NA
#>              Condition#(Cor)
#> nonmem focei              NA
#> 
#> ── Time (sec $time): ──
#> 
#>             setup preprocess postprocess table compress NONMEM
#> elapsed 0.8129762      0.031       0.021 0.045    0.014 320.27
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>               Est. SE %RSE Back-transformed(95%CI) BSV(CV% or SD) Shrink(SD)%
#> tktr       6.24e-7                            1.00           86.5       59.8 
#> tka       -3.01e-6                            1.00           86.5       59.8 
#> tcl          -2.00                           0.135           28.6       1.34 
#> tv            2.05                            7.78           22.8       6.44 
#> prop.err    0.0986                          0.0986                           
#> pkadd.err    0.512                           0.512                           
#> temax         6.42                           0.998        0.00707        100 
#> tec50        0.141                            1.15           45.0       6.06 
#> tkout        -2.95                          0.0522           9.16       32.4 
#> te0           4.57                            96.6           5.24       18.1 
#> pdadd.err     3.72                            3.72                           
#>  
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#>   Information about run found ($runInfo):
#>    • NONMEM terminated due to rounding errors, but reading into nlmixr2/rxode2 anyway 
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     
#> 
#>  WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
#> 
#>  (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
#> 
#>     
#> 0MINIMIZATION TERMINATED
#>  DUE TO ROUNDING ERRORS (ERROR=134)
#>  NO. OF FUNCTION EVALUATIONS USED:     1088
#>  NO. OF SIG. DIGITS UNREPORTABLE
#> 0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
#> 
#>     IPRED relative difference compared to Nonmem IPRED: 0%; 95% percentile: (0%,0%); rtol=6.36e-06
#>     PRED relative difference compared to Nonmem PRED: 0%; 95% percentile: (0%,0%); rtol=6.08e-06
#>     IPRED absolute difference compared to Nonmem IPRED: 95% percentile: (2.53e-06, 0.000502); atol=7.15e-05
#>     PRED absolute difference compared to Nonmem PRED: 95% percentile: (3.79e-07,0.00321); atol=6.08e-06
#>     there are solving errors during optimization (see '$prderr')
#>     nonmem model: 'pk.turnover.emax3-nonmem/pk.turnover.emax3.nmctl' 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 483 × 35
#>   ID     TIME CMT      DV  PRED   RES IPRED   IRES  IWRES eta.ktr eta.ka eta.cl
#>   <fct> <dbl> <fct> <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>   <dbl>  <dbl>  <dbl>
#> 1 1       0.5 cp      0    1.16 -1.16 0.444 -0.444 -0.864  -0.506 -0.506  0.699
#> 2 1       1   cp      1.9  3.37 -1.47 1.45   0.446  0.840  -0.506 -0.506  0.699
#> 3 1       2   cp      3.3  7.51 -4.21 3.96  -0.660 -1.03   -0.506 -0.506  0.699
#> # ℹ 480 more rows
#> # ℹ 23 more variables: eta.v <dbl>, eta.emax <dbl>, eta.ec50 <dbl>,
#> #   eta.kout <dbl>, eta.e0 <dbl>, cp <dbl>, depot <dbl>, gut <dbl>,
#> #   center <dbl>, effect <dbl>, ktr <dbl>, ka <dbl>, cl <dbl>, v <dbl>,
#> #   emax <dbl>, ec50 <dbl>, kout <dbl>, e0 <dbl>, DCP <dbl>, PD <dbl>,
#> #   kin <dbl>, tad <dbl>, dosenum <dbl>
```

Which shows the `preds` and `ipreds` match between `NONMEM` and
`nlmixr2` quite well.

## Optional Step 3: Use `nlmixr2` to help understand why `NONMEM` failed

Since it *is* a `nlmixr2` fit, you can do interesting things with this
fit that you couldn’t do in `NONMEM` or even in another translator. For
example, if you wanted to add a covariance step you can with
[`getVarCov()`](https://rdrr.io/pkg/nlme/man/getVarCov.html):

``` r

getVarCov(f)
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → calculate sensitivities
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → calculate ∂(f)/∂(η)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → finding duplicate expressions in inner model...
#> → finding duplicate expressions in EBE model...
#> → compiling inner model...
#> ✔ done
#> → finding duplicate expressions in FD model...
#> → compiling EBE model...
#> ✔ done
#> → compiling events FD model...
#> ✔ done
#> calculating covariance matrix
#> [====|====|====|====|====|====|====|====|====|====] 0:04:07
#> Warning: "r,s" not usable (S not PD); installed "r"
#> Warning: full R matrix non-positive definite; using s (full)
#> → compress origData in nlmixr2 object, save 27800
#> Updated original fit object f
#>                      tktr           tka           tcl            tv
#> tktr         3.091125e+07 -3.091241e+07 -7.550188e+01  1.739784e+01
#> tka         -3.091241e+07  3.091357e+07  7.547697e+01 -1.736873e+01
#> tcl         -7.550188e+01  7.547697e+01  1.577469e-02 -1.661075e-02
#> tv           1.739784e+01 -1.736873e+01 -1.661075e-02  2.767842e-02
#> prop.err     3.726628e+01 -3.726674e+01 -1.187460e-03  2.549405e-03
#> pkadd.err   -1.858334e+02  1.858292e+02 -1.675658e-03  3.492400e-03
#> temax        8.117110e+03 -8.120875e+03 -2.396442e-01  7.944445e-01
#> tec50       -1.231431e+02  1.231456e+02 -6.967204e-03  1.872698e-02
#> tkout        8.440237e+01 -8.439416e+01 -1.418136e-03  1.232347e-03
#> te0          7.633257e+01 -7.633044e+01 -2.622076e-03  4.097490e-03
#> pdadd.err   -1.524528e+02  1.524549e+02  1.829493e-02 -3.450923e-02
#> om.eta.ktr   2.882030e+07 -2.882110e+07 -2.742672e+02  3.014618e+02
#> om.eta.ka   -2.883673e+07  2.883754e+07  2.743962e+02 -3.015951e+02
#> om.eta.cl    1.840196e+01 -1.839463e+01 -4.479628e-03  5.735584e-03
#> om.eta.v     8.220243e+01 -8.221068e+01  6.905915e-04 -6.308290e-04
#> om.eta.emax -5.234492e+03  5.235207e+03  1.623760e-01 -9.403670e-01
#> om.eta.ec50 -9.893199e+01  9.894168e+01 -1.606492e-04  1.148513e-05
#> om.eta.kout  3.045428e+00 -3.045739e+00  1.129859e-04 -2.994834e-04
#> om.eta.e0    3.728962e+00 -3.728516e+00 -1.177410e-04  1.164968e-04
#>                  prop.err     pkadd.err         temax         tec50
#> tktr         3.726628e+01 -1.858334e+02  8.117110e+03 -1.231431e+02
#> tka         -3.726674e+01  1.858292e+02 -8.120875e+03  1.231456e+02
#> tcl         -1.187460e-03 -1.675658e-03 -2.396442e-01 -6.967204e-03
#> tv           2.549405e-03  3.492400e-03  7.944445e-01  1.872698e-02
#> prop.err     4.263232e-04 -6.675116e-05  1.338641e-01  2.089755e-03
#> pkadd.err   -6.675116e-05  5.593140e-03  5.608511e-01  4.948372e-03
#> temax        1.338641e-01  5.608511e-01  4.372493e+02  2.565099e+00
#> tec50        2.089755e-03  4.948372e-03  2.565099e+00  3.964277e-02
#> tkout        1.267438e-04 -2.019854e-03 -8.282122e-01 -4.931336e-03
#> te0          4.902292e-04  7.195714e-05  6.236262e-02  1.782455e-03
#> pdadd.err   -3.324443e-03 -8.241485e-03 -3.034863e+00 -4.610010e-02
#> om.eta.ktr   4.815713e+01 -1.058084e+02  1.424545e+04  9.183667e+01
#> om.eta.ka   -4.818151e+01  1.058896e+02 -1.425342e+04 -9.187183e+01
#> om.eta.cl    5.292441e-04  8.445913e-04  1.905346e-01  3.797384e-03
#> om.eta.v     4.934218e-05 -3.168388e-05  2.276597e-01  7.826593e-05
#> om.eta.emax -1.163849e-01 -2.425263e-01 -3.437518e+01 -7.617707e-01
#> om.eta.ec50 -1.451874e-04  7.549761e-04 -3.237756e-01  2.782054e-03
#> om.eta.kout -4.949812e-05  4.051283e-05 -1.430888e-05 -4.716310e-04
#> om.eta.e0    1.189533e-05 -9.227448e-05 -1.949027e-02 -5.602469e-05
#>                     tkout           te0     pdadd.err    om.eta.ktr
#> tktr         8.440237e+01  7.633257e+01 -1.524528e+02  2.882030e+07
#> tka         -8.439416e+01 -7.633044e+01  1.524549e+02 -2.882110e+07
#> tcl         -1.418136e-03 -2.622076e-03  1.829493e-02 -2.742672e+02
#> tv           1.232347e-03  4.097490e-03 -3.450923e-02  3.014618e+02
#> prop.err     1.267438e-04  4.902292e-04 -3.324443e-03  4.815713e+01
#> pkadd.err   -2.019854e-03  7.195714e-05 -8.241485e-03 -1.058084e+02
#> temax       -8.282122e-01  6.236262e-02 -3.034863e+00  1.424545e+04
#> tec50       -4.931336e-03  1.782455e-03 -4.610010e-02  9.183667e+01
#> tkout        3.298995e-03  4.801553e-04  1.619224e-03  8.785654e+01
#> te0          4.801553e-04  1.198167e-03 -4.094263e-03  1.148822e+02
#> pdadd.err    1.619224e-03 -4.094263e-03  1.024544e-01 -6.076600e+02
#> om.eta.ktr   8.785654e+01  1.148822e+02 -6.076600e+02  3.324579e+07
#> om.eta.ka   -8.791160e+01 -1.149371e+02  6.079958e+02 -3.326417e+07
#> om.eta.cl    2.015034e-04  7.750589e-04 -8.320624e-03  7.236137e+01
#> om.eta.v    -4.173346e-04  1.352446e-04  2.395645e-04  8.185090e+01
#> om.eta.emax -1.034601e-01 -1.910109e-01  1.818383e+00 -1.755247e+04
#> om.eta.ec50 -4.714040e-04 -1.364992e-04 -2.677171e-03 -9.949625e+01
#> om.eta.kout -6.137740e-05 -4.909978e-05  3.467334e-04  1.467532e+00
#> om.eta.e0    6.428220e-05  4.906887e-05  3.895225e-05  3.623201e+00
#>                 om.eta.ka     om.eta.cl      om.eta.v   om.eta.emax
#> tktr        -2.883673e+07  1.840196e+01  8.220243e+01 -5.234492e+03
#> tka          2.883754e+07 -1.839463e+01 -8.221068e+01  5.235207e+03
#> tcl          2.743962e+02 -4.479628e-03  6.905915e-04  1.623760e-01
#> tv          -3.015951e+02  5.735584e-03 -6.308290e-04 -9.403670e-01
#> prop.err    -4.818151e+01  5.292441e-04  4.934218e-05 -1.163849e-01
#> pkadd.err    1.058896e+02  8.445913e-04 -3.168388e-05 -2.425263e-01
#> temax       -1.425342e+04  1.905346e-01  2.276597e-01 -3.437518e+01
#> tec50       -9.187183e+01  3.797384e-03  7.826593e-05 -7.617707e-01
#> tkout       -8.791160e+01  2.015034e-04 -4.173346e-04 -1.034601e-01
#> te0         -1.149371e+02  7.750589e-04  1.352446e-04 -1.910109e-01
#> pdadd.err    6.079958e+02 -8.320624e-03  2.395645e-04  1.818383e+00
#> om.eta.ktr  -3.326417e+07  7.236137e+01  8.185090e+01 -1.755247e+04
#> om.eta.ka    3.328255e+07 -7.239112e+01 -8.189212e+01  1.756096e+04
#> om.eta.cl   -7.239112e+01  2.159501e-03 -3.671314e-04 -6.620427e-02
#> om.eta.v    -8.189212e+01 -3.671314e-04  1.035927e-03 -4.847621e-02
#> om.eta.emax  1.756096e+04 -6.620427e-02 -4.847621e-02  1.449808e+02
#> om.eta.ec50  9.955873e+01  7.514084e-04 -1.054829e-03  9.582878e-02
#> om.eta.kout -1.467902e+00 -3.587803e-05  5.912877e-05  1.829018e-02
#> om.eta.e0   -3.625640e+00  9.619027e-06 -2.850592e-05  4.693969e-03
#>               om.eta.ec50   om.eta.kout     om.eta.e0
#> tktr        -9.893199e+01  3.045428e+00  3.728962e+00
#> tka          9.894168e+01 -3.045739e+00 -3.728516e+00
#> tcl         -1.606492e-04  1.129859e-04 -1.177410e-04
#> tv           1.148513e-05 -2.994834e-04  1.164968e-04
#> prop.err    -1.451874e-04 -4.949812e-05  1.189533e-05
#> pkadd.err    7.549761e-04  4.051283e-05 -9.227448e-05
#> temax       -3.237756e-01 -1.430888e-05 -1.949027e-02
#> tec50        2.782054e-03 -4.716310e-04 -5.602469e-05
#> tkout       -4.714040e-04 -6.137740e-05  6.428220e-05
#> te0         -1.364992e-04 -4.909978e-05  4.906887e-05
#> pdadd.err   -2.677171e-03  3.467334e-04  3.895225e-05
#> om.eta.ktr  -9.949625e+01  1.467532e+00  3.623201e+00
#> om.eta.ka    9.955873e+01 -1.467902e+00 -3.625640e+00
#> om.eta.cl    7.514084e-04 -3.587803e-05  9.619027e-06
#> om.eta.v    -1.054829e-03  5.912877e-05 -2.850592e-05
#> om.eta.emax  9.582878e-02  1.829018e-02  4.693969e-03
#> om.eta.ec50  8.605450e-03 -1.442539e-04  3.921287e-05
#> om.eta.kout -1.442539e-04  5.181594e-05 -6.005238e-06
#> om.eta.e0    3.921287e-05 -6.005238e-06  8.850947e-06
```

`nlmixr2` is more generous in what constitutes a covariance step. The
`r,s` covariance matrix is the “most” successful covariance step for
`focei`, but the system will fall back to other methods if necessary.

While this covariance matrix is not `r,s`, and should be regarded with
caution, it can still give us some clues on why this things are not
working in `NONMEM`.

When examining the fit, you can see the shrinkage is high for `temax`,
`tktr` and `tka`, so they could be dropped, making things more likely to
converge in `NONMEM`.

## Optional Step 4: Use model piping to get a successful NONMEM run

If we use model piping to remove the parameters, the new run will start
at the last model’s best estimates (saving a bunch of model development
time).

In this case, I specify the output directory `pk.turnover.emax4` with
the control and get the following:

``` r

f2 <- f %>% model(ktr <- exp(tktr)) %>%
  model(ka <- exp(tka)) %>%
  model(emax = expit(temax)) %>%
  nlmixr(data=nlmixr2data::warfarin, est="nonmem",
         control=nonmemControl(readRounding=FALSE,
                               modelName="pk.turnover.emax4"))
#> ! remove between subject variability `eta.ktr`
#> ! remove between subject variability `eta.ka`
#> ! remove between subject variability `eta.emax`
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → finding duplicate expressions in EBE model...
#> → compiling EBE model...
#> ✔ done
#> → Calculating residuals/tables
#> ✔ done
#> → compress origData in nlmixr2 object, save 27800
#> → compress parHistData in nlmixr2 object, save 7512
```

You can see the `NONMEM` run is now successful and validates against the
`rxode2` model below:

``` r

f2
```

``` math
\begin{align*}
{ktr} & = \exp\left({tktr}\right) \\
{ka} & = \exp\left({tka}\right) \\
{cl} & = \exp\left({tcl}+{eta.cl}\right) \\
{v} & = \exp\left({tv}+{eta.v}\right) \\
{emax} & = expit({temax}, {0}, {1}) \\
{ec50} & = \exp\left({tec50}+{eta.ec50}\right) \\
{kout} & = \exp\left({tkout}+{eta.kout}\right) \\
{e0} & = \exp\left({te0}+{eta.e0}\right) \\
{DCP} & = \frac{{center}}{{v}} \\
{PD} & = {1}-\frac{{emax} {\times} {DCP}}{\left({ec50}+{DCP}\right)} \\
effect({0}) & = {e0} \\
{kin} & = {e0} {\times} {kout} \\
\frac{d \: depot}{dt} & = -{ktr} {\times} {depot} \\
\frac{d \: gut}{dt} & = {ktr} {\times} {depot}-{ka} {\times} {gut} \\
\frac{d \: center}{dt} & = {ka} {\times} {gut}-\frac{{cl}}{{v}} {\times} {center} \\
\frac{d \: effect}{dt} & = {kin} {\times} {PD}-{kout} {\times} {effect} \\
{cp} & = \frac{{center}}{{v}} \\
{cp} & \sim prop({prop.err})+add({pkadd.err}) \\
{effect} & \sim add({pdadd.err}){\lor}{pca}
\end{align*}
```

One thing to emphasize: unlike other translators, you will know
immediately if the translation is off because the model will not
validate. Hence you can start this process with confidence - you will
know immediately if something is wrong.

This is related to [converting NONMEM to a nlmixr2
fit](https://nlmixr2.github.io/nonmem2rx/articles/convert-nlmixr2.html).

Since it is a `nlmixr2` object it would be easy to perform a VPC too
(the same is true for NONMEM models):

``` r

v1s <- vpcPlot(f2, show=list(obs_dv=TRUE), scales="free_y") +
  ylab("Warfarin Cp [mg/L] or PCA") +
  xlab("Time [h]")
#> Warning in filter_dv(obs, verbose): No software packages matched for filtering values, not filtering.
#>  Object class: other, data.frame
#>  Available filters: phoenix, nonmem
#> Warning in filter_dv(sim, verbose): No software packages matched for filtering values, not filtering.
#>  Object class: other, nlmixr2vpcSim, data.frame
#>  Available filters: phoenix, nonmem

v2s <- vpcPlot(f2, show=list(obs_dv=TRUE), pred_corr = TRUE, scales="free_y") +
  ylab("Prediction Corrected Warfarin Cp [mg/L] or PCA") +
  xlab("Time [h]")
#> Warning in filter_dv(obs, verbose): No software packages matched for filtering values, not filtering.
#>  Object class: other, data.frame
#>  Available filters: phoenix, nonmem
#> Warning in filter_dv(obs, verbose): No software packages matched for filtering values, not filtering.
#>  Object class: other, nlmixr2vpcSim, data.frame
#>  Available filters: phoenix, nonmem


library()

v1s
```

![](running-nonmem_files/figure-html/unnamed-chunk-2-1.png)

``` r


v2s
```

![](running-nonmem_files/figure-html/unnamed-chunk-2-2.png)
