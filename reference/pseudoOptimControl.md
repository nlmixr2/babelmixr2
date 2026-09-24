# Control for fmeMcmc estimation method in nlmixr2

Control for fmeMcmc estimation method in nlmixr2

## Usage

``` r
pseudoOptimControl(
  npop = NULL,
  numiter = 10000,
  centroid = 3,
  varleft = NULL,
  verbose = FALSE,
  returnPseudoOptim = FALSE,
  stickyRecalcN = 4,
  maxOdeRecalc = 5,
  odeRecalcFactor = 10^(0.5),
  indTolRelax = TRUE,
  useColor = NULL,
  printNcol = NULL,
  print = 1L,
  normType = c("rescale2", "mean", "rescale", "std", "len", "constant"),
  scaleType = c("none", "nlmixr2", "norm", "mult", "multAdd"),
  scaleCmax = 1e+05,
  scaleCmin = 1e-05,
  scaleC = NULL,
  scaleTo = 1,
  rxControl = NULL,
  optExpression = TRUE,
  sumProd = FALSE,
  literalFix = TRUE,
  literalFixRes = TRUE,
  addProp = c("combined2", "combined1"),
  calcTables = TRUE,
  compress = TRUE,
  covMethod = c("r", ""),
  adjObf = TRUE,
  ci = 0.95,
  sigdig = 4,
  sigdigTable = NULL,
  eventSens = c("jump", "fd"),
  ...
)
```

## Arguments

- npop:

  Number of elements in the population. Defaults to max(5\*length(p),50)
  which is calculated from the number of parameters in the model

- numiter:

  Number of iterations to run the optimization. Defaults to 10000. The
  algorithm either stops when `numiter` iterations has been performed or
  when the remaining variation is less than `varleft`.

- centroid:

  Number of elements from which to estimate a new parameter vector. The
  default is 3.

- varleft:

  relative variation remaining; if below this value, the algorithm
  stops. When `NULL` (the default) it is derived from `sigdig`, keeping
  the
  [`FME::pseudoOptim()`](https://rdrr.io/pkg/FME/man/pseudoOptim.html)
  default of `1e-8` at the default `sigdig = 4` and tightening/loosening
  it one order of magnitude per significant digit.

- verbose:

  If TRUE, print information about the optimization from
  [`FME::pseudoOptim`](https://rdrr.io/pkg/FME/man/pseudoOptim.html).
  Default is FALSE.

- returnPseudoOptim:

  return the pseudoOptim output instead of the nlmixr2 fit

- stickyRecalcN:

  The number of bad ODE solves before reducing the atol/rtol for the
  rest of the problem.

- maxOdeRecalc:

  Maximum number of times to reduce the ODE tolerances and try to
  resolve the system if there was a bad ODE solve.

- odeRecalcFactor:

  The ODE recalculation factor when ODE solving goes bad, this is the
  factor the rtol/atol is reduced

- indTolRelax:

  when `TRUE` (default) a subject whose ODE solve had to be retried with
  a relaxed tolerance keeps that relaxed tolerance for the rest of the
  fit instead of resetting it every evaluation

- useColor:

  Logical (or \`NULL\`) emit ANSI bold/color escapes in the iteration
  print. \`NULL\` (default) defers to \[crayon::has_color()\].

- printNcol:

  Integer (or \`NULL\`) parameter columns per row before wrapping.
  \`NULL\` (default) uses \`floor((getOption("width") - 23) / 12)\`.

- print:

  Either a scalar print-frequency (\`0\` = suppress, \`1\` (default) =
  every evaluation, \`N\` = every Nth), OR a pre-built
  \[iterPrintControl()\] object. Equivalent to \`iterPrintControl(every
  = print, ncol = printNcol, useColor = useColor)\`.

- normType:

  Parameter normalization/scaling used to get scaled initial values for
  `scaleType`, of the form `Vscaled = (Vunscaled-C1)/C2` (see [Feature
  Scaling](https://en.wikipedia.org/wiki/Feature_scaling); `rescale2`
  follows the
  [OptdesX](http://apmonitor.com/me575/uploads/Main/optimization_book.pdf)
  manual): `"rescale2"` scales all parameters to (-1, 1); `"rescale"`
  (min-max) scales to (0, 1); `"mean"` centers on the mean with range
  (0, 1); `"std"` standardizes by mean/sd; `"len"` scales to unit
  (Euclidean) length; `"constant"` performs no normalization (`C1=0`,
  `C2=1`).

- scaleType:

  The scaling scheme for nlmixr2: `"nlmixr2"` (default) scales as
  `(current-init)*scaleC[i] + scaleTo`, with `scaleTo` from `normType`
  and scales from `scaleC`; `"norm"` uses the simple scaling from
  `normType`; `"mult"` scales multiplicatively as
  `current/init*scaleTo`; `"multAdd"` scales linearly
  (`(current-init)+scaleTo`) for parameters in an exponential block
  (e.g. `exp(theta)`) and multiplicatively otherwise.

- scaleCmax:

  Maximum value of the scaleC to prevent overflow.

- scaleCmin:

  Minimum value of the scaleC to prevent underflow.

- scaleC:

  Scaling constant used with `scaleType="nlmixr2"`; when not specified,
  chosen by parameter type to keep gradient sizes similar on a log
  scale: \`1\` for exp()-transformed/power/boxCox/ yeoJohnson
  parameters, \`0.5\*abs(est)\` for additive/proportional/ lognormal
  error parameters, \`abs(1/digamma(est+1))\` for factorials, and
  \`log(abs(est))\*abs(est)\` for log-scale parameters. May be set
  explicitly per parameter if these defaults don't apply well.

- scaleTo:

  Scale the initial parameter estimate to this value. By default this
  is 1. When zero or below, no scaling is performed.

- rxControl:

  \`rxode2\` ODE solving options during fitting, created with
  \`rxControl()\`

- optExpression:

  Optimize the rxode2 expression to speed up calculation. By default
  this is turned on.

- sumProd:

  Is a boolean indicating if the model should change multiplication to
  high precision multiplication and sums to high precision sums using
  the PreciseSums package. By default this is `FALSE`.

- literalFix:

  boolean, substitute fixed population values as literals and re-adjust
  ui and parameter estimates after optimization; Default is \`TRUE\`.

- literalFixRes:

  boolean, substitute fixed population values as literals and re-adjust
  ui and parameter estimates after optimization; Default is \`TRUE\`.

- addProp:

  Type of additive-plus-proportional error: \`"combined1"\`, where
  standard deviations add: \$\$y = f + (a + b\times f^c) \times
  \varepsilon\$\$; or \`"combined2"\`, where variances add: \$\$y = f +
  \sqrt{a^2 + b^2\times f^{2\times c}} \times \varepsilon\$\$. Here y =
  observed, f = predicted, a = additive sd, b = proportional/power sd, c
  = power exponent (1 in the proportional case).

- calcTables:

  This boolean is to determine if the foceiFit will calculate tables. By
  default this is `TRUE`

- compress:

  Should the object have compressed items

- covMethod:

  Method for calculating the covariance. `"r,s"` (the default) is the
  sandwich estimator (see below). `"analytic"` uses the exact analytic
  observed-information R-matrix (reported as \\R^{-1}\\) and
  additionally returns the residual and `Omega` standard errors; it
  covers FOCEI/FOCE fits with additive, proportional, or combined error,
  mu-referenced/covariate/other structural parameters (and
  non-mu-referenced etas), and SD-scale inter-occasion variability, and
  emits a message and falls back to the finite-difference Hessian for
  anything out of scope (FO, `nAGQ > 1`, censoring, DV-transformed
  error, bounded-parameter transforms, a structural theta shared by two
  etas, non-SD `iovXform`, or a pure-proportional variance that vanishes
  at a near-zero prediction). The finite-difference methods use R (the
  Hessian) and S (the sum of individual gradient cross-products at the
  empirical Bayes estimates): `"r,s"` sandwich
  (`solve(R)%*%S%*%solve(R)`), `"r"` Hessian-based (`solve(R)`), `"s"`
  cross-product-based (`solve(S)`), or `""` to skip the covariance step.
  `"sa"` (SAEM Louis stochastic-approximation FIM) and `"imp"`
  (importance-sampling Monte-Carlo observed information) are also
  accepted for any method; they are computed post-fit at the converged
  estimates by the decoupled recompute engine.

- adjObf:

  is a boolean to indicate if the objective function should be adjusted
  to be closer to NONMEM's default objective function. By default this
  is `TRUE`

- ci:

  Confidence level for some tables. By default this is 0.95 or 95%
  confidence.

- sigdig:

  Optimization significant digits. One value drives, with a single
  consistent formula, the inner/outer optimizer convergence tolerance
  (`10^-sigdig`), the boundary check tolerance (`5*10^(-sigdig+1)`), and
  the ODE solver tolerances: the `rtol` exponent IS `sigdig` and `atol`
  sits three orders below, so `rtol = 10^-sigdig`,
  `atol = 10^(-sigdig-3)` for every solver (stiff, non-stiff or
  auto-switching). The sensitivity (`atolSens`/`rtolSens`) tolerances
  match the main solve (the outer gradient and covariance are built from
  them); the steady-state (`ssAtol`/`ssRtol`) tolerances run one order
  looser. Keying the optimizer to the same `10^-sigdig` means it
  converges to exactly the precision the solve supports. At the default
  `sigdig = 3` this is `atol = 1e-6`, `rtol = 1e-3`.

- sigdigTable:

  Significant digits in the final output table. If not specified, then
  it matches the significant digits in the \`sigdig\` optimization
  algorithm. If \`sigdig\` is NULL, use 3.

- eventSens:

  method used for the dosing-parameter (alag/F/rate/dur) sensitivities:
  `"jump"` routes them through rxode2's analytic event jumps; `"fd"`
  falls back to Shi2021 finite differences. See
  [`nlmixr2est::nlmControl()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmControl.html).

- ...:

  Ignored parameters

## Value

pseudoOptim control structure

## Author

Matthew L. Fidler

## Examples

``` r

# \donttest{
# A logit regression example with emax model

dsn <- data.frame(i=1:1000)
dsn$time <- exp(rnorm(1000))
dsn$DV=rbinom(1000,1,exp(-1+dsn$time)/(1+exp(-1+dsn$time)))

mod <- function() {
 ini({
   # This estimation method requires all parameters
   # to be bounded:
   E0 <- c(-100, 0.5, 100)
   Em <- c(0, 0.5, 10)
   E50 <- c(0, 2, 20)
   g <- fix(c(0.1, 2, 10))
 })
 model({
   v <- E0+Em*time^g/(E50^g+time^g)
   ll(bin) ~ DV * v - log(1 + exp(v))
 })
}

fit2 <- nlmixr(mod, dsn, est="pseudoOptim")
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → pruning branches (`if`/`else`) of population log-likelihood model...
#> ✔ done
#> → loading into symengine environment...
#> → finding duplicate expressions in population log-likelihood model...
#> ✔ done
#>  
#>  
#> → calculating covariance
#> ✔ done
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → finding duplicate expressions in EBE model...
#> → compiling EBE model...
#>  
#>  
#> ✔ done
#> → Calculating residuals/tables
#> ✔ done
#> → compress origData in nlmixr2 object, save 8360
#> → compress parHistData in nlmixr2 object, save 30136

print(fit2)
#> ── nlmixr² log-likelihood pseudoOptim ──
#> 
#>           OBJF     AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> lPop -696.6074 1147.27 1161.993      -570.6348        399.9608        62.45302
#> 
#> ── Time (sec $time): ──
#> 
#>              setup    optimize covariance preprocess postprocess table compress
#> elapsed 0.01835472 0.002471776  5.069e-06      0.054       0.007 0.021    0.031
#>            other
#> elapsed 1.025168
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>        Est.     SE  %RSE    Back-transformed(95%CI)
#> E0  -0.6405 0.1210 18.90 -0.6405 (-0.8777, -0.4033)
#> Em    5.490  1.192 21.72       5.490 (3.153, 7.827)
#> E50   2.859 0.6150 21.51       2.859 (1.654, 4.065)
#> g     2.000  FIXED FIXED                      2.000
#>  
#>   Covariance Type ($covMethod): r
#>   Some strong fixed parameter correlations exist ($cor) :
#>      cor:Em,E0 cor:E50,E0 cor:E50,Em 
#>     0.353      0.618      0.914  
#>  
#> 
#>   Censoring ($censInformation): No censoring
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 1,000 × 5
#>   ID      TIME    DV  IPRED      v
#>   <fct>  <dbl> <dbl>  <dbl>  <dbl>
#> 1 1     0.0595     0 -0.424 -0.638
#> 2 1     0.0766     1 -1.06  -0.637
#> 3 1     0.0843     0 -0.425 -0.636
#> # ℹ 997 more rows


# }
```
