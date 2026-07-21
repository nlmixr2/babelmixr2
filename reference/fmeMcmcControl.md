# Control for fmeMcmc estimation method in nlmixr2

Control for fmeMcmc estimation method in nlmixr2

## Usage

``` r
fmeMcmcControl(
  jump = NULL,
  prior = NULL,
  niter = 1000L,
  outputlength = niter,
  burninlength = 0,
  updatecov = niter,
  covscale = NULL,
  ntrydr = 1,
  drscale = NULL,
  verbose = FALSE,
  seed = 99,
  returnFmeMcmc = FALSE,
  stickyRecalcN = 4,
  maxOdeRecalc = 5,
  odeRecalcFactor = 10^(0.5),
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
  covMethod = c("mcmc", "r", ""),
  adjObf = TRUE,
  ci = 0.95,
  sigdig = 4,
  sigdigTable = NULL,
  ...
)
```

## Arguments

- jump:

  jump length, either a *number*, a *vector* with length equal to the
  total number of parameters, a *covariance matrix*, or a *function*
  that takes as input the current values of the parameters and produces
  as output the perturbed parameters. See details.

- prior:

  -2\*log(parameter prior probability), either a function that is called
  as `prior(p)` or `NULL`; in the latter case a non-informative prior is
  used (i.e. all parameters are equally likely, depending on `lower` and
  `upper` within min and max bounds).

- niter:

  number of iterations for the MCMC.

- outputlength:

  number of iterations kept in the output; should be smaller or equal to
  `niter`.

- burninlength:

  number of initial iterations to be removed from output.

- updatecov:

  number of iterations after which the parameter covariance matrix is
  (re)evaluated based on the parameters kept thus far, and used to
  update the MCMC jumps.

- covscale:

  scale factor for the parameter covariance matrix, used to perform the
  MCMC jumps.

- ntrydr:

  maximal number of tries for the delayed rejection procedure. It is
  generally not a good idea to set this to a too large value.

- drscale:

  for each try during delayed rejection, the cholesky decomposition of
  the proposal matrix is scaled with this amount; if `NULL`, it is
  assumed to be `c(0.2,0.25, 0.333, 0.333, ...)`

- verbose:

  if `TRUE` or `1`: prints extra output, if numeric value `i > 1`,
  prints status information every `i` iterations.

- seed:

  an integer seed used for the mcmc chain. The chain is fully determined
  by this seed, so a fixed value makes the fit reproducible. The run is
  wrapped in
  [`rxode2::rxWithSeed()`](https://nlmixr2.github.io/rxode2/reference/rxWithSeed.html),
  which restores the prior random number generator state afterward, so
  seeding the fit does not disturb the calling session's stream. Vary
  this to explore different chains.

- returnFmeMcmc:

  return the fmeMcmc output instead of the nlmixr2 fit

- stickyRecalcN:

  The number of bad ODE solves before reducing the atol/rtol for the
  rest of the problem.

- maxOdeRecalc:

  Maximum number of times to reduce the ODE tolerances and try to
  resolve the system if there was a bad ODE solve.

- odeRecalcFactor:

  The ODE recalculation factor when ODE solving goes bad, this is the
  factor the rtol/atol is reduced

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

  Method for calculating the covariance. `"analytic"` (the default) uses
  the exact analytic observed-information R-matrix (reported as
  \\R^{-1}\\) and additionally returns the residual and `Omega` standard
  errors; it covers FOCEI/FOCE fits with additive, proportional, or
  combined error, mu-referenced/covariate/other structural parameters
  (and non-mu-referenced etas), and SD-scale inter-occasion variability,
  and emits a message and falls back to the finite-difference Hessian
  for anything out of scope (FO, `nAGQ > 1`, censoring, DV-transformed
  error, bounded-parameter transforms, a structural theta shared by two
  etas, non-SD `iovXform`, or a pure-proportional variance that vanishes
  at a near-zero prediction). The finite-difference methods use R (the
  Hessian) and S (the sum of individual gradient cross-products at the
  empirical Bayes estimates): `"r,s"` sandwich
  (`solve(R)%*%S%*%solve(R)`), `"r"` Hessian-based (`solve(R)`), `"s"`
  cross-product-based (`solve(S)`), or `""` to skip the covariance step.

- adjObf:

  is a boolean to indicate if the objective function should be adjusted
  to be closer to NONMEM's default objective function. By default this
  is `TRUE`

- ci:

  Confidence level for some tables. By default this is 0.95 or 95%
  confidence.

- sigdig:

  Optimization significant digits; controls the inner/outer optimization
  tolerance (`10^-sigdig`), ODE solver tolerance (`0.5*10^(-sigdig-2)`,
  or `0.5*10^(-sigdig-1.5)` for sensitivity/steady-state with liblsoda),
  and boundary check tolerance (`5*10^(-sigdig+1)`).

- sigdigTable:

  Significant digits in the final output table. If not specified, then
  it matches the significant digits in the \`sigdig\` optimization
  algorithm. If \`sigdig\` is NULL, use 3.

- ...:

  Ignored parameters

## Value

fmeMcmc control structure

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
   E0 <- 0.5
   Em <- 0.5
   E50 <- 2
   g <- fix(2)
 })
 model({
   v <- E0+Em*time^g/(E50^g+time^g)
   ll(bin) ~ DV * v - log(1 + exp(v))
 })
}

fit2 <- nlmixr(mod, dsn, est="fmeMcmc")
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → pruning branches (`if`/`else`) of population log-likelihood model...
#> ✔ done
#> → loading llik model into symengine environment...
#> → finding duplicate expressions in population log-likelihood model...
#> → optimizing duplicate expressions in population log-likelihood model...
#> ✔ done
#>  
#>  
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → finding duplicate expressions in EBE model...
#> → optimizing duplicate expressions in EBE model...
#> → compiling EBE model...
#>  
#>  
#> ✔ done
#> → Calculating residuals/tables
#> ✔ done
#> → compress origData in nlmixr2 object, save 8328
#> → compress parHistData in nlmixr2 object, save 10008

print(fit2)
#> ── nlmixr² log-likelihood fmeMcmc ──
#> 
#>           OBJF      AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> lPop -1275.217 568.6601 583.3834      -281.3301        200.1968        46.25382
#> 
#> ── Time (sec $time): ──
#> 
#>              setup    optimize covariance preprocess postprocess table compress
#> elapsed 0.03216184 0.002070661   5.66e-06      0.056       0.008 0.033    0.022
#>            other
#> elapsed 2.162762
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>        Est.     SE  %RSE   Back-transformed(95%CI)
#> E0  -0.7825 0.2755 35.20 -0.7825 (-1.322, -0.2426)
#> Em    4.844  1.075 22.20      4.844 (2.736, 6.952)
#> E50   2.320 0.3079 13.27      2.320 (1.717, 2.923)
#> g     2.000  FIXED FIXED                     2.000
#>  
#>   Covariance Type ($covMethod): mcmc
#>   Censoring ($censInformation): No censoring
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 1,000 × 5
#>   ID      TIME    DV  IPRED      v
#>   <fct>  <dbl> <dbl>  <dbl>  <dbl>
#> 1 1     0.0529     0 -0.377 -0.780
#> 2 1     0.0603     0 -0.378 -0.779
#> 3 1     0.0710     0 -0.378 -0.778
#> # ℹ 997 more rows

# you can also get the FME modMCMC output with

# fit2$fmeMcmc

# And use it in the summaries from FME, i.e.

summary(fit2$fmeMcmc)
#>              E0        Em      E50
#> mean -0.5654944 3.1879493 1.960768
#> sd    0.2754586 1.0753887 0.307860
#> min  -0.9487418 0.3892937 1.197740
#> max   0.5106981 4.8627463 2.773992
#> q025 -0.7265166 2.7256953 1.770652
#> q050 -0.6304920 3.3302785 1.941923
#> q075 -0.5180678 4.1146109 2.176642

pairs(fit2$fmeMcmc)
#> Warning: argument 1 does not name a graphical parameter
#> Warning: argument 1 does not name a graphical parameter
#> Warning: argument 1 does not name a graphical parameter


# and you can also use the coda package with `as.mcmc()`
coda::raftery.diag(coda::as.mcmc(fit2))
#> 
#> Quantile (q) = 0.025
#> Accuracy (r) = +/- 0.005
#> Probability (s) = 0.95 
#> 
#> You need a sample size of at least 3746 with these values of q, r and s

# }
```
