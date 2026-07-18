# Control for saemix estimation method in nlmixr2

Control for saemix estimation method in nlmixr2

## Usage

``` r
saemixControl(
  map = TRUE,
  fim = TRUE,
  ll.is = TRUE,
  ll.gq = FALSE,
  nbiter.saemix = c(300, 100),
  nbiter.sa = NA,
  nbiter.burn = 5,
  nbiter.map = 5,
  nb.chains = 1,
  fix.seed = TRUE,
  seed = 23456,
  nmc.is = 5000,
  nu.is = 4,
  print.is = FALSE,
  nbdisplay = 100,
  displayProgress = FALSE,
  print = FALSE,
  save = FALSE,
  save.graphs = TRUE,
  directory = "newdir",
  warnings = FALSE,
  nbiter.mcmc = c(2, 2, 2, 0),
  proba.mcmc = 0.4,
  stepsize.rw = 0.4,
  rw.init = 0.5,
  alpha.sa = 0.97,
  nnodes.gq = 12,
  nsd.gq = 4,
  maxim.maxiter = 100,
  nb.sim = 1000,
  nb.simpred = 100,
  ipar.lmcmc = 50,
  ipar.rmcmc = 0.05,
  rxControl = NULL,
  stickyRecalcN = 4,
  maxOdeRecalc = 5,
  odeRecalcFactor = 10^(0.5),
  useColor = NULL,
  printNcol = NULL,
  normType = c("rescale2", "mean", "rescale", "std", "len", "constant"),
  scaleType = c("none", "nlmixr2", "norm", "mult", "multAdd"),
  scaleCmax = 1e+05,
  scaleCmin = 1e-05,
  scaleC = NULL,
  scaleTo = 1,
  addProp = c("combined2", "combined1"),
  calcTables = TRUE,
  compress = TRUE,
  ci = 0.95,
  sigdigTable = NULL,
  sigdig = 4,
  ...
)
```

## Arguments

- map:

  logical, whether to compute MAP estimates (default `TRUE`)

- fim:

  logical, whether to compute FIM (default `TRUE`)

- ll.is:

  logical, whether to compute log-likelihood by Importance Sampling
  (default `TRUE`)

- ll.gq:

  logical, whether to compute log-likelihood by Gaussian Quadrature
  (default `FALSE`)

- nbiter.saemix:

  integer vector of length 2, number of iterations for phase 1
  (exploratory) and phase 2 (smoothing) of SAEM (default `c(300, 100)`)

- nbiter.sa:

  number of iterations for simulated annealing (default `NA`)

- nbiter.burn:

  number of iterations for burn-in (default `5`)

- nbiter.map:

  number of iterations for MAP estimation (default `5`)

- nb.chains:

  number of chains (default `1`)

- fix.seed:

  logical, whether to fix the random seed (default `TRUE`)

- seed:

  random seed (default `23456`)

- nmc.is:

  number of Monte Carlo samples for Importance Sampling (default `5000`)

- nu.is:

  number of degrees of freedom for the student distribution in
  Importance Sampling (default `4`)

- print.is:

  logical, whether to print progress during Importance Sampling (default
  `FALSE`)

- nbdisplay:

  number of iterations between progress displays (default `100`)

- displayProgress:

  logical, whether to display graphical progress (default `FALSE`)

- print:

  whether `saemix` prints algorithm progress (default `FALSE`). Accepts
  the legacy logical/integer scalar or a pre-built
  [`nlmixr2est::iterPrintControl()`](https://nlmixr2.github.io/nlmixr2est/reference/iterPrintControl.html)
  object; internally the control stores a single `iterPrintControl`
  sub-list (matching the `nlmixr2est` iteration-print unification), and
  any nonzero `every` enables the `saemix` progress output

- save:

  logical, whether to save results to files (default `TRUE`)

- save.graphs:

  logical, whether to save graphs (default `TRUE`)

- directory:

  directory where results and graphs are saved (default `"newdir"`)

- warnings:

  logical, whether to show warnings (default `FALSE`)

- nbiter.mcmc:

  integer vector of length 4, number of iterations for MCMC kernel
  updates (default `c(2, 2, 2, 0)`)

- proba.mcmc:

  probability for MCMC kernel selection (default `0.4`)

- stepsize.rw:

  stepsize for random walk kernel (default `0.4`)

- rw.init:

  initial standard deviation for random walk kernel (default `0.5`)

- alpha.sa:

  parameter for simulated annealing (default `0.97`)

- nnodes.gq:

  number of nodes for Gaussian Quadrature (default `12`)

- nsd.gq:

  number of standard deviations for Gaussian Quadrature (default `4`)

- maxim.maxiter:

  maximum number of iterations for maximization step (default `100`)

- nb.sim:

  number of simulations for visual predictive check (default `1000`)

- nb.simpred:

  number of simulations for predictions (default `100`)

- ipar.lmcmc:

  parameter for L-MCMC (default `50`)

- ipar.rmcmc:

  parameter for R-MCMC (default `0.05`)

- rxControl:

  \`rxode2\` ODE solving options during fitting, created with
  \`rxControl()\`

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

- ci:

  Confidence level for some tables. By default this is 0.95 or 95%
  confidence.

- sigdigTable:

  Significant digits in the final output table. If not specified, then
  it matches the significant digits in the \`sigdig\` optimization
  algorithm. If \`sigdig\` is NULL, use 3.

- sigdig:

  Optimization significant digits; controls the inner/outer optimization
  tolerance (`10^-sigdig`), ODE solver tolerance (`0.5*10^(-sigdig-2)`,
  or `0.5*10^(-sigdig-1.5)` for sensitivity/steady-state with liblsoda),
  and boundary check tolerance (`5*10^(-sigdig+1)`).

- ...:

  Ignored parameters

## Value

saemix control structure

## Author

Matthew L. Fidler & Antigravity

## Examples

``` r
# \donttest{
# Example showing control creation:
ctrl <- saemixControl(seed = 123456)
# }
```
