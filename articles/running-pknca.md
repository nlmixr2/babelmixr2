# Obtain initial estimates and unit conversions with PKNCA

## Introduction

Initial estimates for a compartmental population PK model can be
obtained using `babelmixr2` with the `"pknca"` estimation method. Also,
the central compartment scaling factor can be auto-generated based on
units for dosing, concentration measurement, desired volume of
distribution units, and time.

You do not need to perform NCA analysis by hand; the `"pknca"`
estimation method will perform NCA analysis using the
[PKNCA](https://humanpred.github.io/pknca/) package automatically.

The methods used for converting NCA calculations to parameter estimates
are described in the help for
[`nlmixr2Est.pknca()`](https://nlmixr2.github.io/babelmixr2/reference/nlmixr2Est.pknca.md).

## Initial example

Initial model setup is the same as for any other `nlmixr2` model. You
must load the `babelmixr2` library so that the
[`nlmixr()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2.html)
function recognizes `est = "pknca"`.

``` r

library(nlmixr2)
#> ── Attaching packages ───────────────────────────────────────── nlmixr2 6.0.0 ──
#> ★ lotri           1.0.5           ◯ nlmixr2lib      0.3.2      
#> ★ nlmixr2data     2.0.10          ◯ nlmixr2rpt      0.2.2      
#> ★ nlmixr2save     0.1.0           ◯ nlmixr2targets  0.1.0      
#> ★ nlmixr2est      7.0.0           ◯ nonmem2rx       0.1.10     
#> ★ nlmixr2extra    5.1.0           ◯ pmxNODE         0.1.0      
#> ★ nlmixr2plot     5.0.2           ◯ posologyr       1.2.8      
#> ★ rxode2          5.1.4           ◯ shinyMixR       0.5.3      
#> ◯ admixr2         0.2.0           ◯ pmxNODE         0.1.0      
#> ◯ babelmixr2      0.1.11.9000     ◯ FME             1.3.6.4    
#> ◯ ggPMX           1.3.2           ◯ PopED           0.7.0      
#> ◯ monolix2rx      0.0.6           ◯ nlmixr2auto     1.0.0      
#> ◯ nlmixr2auto     1.0.0           ◯ nlmixr2autoinit 1.0.1      
#> ◯ nlmixr2autoinit 1.0.1           ◯ xpose.nlmixr2   0.4.2
#> ── Conflicts ───────────────────────────────────────────── nlmixr2conflicts() ──
#> ✖ rxode2::boxCox()     masks nlmixr2est::boxCox()
#> ✖ coda::traceplot()    masks babelmixr2::traceplot(), nlmixr2plot::traceplot()
#> ✖ rxode2::yeoJohnson() masks nlmixr2est::yeoJohnson()
library(babelmixr2)
one.compartment <- function() {
  ini({
    tka <- log(1.57); label("Ka (1/hr)")
    tcl <- log(2.72); label("Cl (L/hr)")
    tv <- log(31.5); label("V (L)")
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7; label("additive residual error (mg/L)")
  })
  # and a model block with the error specification and model specification
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    vc <- exp(tv + eta.v)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / vc * center
    cp <- center / vc
    cp ~ add(add.sd)
  })
}
```

To use PKNCA to get initial estimates, use `est = "pknca"` instead of
one of the other `nlmixr2` estimation methods.

For unit conversions, provide the units to the
`control = pkncaControl()` argument. Unit conversions are only supported
when the units can be automatically converted; mass/volume can be
converted to any other mass/volume ratio, but mass to molar or molar to
mass cannot because there is not a single mass-to-molar conversion
factor.

``` r

prepared <-
  nlmixr2(
    one.compartment,
    data = theo_sd,
    est = "pknca",
    control = pkncaControl(concu = "ng/mL", doseu = "mg", timeu = "hr", volumeu = "L")
  )
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> Loading required namespace: testthat
#> ℹ change initial estimate (0.89314878960486) and upper/lower bound (-3.50655789731998 to 3.72508541597241) of `tka`
#> → significant model change detected
#> → removed from model: '$getSplitModel'
#> ℹ change initial estimate (8.41044546236311) and upper/lower bound (5.51439905878865 to 10.899462850803) of `tcl`
#> ℹ change initial estimate (10.5377244826318) and upper/lower bound (7.94567233496473 to 13.1050053785005) of `tv`
```

Now, you have the prepared model with updated initial estimates and the
NCA results embedded. You can see the new model and the PKNCA estimates
by looking at the `prepared$ui` (the model as interpreted by rxode2) and
`prepared$nca` (the PKNCAresults object).

Note that in the new model, the fixed effect initial estimates are
changed from their original values. The residual error and
between-subject variability are unchanged.

``` r

prepared$ui
```

``` math
\begin{align*}
{ka} & = \exp\left({tka}+{eta.ka}\right) \\
{cl} & = \exp\left({tcl}+{eta.cl}\right) \\
{vc} & = \exp\left({tv}+{eta.v}\right) \\
\frac{d \: depot}{dt} & = -{ka} {\times} {depot} \\
\frac{d \: center}{dt} & = {ka} {\times} {depot}-\frac{{cl}}{{vc}} {\times} {center} \\
{cp} & = \frac{{1000} {\times} {center}}{{vc}} \\
{cp} & \sim add({add.sd})
\end{align*}
```

``` r


knitr::knit_print(
  summary(prepared$nca)
)
#>  Interval Start Interval End  N AUClast (hr*ng/mL) Cmax (ng/mL)
#>               0           24 12        74.6 [24.3]            .
#>               0          Inf 12                  .  8.65 [17.0]
#>           Tmax (hr) CL (based on AUClast) (mg/(hr*ng/mL))
#>                   .                           4.22 [23.0]
#>  1.14 [0.630, 3.55]                                     .
#>  Vss (based on AUClast) (mg/(ng/mL)) Half-life (hr) AUCinf,obs (hr*ng/mL)
#>                          25.0 [18.5]              .                     .
#>                                    .    8.18 [2.12]            115 [28.4]
#>  Cmax (dose-normalized) ((ng/mL)/mg)
#>                                    .
#>                        0.0274 [18.1]
#> 
#> Caption: AUClast, Cmax, CL (based on AUClast), Vss (based on AUClast), AUCinf,obs, Cmax (dose-normalized): geometric mean and geometric coefficient of variation; Tmax: median and range; Half-life: arithmetic mean and standard deviation; N: number of subjects
```

From the updated model, you can perform estimation on the new model
object, as with any other model that has been created for `nlmixr2`:

``` r

fit <- nlmixr(prepared, data = theo_sd, est = "focei", control = list(print = 0))
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → calculate sensitivities
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → calculate ∂(f)/∂(η)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → calculate ∂(R²)/∂(η)
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → finding duplicate expressions in inner model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → optimizing duplicate expressions in inner model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00 
#> 
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00 
#> 
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → finding duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → optimizing duplicate expressions in EBE model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → compiling inner model...
#> ✔ done
#> → finding duplicate expressions in FD model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → optimizing duplicate expressions in FD model...
#> [====|====|====|====|====|====|====|====|====|====] 0:00:00
#> → compiling EBE model...
#> ✔ done
#> → compiling events FD model...
#> ✔ done
#> calculating covariance matrix
#> [====|====|====|====|====|====|====|====|====|====] 0:00:01 
#> done
#> → Calculating residuals/tables
#> ✔ done

fit
```

``` math
\begin{align*}
{ka} & = \exp\left({tka}+{eta.ka}\right) \\
{cl} & = \exp\left({tcl}+{eta.cl}\right) \\
{vc} & = \exp\left({tv}+{eta.v}\right) \\
\frac{d \: depot}{dt} & = -{ka} {\times} {depot} \\
\frac{d \: center}{dt} & = {ka} {\times} {depot}-\frac{{cl}}{{vc}} {\times} {center} \\
{cp} & = \frac{{1000} {\times} {center}}{{vc}} \\
{cp} & \sim add({add.sd})
\end{align*}
```

## Give PKNCA a different dataset or a completed NCA analysis

To get the initial estimate, `babelmixr2` automatically converts your
modeling dataset to the format the is needed for PKNCA, and the NCA is
automatically performed using all the data.

In some cases (e.g. studies with sparse data), NCA may not be feasible.
In those cases, you can provide a different dataset to PKNCA compared to
the full modeling dataset. Usually, the simplest method is to provide
single-dose, dense-sampling, dose-ranging data (i.e. the
single-ascending dose portion of the first-in-human study) to be
estimated.

To do this, give your data to PKNCA using the `ncaData` argument to
[`pkncaControl()`](https://nlmixr2.github.io/babelmixr2/reference/pkncaControl.md)
as follows:

``` r

# Choose a subset of the full dataset for NCA
dNCA <- theo_sd[theo_sd$ID <= 6, ]

preparedNcaData <-
  nlmixr2(
    one.compartment,
    data = theo_sd,
    est = "pknca",
    control = pkncaControl(concu = "ng/mL", doseu = "mg", timeu = "hr", volumeu = "L", ncaData = dNCA)
  )
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> ℹ change initial estimate (0.929027077269762) and upper/lower bound (-3.50655789731998 to 3.32136703319919) of `tka`
#> → significant model change detected
#> → removed from model: '$getSplitModel'
#> ℹ change initial estimate (8.3955404628088) and upper/lower bound (5.85241523541802 to 10.7637056987378) of `tcl`
#> ℹ change initial estimate (10.5377244826318) and upper/lower bound (7.94370069836702 to 13.1024358787022) of `tv`

preparedNcaData$ui
```

``` math
\begin{align*}
{ka} & = \exp\left({tka}+{eta.ka}\right) \\
{cl} & = \exp\left({tcl}+{eta.cl}\right) \\
{vc} & = \exp\left({tv}+{eta.v}\right) \\
\frac{d \: depot}{dt} & = -{ka} {\times} {depot} \\
\frac{d \: center}{dt} & = {ka} {\times} {depot}-\frac{{cl}}{{vc}} {\times} {center} \\
{cp} & = \frac{{1000} {\times} {center}}{{vc}} \\
{cp} & \sim add({add.sd})
\end{align*}
```

The initial estimates are now based on NCA calculated from the `dNCA`
dataset rather than the full `theo_sd` dataset.

If you already have NCA results calculated by PKNCA with the required
parameters (“tmax”, “cmax.dn”, and “cllast”), you can provide those
instead using the `pkncaControl(ncaResults)` argument.

## Model requirements

To update the initial estimates, the model must have parameters in the
[`model()`](https://nlmixr2.github.io/rxode2/reference/model.html) block
with the names that are expected by `est = "pknca"`. The expected names
are:

- `ka`
- `vc`
- `cl`
- `vp`
- `vp2`
- `q`
- `q2`

Any of those parameter names that are found in the
[`model()`](https://nlmixr2.github.io/rxode2/reference/model.html) block
will be automatically traced back to the initial conditions
([`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html) block),
and the parameter values will be updated. If the parameter is estimated
on the log scale, the updated parameter value will automatically be
converted to the log scale.
