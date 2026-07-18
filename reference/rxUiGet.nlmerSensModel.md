# Compiled nlmer sensitivity model (prediction + gradients)

Returns a list with:

- `thetaGrad`: compiled rxode2 model computing rx_pred\_ and
  rx\_\_sens_rx_pred\_\_BY_THETA_i\_\_\_ for each estimated parameter.
  When `eventSens = "jump"` this model carries rxode2's analytic
  event-jump sensitivities (loaded/activated by
  [`nlmixr2est::.nlmSetupEnv()`](https://nlmixr2.github.io/nlmixr2est/reference/dot-nlmSetupEnv.html)).

- `predOnly`: compiled rxode2 model for prediction only

- `eventTheta`: integer vector flagging `THETA[i]` that need Shi2021
  finite differences for dosing parameters (all 0 under
  `eventSens = "jump"`)

- `paramNames`: character vector of parameter names in `THETA[i]` order

## Usage

``` r
# S3 method for class 'nlmerSensModel'
rxUiGet(x, ...)
```

## Arguments

- x:

  rxUiGet list(ui)

- ...:

  additional arguments (currently ignored)

## Details

The list follows the `modelInfo` contract of
[`nlmixr2est::.nlmSetupEnv()`](https://nlmixr2.github.io/nlmixr2est/reference/dot-nlmSetupEnv.html)
so the nlm C machinery can load it once and keep it resident.
