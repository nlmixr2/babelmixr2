# Unit conversion for pharmacokinetic models

Unit conversion for pharmacokinetic models

## Usage

``` r
modelUnitConversion(
  dvu = NA_character_,
  amtu = NA_character_,
  timeu = NA_character_,
  volumeu = NA_character_
)
```

## Arguments

- dvu, amtu, timeu:

  The units for the DV, AMT, and TIME columns in the data

- volumeu:

  The units for the volume parameters in the model

## Value

A list with names for the units associated with each parameter ("amtu",
"clearanceu", "volumeu", "timeu", "dvu") and the numeric value to
multiply the modeled estimate (for example, `cp`) so that the model is
consistent with the data units.

## See also

Other Unit conversion:
[`simplifyUnit()`](https://nlmixr2.github.io/babelmixr2/reference/simplifyUnit.md)

## Examples

``` r
modelUnitConversion(dvu = "ng/mL", amtu = "mg", timeu = "hr", volumeu = "L")
#> Loading required namespace: testthat
#> $amtu
#> [1] "mg"
#> 
#> $clearanceu
#> [1] "L/h"
#> 
#> $volumeu
#> [1] "L"
#> 
#> $timeu
#> [1] "hr"
#> 
#> $dvu
#> [1] "ng/mL"
#> 
#> $cmtu
#> [1] "mg/L"
#> 
#> $dvConversion
#> [1] 1000
#> 
```
