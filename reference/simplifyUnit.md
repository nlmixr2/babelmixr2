# Simplify units by removing repeated units from the numerator and denominator

Simplify units by removing repeated units from the numerator and
denominator

## Usage

``` r
simplifyUnit(numerator = "", denominator = "")
```

## Arguments

- numerator:

  The numerator of the units (or the whole unit specification)

- denominator:

  The denominator of the units (or NULL if `numerator` is the whole unit
  specification)

## Value

The units specified with units that are in both the numerator and
denominator cancelled.

## Details

`NA` or `""` for `numerator` and `denominator` are considered unitless.

## See also

Other Unit conversion:
[`modelUnitConversion()`](https://nlmixr2.github.io/babelmixr2/reference/modelUnitConversion.md)

## Examples

``` r
simplifyUnit("kg", "kg/mL")
#> [1] "mL"
# units that don't match exactly are not cancelled
simplifyUnit("kg", "g/mL")
#> [1] "kg*mL/g"
```
