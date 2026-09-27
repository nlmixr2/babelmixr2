# Expand a babelmixr2 PopED database

Expand a babelmixr2 PopED database

## Usage

``` r
babel.poped.database(popedInput, ..., optTime = NA)
```

## Arguments

- popedInput:

  The babelmixr2 generated PopED database

- ...:

  other parameters sent to
  [`PopED::create.poped.database()`](https://andrewhooker.github.io/PopED/reference/create.poped.database.html)

- optTime:

  boolean to indicate if the global time indexer inside of babelmixr2 is
  reset if the times are different. By default this is `TRUE`. If
  `FALSE` you can get slightly better run times and possibly slightly
  different results. When `optTime` is `FALSE` the global indexer is
  reset every time the PopED rxode2 is setup for a problem or when a
  poped dataset is created. You can manually reset with
  [`popedMultipleEndpointResetTimeIndex()`](https://nlmixr2.github.io/babelmixr2/reference/popedMultipleEndpointResetTimeIndex.md)

## Value

babelmixr2 PopED database (with \$babelmixr2 in database)

## References

Fidler ML, Denney W, Harrold J, Hooijmaijers R, Papathanasiou T,
Schoemaker R, Taubert M, Trame M, Wilkins J (2024). babelmixr2 and
PopED: Quick Conversion of NONMEM, Monolix and nlmixr2/rxode2 Models to
PopED Optimal Design Analysis. American Conference on Pharmacometrics
(ACoP) 2024. [doi:10.70534/XUMG6226](https://doi.org/10.70534/XUMG6226)

## Author

Matthew L. Fidler
