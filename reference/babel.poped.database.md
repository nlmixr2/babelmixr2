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

## Author

Matthew L. Fidler
