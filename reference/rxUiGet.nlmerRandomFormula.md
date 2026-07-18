# Build the random-effects part (part 3) of the nlmer 3-part formula

Build the random-effects part (part 3) of the nlmer 3-part formula

## Usage

``` r
# S3 method for class 'nlmerRandomFormula'
rxUiGet(x, ...)
```

## Arguments

- x:

  rxUiGet list(ui)

- ...:

  additional arguments (currently ignored)

## Value

Formula: param1 + param2 + ... + (rand1 + rand2 \| ID)
