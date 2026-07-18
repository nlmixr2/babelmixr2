# nlmixr2 estimation using lme4::nlmer

nlmixr2 estimation using lme4::nlmer

## Usage

``` r
# S3 method for class 'nlmer'
nlmixr2Est(env, ...)
```

## Arguments

- env:

  Environment for the nlmixr2 estimation routines.

  This needs to have:

  \- rxode2 ui object in \`\$ui\`

  \- data to fit in the estimation routine in \`\$data\`

  \- control for the estimation routine's control options in \`\$ui\`

- ...:

  Other arguments provided to \`nlmixr2Est()\` provided for flexibility
  but not currently used inside nlmixr
