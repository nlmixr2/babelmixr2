# Creating a New Estimation Method

``` r

library(babelmixr2)
```

This currently talks about the bare minimum needed to create an
estimation method. It will be updated with how to create an estimation
method that returns a nlmixr2 fit object later.

To create a new estimation method, the following steps are required.

## Create a `nlmixr2Est()` method

This method will have an input of the environment of the nlmixr2est UI
object (see
[`?nlmixr2Est`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2Est.html)).
It will output the fit object.

## Create a control method

The control method gives access to any controls that are required for
estimation.
