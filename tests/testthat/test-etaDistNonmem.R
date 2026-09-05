test_that("a declared non-Gaussian random effect translates to NONMEM", {
  # By the time this is reached, nlmixr2est's pre-processing hook has run
  # rxode2::rxEtaDistExpand(), so the declaration is ordinary model code: a
  # latent standard normal with a fixed unit $OMEGA, the copula correlation in
  # ordinary thetas, and a phiU() + inverse CDF line per declared random
  # effect.  phiU() is NONMEM's PHI(x)+DEL idiom and tanh() is DTANH(), so a
  # family with an ELEMENTARY quantile function needs nothing special written.
  .weibull <- function() {
    ini({
      lscale <- log(5)
      lshape <- log(2)
      tv <- 3.45
      tka <- 0.45
      dist(eta.cl) ~ dweibull(shape = exp(lshape), scale = exp(lscale))
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- eta.cl
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - (cl / v) * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }
  .u <- rxode2::rxUiDecompress(rxode2::rxEtaDistExpand(.weibull()))
  .nm <- paste(.u$nonmemModel, collapse = "\n")

  # the latent random effect is a unit-variance FIXED omega
  expect_true(grepl("$OMEGA", .nm, fixed = TRUE))
  expect_match(.nm, "1\\s+FIX")
  # phiU() becomes NONMEM's PHI() with the tail offset that keeps the inverse
  # CDF finite -- phi() saturates at |z| > 8.3 in double precision
  expect_true(grepl("PHI(", .nm, fixed = TRUE))
  expect_true(grepl("1.0E-15", .nm, fixed = TRUE))
})

test_that("a family NONMEM has no expression form for is refused by name", {
  # gamma/beta/Student t reach NONMEM only through `$ABBR FUNCTION
  # GAMMACDFINV(VQ,10)` plus a VQ argument vector filled slot by slot -- a
  # sequence of statements, not an expression, so it cannot come out of the
  # expression translator.  Refusing names the random effect the user wrote and
  # what NONMEM would need, rather than translating into a different model or
  # failing deep inside with "gammapInv is not supported".
  .gamma <- function() {
    ini({
      lclm <- log(5)
      lclrv <- log(0.09)
      tv <- 3.45
      tka <- 0.45
      dist(eta.cl) ~ dgamma(shape = 1 / exp(lclrv),
                            rate = 1 / (exp(lclrv) * exp(lclm)))
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- eta.cl
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - (cl / v) * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }
  .u <- rxode2::rxUiDecompress(rxode2::rxEtaDistExpand(.gamma()))
  expect_error(.u$nonmemModel, "GAMMACDFINV")
})
