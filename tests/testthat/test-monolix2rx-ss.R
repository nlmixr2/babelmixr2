test_that("monolix2rx imports with few steady-state doses get rxode2's minSS/maxSS floor", {
  skip_on_cran()
  skip_if_not_installed("monolix2rx")
  .ss <- function(n) {
    testthat::local_mocked_bindings(.getNbdoses=function(x) n, .package="monolix2rx")
    .c <- suppressMessages(.monolix2rxToFoceiControl(list(etaMat=NULL), NULL))
    c(.c$rxControl$minSS, .c$rxControl$maxSS)
  }
  expect_equal(.ss(3L), c(5L, 7L))
  expect_equal(.ss(5L), c(5L, 7L))
  expect_equal(.ss(7L), c(7L, 8L))
  expect_equal(.ss(10L), c(10L, 11L))
})
