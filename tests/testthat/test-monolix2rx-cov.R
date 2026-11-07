test_that("monolix2rx import with covariates reads Monolix's names (#269)", {
  skip_on_cran()
  skip_if_not_installed("monolix2rx")
  # V has a categorical tAPGAR covariate, k has none
  pb <- suppressMessages(monolix2rx::monolix2rx(
    system.file("cov/phenobarbital_project.mlxtran", package = "monolix2rx")
  ))
  expect_equal(.monolixMuRef(pb), c(V_pop = "V", k_pop = "k"))

  theta <- pb$monolixFullTheta
  expect_equal(theta[["V_pop"]], log(1.24689732935331))
  expect_equal(theta[["k_pop"]], log(0.00479852739901777))
  expect_equal(theta[["beta_V_tAPGAR_High"]], 0.260578209875729)
  expect_equal(theta[["beta_V_tAPGAR_Low"]], 0.127948846160619)
  expect_equal(theta[["b"]], 0.108861379047358)

  omega <- pb$monolixOmega
  expect_equal(omega["omega_V", "omega_V"], 0.401481325663058^2)
  expect_equal(omega["omega_k", "omega_k"], 0.186667045786903^2)

  eta <- pb$monolixEtaObf
  expect_setequal(names(eta), c("ID", "omega_V", "omega_k", "OBJI"))
  expect_equal(eta$omega_V[eta$ID == 1], 0.0889532)
  expect_equal(eta$omega_k[eta$ID == 1], -0.0732039)

  # the covariance keeps the covariate coefficients
  cov <- pb$monolixCovariance
  expect_true(all(
    c("V_pop", "beta_V_tAPGAR_High", "beta_V_tAPGAR_Low", "k_pop") %in% dimnames(cov)[[1]]
  ))

  # a continuous covariate on two parameters; rxode2 records the covariate of
  # beta_V_lw70 as its equation, log(wt/70), so the coefficient is read under
  # its own name
  w2 <- suppressMessages(monolix2rx::monolix2rx(
    system.file("cov/warfarin_covariate2_project.mlxtran", package = "monolix2rx")
  ))
  expect_equal(.monolixMuRef(w2), c(Tlag_pop = "Tlag", ka_pop = "ka", V_pop = "V", Cl_pop = "Cl"))
  theta <- w2$monolixFullTheta
  expect_equal(theta[["V_pop"]], log(8.0607521386391))
  expect_equal(theta[["beta_V_lw70"]], 0.898951301623242)
  expect_equal(theta[["beta_Cl_lw70"]], 0.595223300469559)
  expect_false(anyNA(theta))
  eta <- w2$monolixEtaObf
  expect_equal(eta$omega_Cl[eta$ID == 1], -0.148654)
})
