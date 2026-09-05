# Build a linCmt() model whose linCmt parameters are the given names.
.pharmmlLinCmtUi <- function(pars) {
  .body <- paste0(pars, " <- exp(t", seq_along(pars), ")", collapse = "; ")
  .ini <- paste0("t", seq_along(pars), " <- 1", collapse = "; ")
  .txt <- sprintf(
    "function() { ini({%s; add.sd <- 1}); model({%s; linCmt() ~ add(add.sd)}) }",
    .ini, .body)
  rxode2::rxUiDecompress(eval(parse(text = .txt))())
}

# rxode2's own answer for the structure: solve and read back the state columns.
.pharmmlLinCmtStates <- function(ui, pars) {
  .e <- rxode2::et(amt = 100)
  .e <- rxode2::et(.e, seq(0, 24, 6))
  .p <- setNames(rep(1, length(pars)), paste0("t", seq_along(pars)))
  .s <- rxode2::rxSolve(ui, .e, params = .p, returnType = "data.frame")
  setdiff(names(.s), c("time", pars, "ipredSim", "sim"))
}

.pharmmlLinCmtCases <- list(
  list(pars = c("cl", "v"),                       ncmt = 1L, depot = FALSE, elim = "cl"),
  list(pars = c("cl", "vc"),                      ncmt = 1L, depot = FALSE, elim = "cl"),
  list(pars = c("kel", "v"),                      ncmt = 1L, depot = FALSE, elim = "k"),
  list(pars = c("k", "v"),                        ncmt = 1L, depot = FALSE, elim = "k"),
  list(pars = c("ka", "cl", "v"),                 ncmt = 1L, depot = TRUE,  elim = "cl"),
  list(pars = c("ka", "kel", "v"),                ncmt = 1L, depot = TRUE,  elim = "k"),
  list(pars = c("cl", "v", "q", "vp"),            ncmt = 2L, depot = FALSE, elim = "cl"),
  list(pars = c("cl", "v1", "q", "v2"),           ncmt = 2L, depot = FALSE, elim = "cl"),
  list(pars = c("ka", "cl", "v", "q", "vp"),      ncmt = 2L, depot = TRUE,  elim = "cl"),
  list(pars = c("k", "v", "k12", "k21"),          ncmt = 2L, depot = FALSE, elim = "k"),
  list(pars = c("cl", "v", "q", "vp", "q2", "vp2"), ncmt = 3L, depot = FALSE, elim = "cl"),
  list(pars = c("ka", "cl", "v", "q", "vp", "q2", "vp2"), ncmt = 3L, depot = TRUE, elim = "cl")
)

test_that("the linCmt classifier agrees with rxode2 on structure", {
  skip_on_cran()
  for (.c in .pharmmlLinCmtCases) {
    .ui <- .pharmmlLinCmtUi(.c$pars)
    .got <- .pharmmlLinCmtInfo(.ui)

    expect_equal(.got$ncmt, .c$ncmt, info = paste(.c$pars, collapse = "/"))
    expect_equal(.got$depot, .c$depot, info = paste(.c$pars, collapse = "/"))
    expect_equal(.got$elimination, .c$elim, info = paste(.c$pars, collapse = "/"))

    # cross-check against rxode2's own resolution: the states it allocates when
    # actually solving must match the structure the classifier inferred
    .states <- .pharmmlLinCmtStates(.ui, .c$pars)
    expect_equal("depot" %in% .states, .c$depot,
                 info = paste(.c$pars, collapse = "/"))
    expect_equal(sum(grepl("^(central|peripheral)", .states)), .c$ncmt,
                 info = paste(.c$pars, collapse = "/"))
  }
})

test_that("the classifier reads the depot flag out of linCmtFlg", {
  # linCmtFlg = numSens*100 + nLin*10 + depot; only the units digit is usable
  # at the UI level (see pharmml-plan.md section 4.2)
  for (.c in .pharmmlLinCmtCases) {
    .ui <- .pharmmlLinCmtUi(.c$pars)
    .flg <- rxode2::rxModelVars(.ui)$flags[["linCmtFlg"]]
    expect_equal(.flg %% 10 == 1, .c$depot, info = paste(.c$pars, collapse = "/"))
  }
})

test_that("the classifier rejects parameterisations it cannot map", {
  # Michaelis-Menten and transit absorption reach this code and are refused
  # here, by name.
  expect_error(.pharmmlLinCmtInfo(.pharmmlLinCmtUi(c("cl", "v", "vm", "km"))),
               "'vm' is not supported")
  expect_error(.pharmmlLinCmtInfo(.pharmmlLinCmtUi(c("cl", "v", "ka", "ktr"))),
               "'ktr' is not supported")
})

test_that("rxode2 rejects the parameterisations it owns before this code runs", {
  # These never reach .pharmmlLinCmtInfo(): rxode2 refuses them while building
  # the UI.  Asserted so that a future rxode2 change which starts *accepting*
  # them shows up here as a failure rather than as silently wrong PharmML.
  expect_error(.pharmmlLinCmtUi(c("cl", "vss")), "central volume")
  expect_error(.pharmmlLinCmtUi(c("cl", "v1", "q", "vp")), "volume style")
  expect_error(.pharmmlLinCmtUi(c("alpha", "beta", "v")), "b")
})


# ---------------------------------------------------------------------------
# Numerical validation of the classifier's *interpretation*.
#
# PharmML cannot be executed here, so the emitted XML cannot be solved and
# compared directly.  What can be checked is the step that carries the real
# risk: that the classifier reads cl/v/q/vp (etc.) as the compartmental
# structure rxode2 means by them.  Building the equivalent explicit-ODE model
# from the classifier's output and solving both is a genuine test of that
# reading -- misreading `kel` as `cl`, or `vp` as the central volume, moves the
# curve.  XML emission from a validated interpretation is then mechanical, and
# is covered by schema validation and the structural comparisons below.
#
# These models are built with plain `rxode2::rxode2()` rather than through an
# nlmixr2 ui.  That is deliberate: a ui-based harness was tried first and could
# not distinguish two structurally different models, so it would have passed
# whatever the classifier said.  The plain path demonstrably responds to the
# difference (see the mutation guard at the end of this file).
# ---------------------------------------------------------------------------

# Explicit-ODE equivalent of a linCmt() model, as rxode2 model text.
.pharmmlLinCmtOdeText <- function(info) {
  .v <- info$v
  .ke <- if (info$elimination == "cl") paste0("(", info$cl, ")/(", .v, ")") else info$k
  .lines <- character(0)
  .central <- paste0("-(", .ke, ")*central")

  for (.i in seq_along(info$peripheral)) {
    .p <- info$peripheral[[.i]]
    .st <- paste0("peripheral", .i)
    if (.p$style == "q") {
      .k1i <- paste0("(", .p$q, ")/(", .v, ")")
      .ki1 <- paste0("(", .p$q, ")/(", .p$v, ")")
    } else {
      .k1i <- .p$k1i
      .ki1 <- .p$ki1
    }
    .central <- paste0(.central, "-(", .k1i, ")*central+(", .ki1, ")*", .st)
    .lines <- c(.lines,
                sprintf("d/dt(%s) <- (%s)*central-(%s)*%s", .st, .k1i, .ki1, .st))
  }

  if (info$depot) {
    .lines <- c(sprintf("d/dt(depot) <- -(%s)*depot", info$ka),
                sprintf("d/dt(central) <- (%s)*depot%s", info$ka, .central),
                .lines)
  } else {
    .lines <- c(sprintf("d/dt(central) <- %s", .central), .lines)
  }
  paste(c(.lines, sprintf("cp <- central/(%s)", .v)), collapse = "\n")
}

.pharmmlOdeModel <- function(info, pars) {
  .txt <- paste(c(sprintf("param(%s)", paste(pars, collapse = ", ")),
                  .pharmmlLinCmtOdeText(info)), collapse = "\n")
  rxode2::rxode2(.txt)
}

# A plain linCmt() model over the given parameter names.  A `param()` block
# declares them, which is what lets rxode2 resolve the parameterisation; an
# `x <- x` self-assignment does not (rxode2 reports "Ambiguous 'kel'").
.pharmmlLinCmtPlain <- function(pars) {
  .txt <- paste(c(sprintf("param(%s)", paste(pars, collapse = ", ")),
                  "cp <- linCmt()"), collapse = "\n")
  rxode2::rxode2(.txt)
}

.pharmmlLinCmtParVals <- function(pars) {
  setNames(seq(0.7, by = 0.37, length.out = length(pars)), pars)
}

test_that("the classifier's reading of linCmt() matches rxode2 numerically", {
  skip_on_cran()
  .e <- rxode2::et(amt = 100)
  .e <- rxode2::et(.e, seq(0, 48, 2))

  for (.c in .pharmmlLinCmtCases) {
    .lab <- paste(.c$pars, collapse = "/")
    .info <- .pharmmlLinCmtInfo(.pharmmlLinCmtUi(.c$pars))
    .p <- .pharmmlLinCmtParVals(.c$pars)

    .lin <- rxode2::rxSolve(.pharmmlLinCmtPlain(.c$pars), .p, .e,
                            returnType = "data.frame")
    .ode <- rxode2::rxSolve(.pharmmlOdeModel(.info, .c$pars), .p, .e,
                            returnType = "data.frame")

    expect_lt(max(abs(.lin$cp - .ode$cp)) / max(abs(.ode$cp)), 1e-5)
  }
})

test_that("the numerical check can actually fail", {
  # Guard against a vacuous comparison: corrupt the classifier's output and the
  # equivalence test above must reject it.  Without this, a harness that solved
  # the same model twice would pass silently.
  skip_on_cran()
  .e <- rxode2::et(amt = 100)
  .e <- rxode2::et(.e, seq(0, 48, 2))
  .pars <- c("cl", "v", "q", "vp")
  .p <- .pharmmlLinCmtParVals(.pars)
  .info <- .pharmmlLinCmtInfo(.pharmmlLinCmtUi(.pars))
  .lin <- rxode2::rxSolve(.pharmmlLinCmtPlain(.pars), .p, .e, returnType = "data.frame")

  .mutations <- list(
    "peripheral volume read as central" = function(i) { i$peripheral[[1]]$v <- i$v; i },
    "clearance read as a rate constant" = function(i) { i$elimination <- "k"; i$k <- i$cl; i },
    "peripheral compartment dropped"    = function(i) { i$peripheral <- list(); i })

  for (.nm in names(.mutations)) {
    .bad <- rxode2::rxSolve(.pharmmlOdeModel(.mutations[[.nm]](.info), .pars), .p, .e,
                            returnType = "data.frame")
    expect_gt(max(abs(.lin$cp - .bad$cp)) / max(abs(.bad$cp)), 1e-3)
  }
})
