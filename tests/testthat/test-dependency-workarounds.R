# Work-arounds for rxode2 features newer than the oldest rxode2 babelmixr2
# supports (DESCRIPTION).  Off CRAN, each block warns once the installed rxode2
# is a release that has its feature, and names what to remove then; before
# that it is quiet.  A new work-around adds a block here, and a row to
# .depBlocks so the files it touches are checked against what it names.

# The warning a work-around's block gives, or NULL unless the block's probe
# found the feature (`present`) and the installed version is a release after
# `dev`, the development version that added it.  A development build
# (x.y.z.9000 and later) is never a release, and a release without the
# feature still needs the work-around.
.depWorkaroundWarning <- function(installed, dev, remove, present) {
  if (!isTRUE(present) || !.depIsRelease(installed) || installed <= dev) {
    return(NULL)
  }
  paste0(
    "rxode2 ",
    installed,
    " is a release after ",
    dev,
    ", so this work-around can go: ",
    remove,
    "; and raise DESCRIPTION to rxode2 (>= ",
    installed,
    ")"
  )
}

# A release rather than a development version (a fourth component of 9000 or
# more)
.depIsRelease <- function(version) {
  .p <- unlist(numeric_version(version))
  length(.p) < 4L || .p[[4L]] < 9000L
}

# Does the file use any of `symbols`?
.depFileUses <- function(file, symbols) {
  .l <- readLines(file, warn = FALSE)
  any(vapply(symbols, .depLinesUse, logical(1), lines = .l))
}

.depLinesUse <- function(symbol, lines) {
  any(grepl(symbol, lines, fixed = TRUE))
}

# The files under R/ and src/ (when the source tree is there) and
# tests/testthat that use any of `symbols`, other than this file
.depFilesUsing <- function(symbols) {
  .root <- testthat::test_path("..", "..")
  .dirs <- c(file.path(.root, c("R", "src")), testthat::test_path())
  .f <- unlist(lapply(
    .dirs[dir.exists(.dirs)],
    list.files,
    pattern = "\\.(R|c|cpp|h|hpp)$",
    full.names = TRUE
  ))
  .f <- .f[basename(.f) != "test-dependency-workarounds.R"]
  basename(.f[vapply(.f, .depFileUses, logical(1), symbols = symbols)])
}

# stand-ins for the version mocks
.depVersion519 <- function() numeric_version("5.1.9")
.depVersionDev <- function() numeric_version("5.1.9.9000")

.depSlot104Remove <- paste(
  "rxode2's setIndNewind() (slot 104): remove bmSetIndNewind from src/rxode2Optional.h",
  "and its row of rxode2OptionalSlots in src/poped.cpp, the NULL check in popedCalcLhs()",
  "(src/poped.cpp, which then calls setIndNewind()), the skips in test-poped-newind.R,",
  "and this block; if no entry is left, also the rest of the rxode2Optional machinery:",
  "src/rxode2Optional.h, iniRxode2Optional() and rxode2Optional() in src/poped.cpp,",
  "their registration in src/init.c, .rxode2Optional() and its call in R/zzz.R, and",
  "its test in test-dependency-workarounds.R"
)

# The slot-104 block's warning for the installed rxode2
.depSlot104Warning <- function(present = .rxode2Optional()[["setIndNewind"]]) {
  .depWorkaroundWarning(.rxode2Version(), "5.1.8.9000", .depSlot104Remove, present = present)
}

# Each block: what its warning names to remove, the symbols of its probe,
# skips and fallbacks, and a file the scan must find (anchor), so a scan run
# from the wrong root fails instead of passing on what little it saw
.depBlocks <- list(
  slot104 = list(
    remove = .depSlot104Remove,
    anchor = "zzz.R",
    symbols = c("bmSetIndNewind", "rxode2Optional")
  )
)

test_that("the work-around warning is given only for a release with the feature", {
  expect_null(.depWorkaroundWarning(numeric_version("5.1.8"), "5.1.8.9000", "x", present = TRUE))
  expect_null(.depWorkaroundWarning(numeric_version("5.1.8.9000"), "5.1.8.9000", "x", present = TRUE))
  # a development build after the development version is not a release
  expect_null(.depWorkaroundWarning(numeric_version("5.1.8.9001"), "5.1.8.9000", "x", present = TRUE))
  expect_null(.depWorkaroundWarning(numeric_version("5.1.9.9000"), "5.1.8.9000", "x", present = TRUE))
  # a release without the feature still needs the work-around
  expect_null(.depWorkaroundWarning(numeric_version("5.1.9"), "5.1.8.9000", "x", present = FALSE))
  expect_null(.depWorkaroundWarning(numeric_version("5.1.9"), "5.1.8.9000", "x", present = NA))
  expect_match(
    .depWorkaroundWarning(numeric_version("5.1.9"), "5.1.8.9000", "x", present = TRUE),
    "rxode2 5.1.9 is a release after 5.1.8.9000.*x; and raise DESCRIPTION to rxode2 \\(>= 5.1.9\\)"
  )
  expect_true(.depIsRelease("5.1.9"))
  expect_true(.depIsRelease("6.0.0.1"))
  expect_false(.depIsRelease("5.1.8.9000"))
})

test_that(".rxode2Optional() reports setIndNewind from rxode2's pointer table", {
  .p <- rxode2::.rxode2ptrs()
  .has <- length(.p) > 104L && names(.p)[105L] == "rxode2setIndNewind"
  expect_identical(.rxode2Optional(), c(setIndNewind = .has))
  on.exit(.Call(`_babelmixr2_iniRxode2Optional`, .p, PACKAGE = "babelmixr2"), add = TRUE)
  # a table that ends before slot 104, or names another entry there, loads none
  .Call(`_babelmixr2_iniRxode2Optional`, .p[1:99], PACKAGE = "babelmixr2")
  expect_identical(.rxode2Optional(), c(setIndNewind = FALSE))
  if (.has) {
    .q <- .p
    names(.q)[105L] <- "rxode2other"
    .Call(`_babelmixr2_iniRxode2Optional`, .q, PACKAGE = "babelmixr2")
    expect_identical(.rxode2Optional(), c(setIndNewind = FALSE))
    .Call(`_babelmixr2_iniRxode2Optional`, .p, PACKAGE = "babelmixr2")
    expect_identical(.rxode2Optional(), c(setIndNewind = TRUE))
  }
})

test_that("slot-104 work-around: warns only for a release with setIndNewind()", {
  local_mocked_bindings(.rxode2Version = .depVersion519)
  expect_null(.depSlot104Warning(present = FALSE))
  expect_match(.depSlot104Warning(present = TRUE), "rxode2 5.1.9 is a release after 5.1.8.9000", fixed = TRUE)
  local_mocked_bindings(.rxode2Version = .depVersionDev)
  expect_null(.depSlot104Warning(present = TRUE))
})

test_that("slot-104 work-around: rxode2's setIndNewind()", {
  skip_on_cran()
  .w <- .depSlot104Warning()
  if (!is.null(.w)) {
    warning(.w, call. = FALSE)
  }
  succeed()
})

# after every block has added its row
test_that("each work-around names every file that uses its probe, skips and fallbacks", {
  # tests/ alone (an installed-package check) cannot show what R/ uses
  skip_if_not(
    dir.exists(file.path(testthat::test_path("..", ".."), "R")),
    "no R/ beside tests/: the file-name check needs the source tree"
  )
  for (.b in names(.depBlocks)) {
    .files <- .depFilesUsing(.depBlocks[[.b]]$symbols)
    expect_true(.depBlocks[[.b]]$anchor %in% .files, label = paste(.b, "scan finds its anchor"))
    for (.f in .files) {
      expect_true(grepl(.f, .depBlocks[[.b]]$remove, fixed = TRUE), label = paste(.b, "names", .f))
    }
  }
})
