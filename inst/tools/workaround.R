## This is only for nlmixr
for (f in c("src/RcppExports.cpp")) {
  l <- readLines(f)
  w <- which(regexpr("^[#]include <RcppArmadillo.h>", l) != -1)
  if (length(w) == 1) {
    l <- l[-w]
    message("Excluding RcppArmadillo from ", f)
    writeLines(l, f)
  }
}
.in <- suppressWarnings(readLines("src/Makevars.in"))
if (.Platform$OS.type == "windows") {
  .makevars <- file("src/Makevars.win", "wb")
  .i <- "I"
} else {
  .makevars <- file("src/Makevars", "wb")
  if (file.exists("/etc/os-release")) {
    .os <- readLines("/etc/os-release")
    if (any(grepl("Pop!_OS", .os, fixed=TRUE))) {
      .i <- "isystem"
    } else {
      .i <- "I"
    }
  } else {
    .i <- "I"
  }
}

writeLines(gsub("@ISYSTEM@", .i, .in),
           .makevars)
close(.makevars)
