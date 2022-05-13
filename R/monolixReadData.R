#' @export
rxUiGet.monolixOutputVersion <- function(x, ...) {
  .exportPath <- rxUiGet.monolixExportPath(x, ...)
  .summary <- file.path(.exportPath, "summary.txt")
  if (file.exists(.summary)) {
    .lines <- readLines(.summary, n=5)
    .w <- which(regexpr(".*[vV]ersion *: *[^ ]*.*", .lines) != -1)
    if (length(.w) == 1) {
      .line <- .lines[.w]
      # 2019r is 5.1.1
      return(sub(".*[vV]ersion *: *([^ ]*).*", "\\1", .line))
    }
  }
  NULL
}
