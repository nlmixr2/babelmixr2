#' @export
rxUiGet.monolixJacobian <- function(x, ...) {
  # jacobian for lognormal
  #
  # D(S("exp(x)"),"x") = exp(x)
  # If the parameter comes from the rxode2 ui, then this would be:
  # exp(x); from monolix it would be exp(log(x)) = x

  # jacobian for expit()
  # exp(-x)*(hi - low)/(1 + exp(-x))^2

  # probit for probitNormal
  # 0.5*sqrt(2)*exp((-1/2)*x^2)*(hi - low)/sqrt(pi)
  #
  # This can be used to transform Monolix's covariance matrix to
  # nlmixr's style of covariance matrix
  #
  .ui <- x[[1]]
  .muRef <- .ui$muRefCurEval
  .iniDf <- .ui$iniDf
  .vars <- .iniDf[!is.na(.iniDf$ntheta) & is.na(.iniDf$err), c("est", "name")]
  .name <- .vars$name
  .est <- .vars$est
  .ret <- vapply(seq_along(.est), function(i) {
    .x <- .est[i]
    .n <- .name[i]
    .w <- which(.muRef$parameter == .n)
    if (length(.w) != 1) return(1)
    .low <- .muRef$low[.w]
    if (is.na(.low)) .low <- 0
    .hi <- .muRef$hi[.w]
    if (is.na(.hi)) .hi <- 1
    .v <- paste(.muRef$curEval[.w])
    switch(.v,
           exp=exp(.x),
           expit=exp(-.x)*(.hi - .low)/(1 + exp(-.x))^2,
           probitInv=0.5*sqrt(2)*exp((-1/2)*.x^2)*(.hi - .low)/sqrt(pi),
           1)
  }, double(1), USE.NAMES=FALSE)
  .ret <- diag(.ret)
  # This is used according to the table https://monolix.lixoft.com/tasks/result-files-generated-monolix/
  #
  # fim type | 2018        | 2019        | 2020        | 2021      |
  # ---------+-------------+-------------+-------------+-----------+
  # SA       | untransform | untransform | transform   | transform |
  # Lin      | untransform | untransform | untransform | transform |
  # ---------+-------------+-------------+-------------+-----------+
  dimnames(.ret) <- list(.name, .name)
  .ret
}
attr(rxUiGet.monolixJacobian, "rstudio") <- "monolixJacobian"
