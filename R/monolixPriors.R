## Priors from the `ini({})` block, written as Monolix MAP estimation
##
## Monolix puts a prior on a population parameter by estimating it with
## `method=MAP` and describing the prior as a distribution of that
## parameter in the `[POPULATION]` section of `<MODEL>`:
##
##   [POPULATION]
##   DEFINITION:
##   ka_pop = {distribution=logNormal, typical=1.5, sd=0.5}
##
## The transform decision: Monolix's prior on a typical value has the same
## distribution as the parameter itself, and its `sd` is the standard
## deviation in the Gaussian space (ie `log(ka_pop) ~ N(log(1.5), 0.5^2)`).
## nlmixr2's `tka <- log(1.5)` with `prior(tka) ~ dnorm(log(1.5), 0.5)`
## is a normal prior on exactly that Gaussian-space value, so the prior
## *mean* is back-transformed like the estimate it describes (`exp()`,
## `expit()`, `probitInv()`) while the prior *sd* is written unchanged.
## Nothing is approximated -- there is no delta method, the two priors
## are the same distribution.  Covariate effects and residual error
## parameters are untransformed in both programs, so they are written as
## `distribution=normal` with the mean and sd as given.
##
## What Monolix cannot express is refused rather than dropped:
##
## - a prior on an omega or correlation element (Monolix's own MAP prior
##   on a standard deviation is an inverse Wishart, not a normal on the
##   value nlmixr2 describes)
## - a multivariate normal prior (Monolix priors are one per parameter)
## - a prior on a residual error parameter nlmixr2 estimates as a
##   variance, since Monolix estimates its square root and the normal
##   prior does not survive that transform
## - any prior that is not normal (`dcauchy()`, `invWishart()`, ...),
##   which the `"tnpri"` level also has nlmixr2est refuse up front

#' Stored priors of a ui as a data frame
#'
#' @param ui rxode2 ui
#' @return the `iniDf` rows that carry a prior (possibly none)
#' @noRd
#' @author Matthew L. Fidler
.mlxtranPriorRows <- function(ui) {
  .iniDf <- ui$iniDf
  if (!any(names(.iniDf) == "prior")) return(.iniDf[0, , drop=FALSE])
  .iniDf[!is.na(.iniDf$prior), , drop=FALSE]
}

#' Read the mean and sd of a univariate normal prior
#'
#' `lotri` stores normal priors positionally (`dnorm(mean, sd)`) or as
#' `stdNormal()`; the arguments may still be expressions like `log(2)`.
#'
#' @param prior prior as stored in the `prior` column
#' @param name parameter name (for the error)
#' @return numeric vector `c(mean, sd)`
#' @noRd
#' @author Matthew L. Fidler
.mlxtranPriorNormal <- function(prior, name) {
  .e <- try(str2lang(prior), silent=TRUE)
  .fn <- if (!inherits(.e, "try-error") && is.call(.e)) as.character(.e[[1]]) else ""
  if (length(.fn) != 1L) .fn <- ""
  if (.fn %in% c("stdNormal", "std_normal") && length(.e) == 1L) {
    return(c(0, 1))
  }
  if (.fn %in% c("dnorm", "normal") && length(.e) == 3L &&
        all(names(as.list(.e))[-1] %in% c("", NA_character_))) {
    .v <- try(vapply(as.list(.e)[-1], function(a) {
      as.double(eval(a, envir=baseenv()))
    }, double(1), USE.NAMES=FALSE), silent=TRUE)
    if (!inherits(.v, "try-error") && all(is.finite(.v)) && .v[2] > 0) {
      return(.v)
    }
  }
  stop("the prior on '", name, "' (", prior, ") cannot be written for Monolix; ",
       "only a univariate normal prior (dnorm(mean, sd) or stdNormal()) can be used with est=\"monolix\"",
       call.=FALSE)
}

#' Monolix prior information for every parameter that carries a prior
#'
#' @param ui rxode2 ui (decompressed)
#' @return data frame with the nlmixr2 `name`, Monolix parameter `par`,
#'   the `definition` line for the `[POPULATION]` section and `fix`; zero
#'   rows when there are no priors.  Errors for a prior Monolix cannot
#'   represent.
#' @noRd
#' @author Matthew L. Fidler
.mlxtranPriorInfo <- function(ui) {
  .ui <- ui
  .p <- .mlxtranPriorRows(.ui)
  .ret <- data.frame(name=character(0), par=character(0),
                     definition=character(0), fix=logical(0),
                     stringsAsFactors=FALSE)
  if (length(.p$name) == 0L) return(.ret)
  .omega <- which(!is.na(.p$neta1))
  if (length(.omega) > 0L) {
    stop("the model puts a prior on the omega parameter(s) ",
         paste0("'", .p$name[.omega], "'", collapse=", "),
         ", which cannot be written for Monolix (its MAP prior on a ",
         "standard deviation is not a normal prior on the omega)",
         call.=FALSE)
  }
  .split <- .ui$getSplitMuModel
  .muRef <- c(.split$pureMuRef, .split$taintMuRef)
  .covDataFrame <- .ui$saemMuRefCovariateDataFrame
  .curEval <- .ui$muRefCurEval
  .def <- vapply(seq_along(.p$name), function(i) {
    .name <- .p$name[i]
    .mv <- .mlxtranPriorNormal(.p$prior[i], .name)
    .w <- which(.covDataFrame$covariateParameter == .name)
    if (!is.na(.p$err[i])) {
      if (.mlxtranIsVarianceErr(.ui, .name)) {
        stop("the prior on '", .name, "' is on a variance, but Monolix ",
             "estimates the standard deviation, so it cannot be written for Monolix",
             call.=FALSE)
      }
      .par <- eval(str2lang(paste0("rxToMonolix(", .name, ", ui=.ui)")))
      .dist <- "distribution=normal"
      .typical <- .mv[1]
    } else if (length(.w) == 1L) {
      .par <- paste0("beta_", .muRef[.covDataFrame$theta[.w]], "_",
                     .covDataFrame$covariate[.w])
      .dist <- "distribution=normal"
      .typical <- .mv[1]
    } else {
      .var <- .muRef[.name]
      if (is.na(.var)) stop("babelmixr2 can't figure out '", .name, "' for monolix",
                            call.=FALSE)
      .par <- paste0(.var, "_pop")
      .wc <- which(.curEval$parameter == .name)
      .ce <- if (length(.wc) == 1L) paste(.curEval$curEval[.wc]) else ""
      .low <- if (length(.wc) == 1L) .curEval$low[.wc] else NA_real_
      .hi <- if (length(.wc) == 1L) .curEval$hi[.wc] else NA_real_
      if (is.na(.low) && !is.na(.hi)) .low <- 0
      .dist <- paste(c(.mlxTranCurEvalToDistribution(.ce),
                       .mlxTranGetLimits(.ce, .low, .hi)), collapse=", ")
      .typical <- .getNonMonolixParameterIni(.mv[1], .name, .curEval, .ui)
    }
    paste0(.par, " = {", .dist, ", typical=", .typical, ", sd=", .mv[2], "}")
  }, character(1), USE.NAMES=FALSE)
  data.frame(name=.p$name, par=sub(" = .*$", "", .def),
             definition=.def, fix=.p$fix, stringsAsFactors=FALSE)
}

#' @export
rxUiGet.mlxtranModelPopulation <- function(x, ...) {
  .ui <- x[[1]]
  .info <- .mlxtranPriorInfo(.ui)
  .info <- .info[!.info$fix, , drop=FALSE]
  if (length(.info$name) == 0L) return("")
  paste0("[POPULATION]\n",
         "DEFINITION:\n",
         paste(.info$definition, collapse="\n"), "\n\n")
}
attr(rxUiGet.mlxtranModelPopulation, "rstudio") <- "character"
