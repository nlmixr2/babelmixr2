## NONMEM `$PRIOR NWPRI` records from the `ini({})` block priors
##
## NWPRI gives a multivariate normal prior to the population
## parameters (`$THETAP` means, `$THETAPV` variance) and an inverse
## Wishart prior to each omega block (`$OMEGAP` scale, `$OMEGAPD`
## degrees of freedom).  It applies them to the *first* THETAs and the
## *first* omega blocks of the problem, so a model whose priors are not
## on its leading parameters is refused rather than given a prior on the
## wrong parameter.
##
## A prior that cannot be expressed this way (a Cauchy prior, a normal
## prior on an omega element, an inverse Wishart with its own scale
## matrix, ...) is refused here as well, so it is never silently dropped.

#' Environment a stored prior's arguments are evaluated in
#'
#' @return environment
#' @noRd
#' @author Matthew L. Fidler
.nonmemPriorEvalEnv <- function() {
  new.env(parent=baseenv())
}

#' Complain about a prior that cannot be written as NONMEM's NWPRI
#'
#' @param name parameter name(s) the prior is on
#' @param prior prior as stored in the `prior` column
#' @param why reason
#' @return nothing, called for the error
#' @noRd
#' @author Matthew L. Fidler
.nonmemPriorStop <- function(name, prior, why) {
  stop("cannot translate the prior on '", paste(name, collapse="', '"),
       "' (", prior, ") to NONMEM's $PRIOR NWPRI: ", why,
       call.=FALSE)
}

#' Evaluate one argument of a stored prior
#'
#' @param arg language object
#' @param name parameter name(s), for errors
#' @param prior prior, for errors
#' @return numeric value
#' @noRd
#' @author Matthew L. Fidler
.nonmemPriorEval <- function(arg, name, prior) {
  .v <- try(eval(arg, envir=.nonmemPriorEvalEnv()), silent=TRUE)
  if (inherits(.v, "try-error") || !is.numeric(.v) || any(!is.finite(.v))) {
    .nonmemPriorStop(name, prior, "its parameters could not be read back")
  }
  as.double(.v)
}

#' Does this model carry any `ini({})` prior?
#'
#' @param ui rxode2 ui
#' @return logical
#' @noRd
#' @author Matthew L. Fidler
.nonmemHasPriors <- function(ui) {
  .iniDf <- ui$iniDf
  if (is.null(.iniDf) || !any(names(.iniDf) == "prior")) return(FALSE)
  any(!is.na(.iniDf$prior))
}

#' The omega blocks as they are written to `$OMEGA`
#'
#' @param ui rxode2 ui
#' @return list of named matrices, in eta order
#' @noRd
#' @author Matthew L. Fidler
.nonmemOmegaBlocks <- function(ui) {
  lotri::lotriMatInv(lotri::lotriEst(lotri::as.lotri(ui$iniDf), drop=TRUE))
}

#' The `ini({})` priors as NONMEM's NWPRI needs them
#'
#' @param ui rxode2 ui
#' @return `NULL` when the model has no prior; otherwise a list with
#'
#' - `theta` the names of the THETAs with a prior (always
#'   `THETA(1)`..`THETA(k)`, possibly empty)
#'
#' - `thetaMean` the prior means of those THETAs
#'
#' - `thetaVar` the prior variance matrix of those THETAs
#'
#' - `omega` a list of the omega blocks with a prior (always the first
#'   blocks written to `$OMEGA`, possibly empty); the scale of each is
#'   the block's own initial estimate
#'
#' - `omegaDf` the degrees of freedom of each of those blocks
#'
#' @noRd
#' @author Matthew L. Fidler
.nonmemPriorSpec <- function(ui) {
  if (!.nonmemHasPriors(ui)) return(NULL)
  .p <- rxode2::rxUiPriors(ui)
  if (length(.p$name) == 0L) return(NULL)
  .iniDf <- ui$iniDf
  .theta <- .iniDf[!is.na(.iniDf$ntheta), ]
  .theta <- .theta[order(.theta$ntheta), ]
  .nt <- length(.theta$name)
  .mean <- setNames(rep(NA_real_, .nt), .theta$name)
  .var <- matrix(0, .nt, .nt, dimnames=list(.theta$name, .theta$name))
  .blocks <- .nonmemOmegaBlocks(ui)
  .blockNames <- lapply(.blocks, function(b) dimnames(b)[[1]])
  .blockDf <- rep(NA_real_, length(.blocks))
  # the prior each parameter is already covered by; a block or joint
  # prior is stored on every member, but a second, different prior on a
  # covered parameter would otherwise be silently dropped
  .done <- setNames(character(0), character(0))
  .mark <- function(names, prior) {
    .done[names] <<- prior
  }
  for (.i in seq_along(.p$name)) {
    .name <- .p$name[.i]
    .prior <- .p$prior[.i]
    if (.name %in% names(.done)) {
      if (!identical(.done[[.name]], .prior)) {
        .nonmemPriorStop(.name, .prior,
                         paste0("it is already covered by the prior '", .done[[.name]], "'"))
      }
      next
    }
    .e <- try(str2lang(.prior), silent=TRUE)
    if (inherits(.e, "try-error") || !is.call(.e)) {
      .nonmemPriorStop(.name, .prior, "the prior could not be parsed")
    }
    .fn <- deparse1(.e[[1]])
    .args <- as.list(.e)[-1]
    if (!is.na(.p$neta1[.i])) {
      # A prior on an omega element: only degrees of freedom on a
      # whole omega block has an NWPRI analogue
      if (.fn != "invWishart") {
        .nonmemPriorStop(.name, .prior,
                         "only 'invWishart(nu)' can be given to an omega block")
      }
      if (length(.args) != 1L) {
        .nonmemPriorStop(.name, .prior,
                         paste0("NWPRI uses the block's own initial estimates as ",
                                "the inverse Wishart scale, so only 'invWishart(nu)' ",
                                "is supported"))
      }
      .b <- which(vapply(.blockNames, function(n) n[1] == .name,
                         logical(1), USE.NAMES=FALSE))
      if (length(.b) != 1L) {
        .nonmemPriorStop(.name, .prior,
                         "it has to be on the first eta of an omega block")
      }
      .nu <- .nonmemPriorEval(.args[[1]], .name, .prior)
      .dim <- length(.blockNames[[.b]])
      if (length(.nu) != 1L || .nu <= .dim - 1) {
        .nonmemPriorStop(.blockNames[[.b]], .prior,
                         paste0("an inverse Wishart on a ", .dim, "x", .dim,
                                " block needs degrees of freedom greater than ",
                                .dim - 1))
      }
      .blockDf[.b] <- .nu
      .mark(.blockNames[[.b]], .prior)
      next
    }
    if (.fn == "dnorm") {
      .mean[.name] <- .nonmemPriorEval(.args[[1]], .name, .prior)
      .var[.name, .name] <- .nonmemPriorEval(.args[[2]], .name, .prior)^2
      .mark(.name, .prior)
    } else if (.fn == "stdNormal") {
      .mean[.name] <- 0
      .var[.name, .name] <- 1
      .mark(.name, .prior)
    } else if (.fn == "multiNormal") {
      .cov <- .args[[2]]
      if (!(is.call(.cov) && identical(.cov[[1]], quote(`lotri`)))) {
        .nonmemPriorStop(.name, .prior, "the covariance could not be read back")
      }
      .cov <- try(eval(.cov, envir=list(lotri=lotri::lotri),
                       enclos=.nonmemPriorEvalEnv()), silent=TRUE)
      if (inherits(.cov, "try-error") || !is.matrix(.cov)) {
        .nonmemPriorStop(.name, .prior, "the covariance could not be read back")
      }
      .nm <- dimnames(.cov)[[1]]
      if (!all(.nm %in% .theta$name)) {
        .nonmemPriorStop(.nm, .prior,
                         "a joint normal prior can only be on population parameters")
      }
      .mu <- .nonmemPriorEval(.args[[1]], .nm, .prior)
      .mean[.nm] <- rep_len(.mu, length(.nm))
      .var[.nm, .nm] <- .cov[.nm, .nm]
      .mark(.nm, .prior)
    } else {
      .nonmemPriorStop(.name, .prior,
                       paste0("'", .fn, "()' is not a normal prior, which is ",
                              "all NWPRI gives population parameters"))
    }
  }
  # NWPRI puts its priors on the first THETAs and the first omega
  # blocks, so those have to be the ones with a prior
  .wt <- unname(which(!is.na(.mean)))
  if (length(.wt) > 0L && !identical(.wt, seq_len(max(.wt)))) {
    .miss <- setdiff(seq_len(max(.wt)), .wt)
    stop("NONMEM's $PRIOR NWPRI gives priors to the first THETAs; ",
         "move the population parameter(s) with a prior (",
         paste0("'", names(.mean)[.wt], "'", collapse=", "),
         ") before the one(s) without (",
         paste0("'", names(.mean)[.miss], "'", collapse=", "),
         ") in the `ini({})` block",
         call.=FALSE)
  }
  .wo <- unname(which(!is.na(.blockDf)))
  if (length(.wo) > 0L && !identical(.wo, seq_len(max(.wo)))) {
    .miss <- setdiff(seq_len(max(.wo)), .wo)
    stop("NONMEM's $PRIOR NWPRI gives priors to the first omega blocks; ",
         "move the eta(s) with a prior (",
         paste0("'", unlist(.blockNames[.wo]), "'", collapse=", "),
         ") before the one(s) without (",
         paste0("'", unlist(.blockNames[.miss]), "'", collapse=", "),
         ") in the `ini({})` block",
         call.=FALSE)
  }
  .tn <- names(.mean)[.wt]
  list(theta=.tn,
       thetaMean=.mean[.tn],
       thetaVar=.var[.tn, .tn, drop=FALSE],
       omega=.blocks[.wo],
       omegaDf=.blockDf[.wo])
}

#' @export
rxUiGet.nonmemPriorSpec <- function(x, ...) {
  .nonmemPriorSpec(x[[1]])
}
attr(rxUiGet.nonmemPriorSpec, "rstudio") <- list()

#' @export
rxUiGet.nonmemPrior <- function(x, ...) {
  if (is.null(rxUiGet.nonmemPriorSpec(x, ...))) return("")
  "$PRIOR NWPRI\n\n"
}
attr(rxUiGet.nonmemPrior, "rstudio") <- "nonmemPrior"

#' @export
rxUiGet.nonmemPriorRecords <- function(x, ...) {
  .spec <- rxUiGet.nonmemPriorSpec(x, ...)
  if (is.null(.spec)) return("")
  .ui <- x[[1]]
  .sigdig <- rxode2::rxGetControl(.ui, "iniSigDig", 5)
  .ret <- character(0)
  .nt <- length(.spec$theta)
  if (.nt > 0L) {
    .iniDf <- .ui$iniDf
    .num <- .iniDf$ntheta[match(.spec$theta, .iniDf$name)]
    .t0 <- .nonmemThetaPad(c("$THETAP", rep("", .nt - 1L)))
    .t1 <- .nonmemThetaPad(paste0(signif(.spec$thetaMean, .sigdig), " FIX"))
    .t2 <- .nonmemThetaPad(paste(.num), TRUE)
    .t3 <- .nonmemThetaPad(.spec$theta)
    .ret <- c(.ret,
              paste0(paste(paste0(.t0, " (", .t1, ") ; ", .t2, " - ", .t3),
                           collapse="\n"), "\n\n"),
              paste0(.nonmemHandleOneOmega(.spec$thetaVar, .ui, rec="$THETAPV",
                                           fix=TRUE), "\n"))
  }
  .no <- length(.spec$omega)
  if (.no > 0L) {
    # each block's $OMEGAP is followed by its own $OMEGAPD, the way the
    # NONMEM examples write them
    .ret <- c(.ret,
              paste0(paste(vapply(seq_len(.no), function(i) {
                .b <- .spec$omega[[i]]
                paste0(.nonmemHandleOneOmega(.b, .ui, rec="$OMEGAP", fix=TRUE),
                       "$OMEGAPD (", signif(.spec$omegaDf[i], .sigdig), " FIX) ; ",
                       paste(dimnames(.b)[[1]], collapse=" "), "\n")
              }, character(1), USE.NAMES=FALSE), collapse=""), "\n"))
  }
  paste(.ret, collapse="")
}
attr(rxUiGet.nonmemPriorRecords, "rstudio") <- "nonmemPriorRecords"
