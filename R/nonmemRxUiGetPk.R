.nonmemGetMuName <- function(theta, ui) {
  .muRefDf <- ui$muRefDataFrame
  .iniDf <- ui$iniDf
  .w <- which(.muRefDf$theta == theta)
  if (length(.w) != 1) return(NA_character_)
  .muRefDf$eta[.w]
}

.nonmemGetMuNum0 <- function(theta, ui) {
  .eta <- .nonmemGetMuName(theta, ui)
  if (is.na(.eta)) return(NA_real_)
  .iniDf <- ui$iniDf
  .w <- which(.iniDf$name == .eta)
  if (length(.w) != 1) return(NA_real_)
  .iniDf$neta1[.w]
}

.nonmemGetMuNum <- function(theta, ui) {
  .neta <- .nonmemGetMuNum0(theta, ui)
  if (is.na(.neta)) return(NA_character_)
  .muRef <- rxode2::rxGetControl(ui, "muRef", TRUE)
  paste0(ifelse(.muRef, "MU_", "UM_"), .neta)
}

.nonmemGetThetaNum <- function(theta, ui) {
  .iniDf <- ui$iniDf
  .w <- which(.iniDf$name == theta)
  if (length(.w) != 1) return(NA_character_)
  if (is.na(.iniDf$ntheta[.w])) return(NA_character_)
  paste0("THETA(", .iniDf$ntheta[.w], ")")
}

.nonmemGetEtaNum <- function(eta, ui) {
  .iniDf <- ui$iniDf
  .w <- which(.iniDf$name == eta)
  if (length(.w) != 1) return(NA_character_)
  if (is.na(.iniDf$neta1[.w])) return(NA_character_)
  paste0("ETA(", .iniDf$neta1[.w], ")")
}

.nonmemGetThetaMuCov <- function(theta, ui, covRefDf) {
  .w <- which(covRefDf$theta == theta)
  if (length(.w) == 0) return(NA_character_)
  # the covariate uses its NONMEM name (as in $INPUT)
  paste(paste0(vapply(covRefDf$covariate[.w], .nmGetVar, character(1),
                      ui=ui, USE.NAMES=FALSE), "*",
               vapply(covRefDf$covariateParameter[.w], .nonmemGetThetaNum, character(1),
                      ui=ui)),
        collapse="+")
}

#'@export
rxUiGet.nonmemThetaRep <- function(x, ...) {
  .ui <- x[[1]]
  .split <- .ui$getSplitMuModel
  .muRef <- c(.split$pureMuRef, .split$taintMuRef)
  .thetas <- names(.muRef)
  .covRefDf <- .ui$saemMuRefCovariateDataFrame
  .ret <- data.frame(theta=.thetas,
                     nmTheta=vapply(.thetas, .nonmemGetThetaNum, character(1), ui=.ui,
                                    USE.NAMES=FALSE),
             mu=vapply(.thetas, .nonmemGetMuNum, character(1), ui=.ui,
                       USE.NAMES=FALSE),
             cov=vapply(.thetas, .nonmemGetThetaMuCov, character(1),
                        ui=.ui, covRefDf=.covRefDf, USE.NAMES=FALSE))
  .ret$nmEta <- ifelse(is.na(.ret$mu), NA_character_,
                       paste0("ETA(",substr(.ret$mu,4, 10),")"))
  .ret
}
attr(rxUiGet.nonmemThetaRep, "rstudio") <- "nonmemThetaRep"

#'@export
rxUiGet.nonmemPkDesErr0 <- function(x, ...) {
  .ui <- x[[1]]
  .bblLinCmtAssertOde(.ui, "nonmem")
  rxode2::rxAssignControlValue(.ui, ".nmVarResNum", 1)
  rxode2::rxAssignControlValue(.ui, ".nmGetVarReservedDf",
                               data.frame(var=character(0),
                                          nm=character(0)))
  .advan <- .nonmemLinCmtAdvan(.ui)
  rxode2::rxAssignControlValue(.ui, ".nmLinCmtReserved", .advan$reserved)
  .split <- .ui$getSplitMuModel
  .mu <- rxUiGet.nonmemThetaRep(x, ...)
  .ret <- vapply(seq_along(.mu$mu), function(i) {
    if (is.na(.mu$mu[i])) return(NA_character_)
    paste0("  ", .mu$mu[i], "=", .mu$nmTheta[i],
           ifelse(is.na(.mu$cov[i]), "",
                  paste0("+", .mu$cov[i])))
  }, character(1), USE.NAMES=FALSE)
  .ret <- paste(.ret[!is.na(.ret)], collapse="\n")
  .mu2 <- setNames(ifelse(is.na(.mu$mu),
                          .mu$nmTheta,
                          paste0(.mu$mu, "+", .mu$nmEta)),
                   .mu$theta)
  assign(".thetaMu", .mu2, envir=.ui)
  on.exit({
    if (exists(".thetaMu", envir=.ui)) {
      rm(".thetaMu", envir=.ui)
    }
  })
  .isPred <- (length(rxode2::rxState(.ui)) == 0)
  .muRefDef <- .split$muRefDef
  .pk <- paste0(ifelse(.isPred,"$PRED\n","$PK\n"),
                 .ret,"\n",
                 paste(vapply(seq_along(.muRefDef),
                              function(i) {
                                x <-.rxToNonmem(.muRefDef[[i]], ui=.ui)
                                x
                              }, character(1), USE.NAMES=FALSE),
                       collapse="\n"))
  rm(".thetaMu", envir=.ui)

  .desModel <- .split$modelWithDrop[-.ui$predDf$line]
  .rmModel <- which(vapply(seq_along(.desModel),
                           function(i) {
                             identical(.desModel[[i]], quote(`_drop`))
                           }, logical(1), USE.NAMES=FALSE))
  if (length(.rmModel) > 0L) .desModel <- .desModel[-.rmModel]
  if (!is.null(.advan)) {
    # closed-form ADVAN: the state independent lines are calculated in
    # $PK before the micro-constants; the rest are in $ERROR
    .dep <- .bblLinCmtStateDep(.desModel, rxode2::rxState(.ui))
    # compartment properties (like alag(depot)) are written in $PK, but
    # still need translating to be recorded; the ODEs are kept so the
    # compartments are defined, and their DADT lines are dropped later
    .pkModel <- .desModel[!.dep$dep]
  }

  .mainModel <- rxode2::rxCombineErrorLines(.ui,
                                            errLines=nmGetDistributionNonmemLines(.ui),
                                            paramsLine=NA,
                                            modelVars=TRUE,
                                            cmtLines=FALSE,
                                            dvidLine=FALSE,
                                            lstExpr=.split$modelWithDrop,
                                            useIf=FALSE)
  .mv <- rxode2::rxModelVars(paste(vapply(seq_along(.desModel),
                                          function(i) {
                                            deparse1(.desModel[[i]])
                                          }, character(1), USE.NAMES=FALSE),
                                   collapse="\n"))
  .normMain <- strsplit(rxode2::rxNorm(eval(.mainModel)), "\n")[[1]]
  .normMainL <- vapply(seq_along(.normMain),
                       function(i) {
                         regexpr("^((alag|f|F|rate|dur|lag)[(][^)]+[)]|[^(]+[(]0[)]|d[/]dt[(][^(]+[)])=", .normMain[i]) == -1
                       }, logical(1), USE.NAMES=FALSE)

  .normMain <- paste(.normMain[.normMainL], collapse="\n")
  .lhs <- vapply(.mv$lhs,
                 function(v) {
                   paste0("RXE_", .rxToNonmemHandleNamesOrAtomic(str2lang(v), .ui))
                 }, character(1), USE.NAMES=TRUE)
  .ini <- vapply(names(.mv$ini[!is.na(.mv$ini)]),
                 function(v) {
                   paste0("RXE_", .rxToNonmemHandleNamesOrAtomic(str2lang(v), .ui))
                 }, character(1), USE.NAMES=TRUE)
  assign(".thetaMu", c(.lhs, .ini), envir=.ui)
  .nonmemResetUi(.ui, "E")
  #rxode2::rxAssignControlValue(.ui, ".nmVarExtra", "E")
  .err <- rxToNonmem(.normMain, .ui)
  .nonmemResetUi(.ui, "")
  #rxode2::rxAssignControlValue(.ui, ".nmVarExtra", "")
  rm(".thetaMu", envir=.ui)
  if (is.null(.advan)) {
    .norm <- rxode2::rxNorm(.mv)
    .des <- rxToNonmem(.norm, ui=.ui)
  } else {
    .des <- .nonmemLinCmtPk(.pkModel, .advan, .ui)
  }
  .prop <- .nonmemGetCmtProperties(.ui)
  .pk2 <- vapply(seq_along(.prop$cmt),
                 function(i) {
                   .cmt <- .prop$cmt[i]
                   .ret <- NULL
                   if (!is.na(.prop$f[i])) {
                     .ret <- c(.ret,
                               paste0("  F", .cmt, "=", .prop$f[i]))
                   }
                   if (!is.na(.prop$rate[i])) {
                     .ret <- c(.ret,
                               paste0("  R", .cmt, "=", .prop$rate[i]))
                   }
                   if (!is.na(.prop$dur[i])) {
                     .ret <- c(.ret,
                               paste0("  D", .cmt, "=", .prop$dur[i]))
                   }
                   if (!is.na(.prop$lag[i])) {
                     .ret <- c(.ret,
                               paste0("  ALAG", .cmt, "=", .prop$lag[i]))
                   }
                   if (!is.na(.prop$init[i])) {
                     .ret <- c(.ret,
                               paste0("  A_0(", .cmt, ")=", .prop$init[i]))
                   }
                   if (is.null(.ret)) return(NA_character_)
                   paste(.ret, collapse="\n")
                 }, character(1), USE.NAMES=FALSE)
  .pk2 <- .pk2[!is.na(.pk2)]
  .pk2 <- ifelse(length(.pk2) > 0, paste0("\n", paste(.pk2, collapse="\n")), "")
  if (!is.null(.advan)) {
    return(paste0(.pk, .des, .pk2,
                  "\n\n$ERROR\n  ;Redefine LHS in $PK by prefixing with on RXE_ for $ERROR\n",
                  .err))
  }
  paste0(.pk, .pk2,
         ifelse(.isPred, "\n", "\n\n$DES\n"),
         .des,
         ifelse(.isPred, "\n", "\n\n$ERROR\n  ;Redefine LHS in $DES by prefixing with on RXE_ for $ERROR\n"),
         .err)
}
attr(rxUiGet.nonmemPkDesErr0, "rstudio") <- "nonmemPkDesErr0"

#' NONMEM closed-form ADVAN for a linCmt() model
#'
#' @param ui rxode2 ui (the ODE version of the linCmt() model) with the
#'   micro-constants in the `.linCmtMicro` control value
#' @return `NULL` when the model is solved with ODEs, otherwise a list
#'   with `advan` (the ADVAN name), `par` (named list of TRANS1 NONMEM
#'   parameter to micro-constant expression) and `reserved` (TRANS1
#'   names model variables cannot use)
#' @noRd
#' @author Matthew L. Fidler
.nonmemLinCmtAdvan <- function(ui) {
  .micro <- rxode2::rxGetControl(ui, ".linCmtMicro", NULL)
  if (is.null(.micro)) return(NULL)
  .oral <- .micro$oral0 == 1L
  if (.micro$ncmt == 1L) {
    .advan <- ifelse(.oral, "ADVAN2", "ADVAN1")
    .par <- list(K=.micro$k)
  } else if (.micro$ncmt == 2L) {
    if (.oral) {
      .advan <- "ADVAN4"
      .par <- list(K=.micro$k, K23=.micro$k12, K32=.micro$k21)
    } else {
      .advan <- "ADVAN3"
      .par <- list(K=.micro$k, K12=.micro$k12, K21=.micro$k21)
    }
  } else {
    if (.oral) {
      .advan <- "ADVAN12"
      .par <- list(K=.micro$k, K23=.micro$k12, K32=.micro$k21,
                   K24=.micro$k13, K42=.micro$k31)
    } else {
      .advan <- "ADVAN11"
      .par <- list(K=.micro$k, K12=.micro$k12, K21=.micro$k21,
                   K13=.micro$k13, K31=.micro$k31)
    }
  }
  if (.oral) .par$KA <- .micro$ka
  # a model variable that is already the TRANS1 parameter (like ka ->
  # KA) needs no extra line; other model variables with these names
  # are renamed
  .same <- vapply(names(.par), function(n) {
    .e <- .par[[n]]
    is.name(.e) && identical(gsub(".", "_", toupper(as.character(.e)), fixed=TRUE), n)
  }, logical(1), USE.NAMES=FALSE)
  # NONMEM also knows the rate constants by compartment number (like
  # K12 for KA and K20 for K in ADVAN4), so a model variable cannot use
  # any of those names either
  .alias <- c("K", "KA", "K10", "K20", "K30", "K40", "K12", "K21", "K13",
              "K31", "K23", "K32", "K24", "K42", "K34", "K43")
  list(advan=.advan, par=.par[!.same],
       reserved=setdiff(.alias, names(.par)[.same]))
}

#' Write the $PK lines for a closed-form NONMEM ADVAN
#'
#' @param pkModel state independent model lines (calculated in $PK)
#'   and the ODEs (which are dropped)
#' @param advan list from `.nonmemLinCmtAdvan()`
#' @param ui rxode2 ui
#' @return NONMEM $PK lines (starting with a new line) with the
#'   model lines and the TRANS1 micro-constants
#' @noRd
#' @author Matthew L. Fidler
.nonmemLinCmtPk <- function(pkModel, advan, ui) {
  .ret <- ""
  if (length(pkModel) > 0L) {
    .mv <- rxode2::rxModelVars(paste(vapply(seq_along(pkModel),
                                            function(i) {
                                              deparse1(pkModel[[i]])
                                            }, character(1), USE.NAMES=FALSE),
                                     collapse="\n"))
    .norm <- rxode2::rxNorm(.mv)
    .pk <- strsplit(rxToNonmem(.norm, ui=ui), "\n")[[1]]
    .pk <- .pk[!grepl("^ *DADT[(]", .pk)]
    if (length(.pk) > 0L) .ret <- paste0("\n", paste(.pk, collapse="\n"))
  }
  if (length(advan$par) > 0L) {
    .ret <- paste0(.ret, "\n",
                   paste(vapply(names(advan$par), function(n) {
                     .e <- advan$par[[n]]
                     paste0("  ", n, "=", .rxToNonmem(.e, ui=ui),
                            .babelmixr2Deparse(.e))
                   }, character(1), USE.NAMES=FALSE),
                   collapse="\n"))
  }
  .ret
}
