.popedInfo <- new.env(parent=emptyenv())
#' get the bpop number (which is a theta in PopED)
#'
#' @param theta name of the population parameter
#' @param ui rxode2 ui object
#' @return bpop[#] where # is the theta number
#' @noRd
#' @author Matthew L. Fidler
.popedGetBpopNum <- function(theta, ui) {
  .iniDf <- ui$iniDf
  .w <- which(.iniDf$name == theta)
  if (length(.w) != 1) return(NA_character_)
  if (is.na(.iniDf$ntheta[.w])) return(NA_character_)
  paste0("bpop[", .iniDf$ntheta[.w], "]")
}

#' @export
rxUiGet.popedBpopRep <- function(x, ...) {
  .ui <- x[[1]]
  .split <- .ui$getSplitMuModel
  .muRef <- c(.split$pureMuRef, .split$taintMuRef)
  .thetas <- names(.muRef)
  .covRefDf <- .ui$saemMuRefCovariateDataFrame
  .ret <- data.frame(theta=.thetas,
                     bpop=vapply(.thetas, .popedGetBpopNum,
                                    character(1), ui=.ui, USE.NAMES=FALSE),
                     mu=vapply(.thetas, .nonmemGetMuNum, character(1), ui=.ui,
                               USE.NAMES=FALSE),
                     cov=vapply(.thetas, .nonmemGetThetaMuCov, character(1),
                                ui=.ui, covRefDf=.covRefDf, USE.NAMES=FALSE))
  .ret$b <- ifelse(is.na(.ret$mu), NA_character_,
                        paste0("b[",substr(.ret$mu,4, 10),"]"))
  .ret
}
attr(rxUiGet.popedBpopRep, "desc") <- "PopED data frame for replacements used in $popedFgFun"


#' @export
rxUiGet.popedFgFun <- function(x, ...) {
  # function(x, a, bpop, b, bocc)
  # x=?
  # a=covariates (could be dose, tau etc)
  # bpop = population variables
  # b = eta variables
  # bocc = occasion variables
  .ui <- x[[1]]
  .split <- .ui$getSplitMuModel

  .mu <- rxUiGet.popedBpopRep(x, ...)
  .ret <- vapply(seq_along(.mu$mu), function(i) {
    if (is.na(.mu$mu[i])) return(NA_character_)
    paste0(.mu$mu[i], " <- ", .mu$bpop[i],
           ifelse(is.na(.mu$cov[i]), "",
                  paste0("+", .mu$cov[i])))
  }, character(1), USE.NAMES=FALSE)
  .ret <- .ret[!is.na(.ret)]
  .mu2 <- ifelse(is.na(.mu$mu),
                          paste0(.mu$theta, "<-", .mu$bpop),
                 paste0(.mu$theta, "<-", .mu$mu, "+", .mu$b))
  .covDef <- .ui$allCovs
  .covDefLst <- lapply(seq_along(.covDef),
                       function(i) {
                         str2lang(paste0(.covDef[i], "<- a[", i + 1, "]"))
                       })
  .v <- c(.split$pureMuRef, .split$taintMuRef, .covDef)
  .allCovs <- .ui$allCovs
  .body1 <- c(list(quote(`{`)),
              lapply(c(.ret, .mu2),
                     function(x) {
                       str2lang(x)
                     }),
              .covDefLst,
              .split$muRefDef,
              list(str2lang(paste("c(ID=a[1],", paste(paste0(.v, "=", .v), collapse=","),
                                  ")"))))
  .body1 <- as.call(.body1)
  .f <- function(x, a, bpop, b, bocc) {}
  body(.f) <- .body1
  .f
}
attr(rxUiGet.popedFgFun, "desc") <- "PopED parameter model (fg_fun)"

.poped <- new.env(parent=emptyenv())
.poped$s <- NULL
.poped$minVal <- NULL

#' Assign the rxode2 solved value to .poped$s for calling in the error function
#'
#' This should not be called directly
#'
#' @param s rxode2 solved dataset
#' @return nothing, called for side effects
#' @export
#' @keywords internal
#' @author Matthew L. Fidler
.popedA <- function(s) {
  .poped$s <- s
}
#' Get the weight from the rxode2 solve
#'
#' This shouldn't be called directly
#'
#' @return rxode2 weights for the poped error function
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.popedW <- function() {
  sqrt(.poped$s[["rx_r_"]])
}
#' Get the function value from the rxode2 solve
#'
#' This shouldn't be called directly
#'
#' @return rxode2 weights for the poped error function
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.popedF <- function() {
  .poped$s[["rx_pred_"]]
}

#'@export
rxUiGet.popedFErrorFun  <- function(x, ...) {
  function(model_switch, xt, parameters, epsi, poped.db) {
    # Will assign in babelmixr2 namespace
    do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))
    y <- .popedF() + .popedW() * eps[, 1]
    return(list(y=y,poped.db=poped.db))
  }
}
attr(rxUiGet.popedFErrorFun, "desc") <- "PopED error function(fError_fun)"

#' @export
rxUiGet.popedRxmodelBase <- function(x, ...) {
  .ui <- x[[1]]
  .ui <- rxode2::rxUiDecompress(.ui)
  .split <- .ui$getSplitMuModel
  .model2 <- .split$modelWithDrop
  .w <- which(vapply(seq_along(.model2),
                     function(i) {
                       !identical(.model2[[i]], quote(`_drop`))
                     }, logical(1), USE.NAMES=FALSE))
  .model2 <- .model2[.w]
  .model2 <- vapply(seq_along(.model2),
                    function(i) {
                      deparse1(.model2[[i]])
                    }, character(1), USE.NAMES = FALSE)
  .ui2 <- rxode2::rxUiDecompress(.ui)
  # test edge cases of all thetas and no thetas
  .theta <- .ui$iniDf[!is.na(.ui$iniDf$ntheta), "name"]
  .eta <- .ui$iniDf[which(.ui$iniDf$neta1 == .ui$iniDf$neta2), "name"]
  .pars <- c(.theta, .eta)
  .pars <- paste0("rxDummy", .pars, " <- ", .pars)
  .model2 <- c(.pars,
               .model2)
  suppressMessages(suppressWarnings(model(.ui2) <- .model2))
  .ui2 <- rxode2::rxUiDecompress(.ui2)
  .mod <- rxode2::rxCombineErrorLines(.ui2, errLines=nlmixr2est::rxGetDistributionFoceiLines(.ui2),
                                      cmtLines=FALSE, dvidLine=FALSE, useIf = FALSE)
  ## .mod <- rxode2::rxCombineErrorLines(.ui2, errLines=popedErr(.ui2),
  ##                                     cmtLines=FALSE, dvidLine=FALSE, useIf = FALSE)
  .mod <- .mod[[2]]
  .mod <- as.list(.mod[seq_along(.mod)[-1]])
  # drop dummy and params()
  .w <- which(vapply(seq_along(.mod),
                     function(i) {
                       if (i == 1) return(FALSE)
                       .cur <- .mod[[i]]
                       if (identical(.cur[[1]], quote(`<-`))) {
                         .lhs <- deparse1(.cur[[2]])
                         if (substr(.lhs, 1, 7) == "rxDummy") {
                           return(FALSE)
                         }
                       }
                       TRUE
                     },
                     logical(1), USE.NAMES = FALSE))
  .mod <- .mod[.w]
  .mod <- lapply(seq_along(.mod),
                 function(i) {
                   .cur <- .mod[[i]]
                   if (identical(.cur[[1]], quote(`<-`)) ||
                         identical(.cur[[1]], quote(`=`))) {
                     .cur[[1]] <- quote(`~`)
                   }
                   if (identical(.cur[[1]], quote(`~`))) {
                     if (identical(.cur[[2]], quote(`rx_pred_`)) ||
                           identical(.cur[[2]], quote(`rx_r_`)))
                       .cur[[1]] <- quote(`<-`)
                   }
                   .cur
                 })
  .mod
}
attr(rxUiGet.popedRxmodelBase, "desc") <- "This gets the base rxode2 model for PopED"
#' Get PopEd's rxode2 model (unevaluated)
#'
#' This is an internal function and really shouldn't be called directly
#'
#' @param ui User interface
#'
#' @param maxNumTime Maximum number of evaluation times used in model
#'   (usually these are design points).  They are implemented as
#'   modeling times.  Note when there are modeling times present, you
#'   cannot use the faster solving method (though slow solving still
#'   works)
#'
#' @return unevaluated rxode2 model with 2 attributes:
#'
#'  - First attribute is "mtime", which tells if the modeled times
#'    represent design points (boolean)
#'
#'  - Second attribute is "maxNumTime", which tells the maximum number
#'    of design points present in this model (integerish)
#'
#' @noRd
#' @keywords internal
#' @author Matthew L. Fidler
.popedRxModel <- function(ui, maxNumTime=2) {
  checkmate::testIntegerish(maxNumTime, lower=1, len=1)
  .base <- rxUiGet.popedRxmodelBase(list(ui))
  .base2 <- eval(as.call(c(list(quote(`rxModelVars`)),
                           as.call(c(list(quote(`{`)), .base)))))
  if (.base2$nMtime > 0L) {
    .minfo("mtime() detected in model, slower solving used")
    .ret <- as.call(c(list(quote(`rxode2`)),
                        as.call(c(list(quote(`{`)), .base))))
    attr(.ret, "mtime") <- FALSE
    attr(.ret, "maxNumTime") <- 0L
    return(.ret)
  }
  .split <- ui$getSplitMuModel
  .covDef <- ui$allCovs
  .param <- str2lang(paste0("params(", paste(c(.split$pureMuRef, .split$taintMuRef, .covDef), collapse=","), ")"))
  .ret <- as.call(c(list(quote(`rxode2`)),
                    as.call(c(list(quote(`{`), .param),
                              .base,
                              lapply(seq_len(maxNumTime),
                                     function(i) {
                                       str2lang(sprintf("mtime(rxXt_%d_v) ~ rxXt_%d", i, i))
                                     })))))
  attr(.ret, "mtime") <- TRUE
  attr(.ret, "maxNumTime") <- maxNumTime
  .ret
}



#' @export
rxUiGet.popedFfFun <- function(x, ...) {
  .ui <- x[[1]]

}

#' @export
rxUiGet.popedBpop <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .bpop <- .iniDf[!is.na(.ui$iniDf$ntheta), ]
  setNames(.bpop$est, .bpop$name)
}
attr(rxUiGet.popedBpop, "desc") <- "Get PopED's $bpop"

#' @export
rxUiGet.popedNotfixedBpop <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .bpop <- .iniDf[!is.na(.ui$iniDf$ntheta), ]
  1 - .bpop$fix * 1
}
attr(rxUiGet.popedBpop, "desc") <- "Get PopED's $notfixed_bpop"

## Omega matrix specification
#' @export
rxUiGet.popedD <- function(x, ...) {
  .ui <- x[[1]]
  .omega <- .ui$omega
  # does it have to have the mu-referenced named (just like examples?)
  diag(.omega)
}
attr(rxUiGet.popedD, "desc") <- "Get PopED's $d parameter"

#' @export
rxUiGet.popedCovd <- function(x, ...){
  .ui <- x[[1]]
  .omega <- .ui$omega
  .ret <- .omega[lower.tri(.omega)]
  names(.ret) <- NULL
  .ret
}
attr(rxUiGet.popedCovd, "desc") <- "Get PopED's $covd parameter"

#' @export
rxUiGet.popedNotfixedD <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .iniDf <- .iniDf[which(is.na(.iniDf$ntheta) & .iniDf$neta1 == .iniDf$neta2), ]
  1 - .iniDf$fix * 1
}
attr(rxUiGet.popedNotfixedD, "desc") <- "Get PopED's $notfixed_d"

#' @export
rxUiGet.popedNotfixedCovd <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .lotri <- lotri::as.lotri(.iniDf)
  .lotriFix <- attr(.lotri,"lotriFix")
  if (is.null(.lotriFix)) {
    .ret <- .lotri[lower.tri(.lotri)]
    return(rep(1, length(.ret)))
  }
  .ret <- .lotriFix[lower.tri(.lotriFix)]
  names(.ret) <- NULL
  1 - .ret * 1
}
attr(rxUiGet.popedNotfixedCovd, "desc") <- "Get PopED's $notfixed_covd"

## Currently nlmixr2 only has sigma fixed to 1.
#' @export
rxUiGet.popedSigma <- function(x, ...) {
  return(1.0)
}
attr(rxUiGet.popedSigma, "desc") <- "PopED database $sigma"

#' @export
rxUiGet.popedNotfixedSigma <- function(x, ...) {
  return(0L)
}
attr(rxUiGet.popedNotfixedSigma, "desc") <- "PopED database $notfixed_sigma"

#' Create a babelmixr2/nlmixr2 design space based on a data frame
#'
#' From the dataset construct `xt`, `minxt`, `maxxt` and `a`
#'
#' This does not work with time-varying covariates
#'
#' @param ui rxode2 ui object
#' @param data babelmixr2 design data, very similar to a rxode2 modeling dataset
#' @param time string that represents the time in the dataset (ie xt)
#' @param timeLow string that represents the lower design time (ie minxt)
#' @param timeHi string that represents the upper design time (ie maxmt)
#' @param id The id variable
#' @inheritParams PopED::create_design
#' @inheritParams PopED::create_design_space
#' @return PopED design space
#' @export
#' @author Matthew L. Fidler
popedDataToDesignSpace <- function(ui, data, groupsize, time="time", timeLow="low", timeHi="high",
                                   id="id", m = NULL, x = NULL, ni = NULL,
                                   model_switch = NULL,
                                   maxni = NULL,
                                   minni = NULL,
                                   maxtotni = NULL,
                                   mintotni = NULL,
                                   maxgroupsize = NULL,
                                   mingroupsize = NULL,
                                   maxtotgroupsize = NULL,
                                   mintotgroupsize = NULL,
                                   xt_space = NULL,
                                   maxa = NULL,
                                   mina = NULL,
                                   a_space = NULL,
                                   x_space = NULL,
                                   use_grouped_xt = FALSE,
                                   grouped_xt = NULL,
                                   use_grouped_a = FALSE,
                                   grouped_a = NULL,
                                   use_grouped_x = FALSE,
                                   grouped_x = NULL,
                                   our_zero = NULL) {
  rxode2::rxReq("PopED")
  ui <- rxode2::assertRxUi(ui)
  .et <- rxode2::etTrans(data, ui)
  .tmp <- attr(class(.et),".rxode2.lst")
  class(.tmp) <- NULL
  .a <- as.matrix(.tmp$cov1)
  .a <- .a[, c("ID", ui$allCovs)]
  .nd <- tolower(names(data))
  .data <- data
  .wevid <- which(.nd == "evid")
  if (length(.wevid) == 0L) {
    .minfo("could not find evid, assuming all are design points")
    .data$evid <- 0
  }
  .wtime <- which(.nd == time)
  if (length(.wtime) != 1L) {
    stop(paste0("could not find the required data item: ", time),
         call.=FALSE)
  }
  .wid <- which(.nd == id)
  if (length(.wid) != 1L) {
    stop(paste0("could not find the required data item: ", id),
         call.=FALSE)
  }
  .env <- new.env(parent=emptyenv())
  .env$mt <- -Inf
  .design <- PopED::create_design(xt=lapply(unique(.data[[.wid]]),
                                            function(id) {
                                              .data <- .data[.data[[.wid]] == id &
                                                               .data[[.wevid]] == 0, ]
                                              .ret <- .data[[.wtime]]
                                              .env$mt <- max(c(.ret, .env$mt))
                                              .ret
                                            }),
                                  groupsize=groupsize,
                                  m = m, x = x, a = .a, ni = ni,
                                  model_switch = model_switch)
  .wlow <- which(.nd == timeLow)
  .minxt <- NULL
  if (length(.wlow) == 1L) {
    .minxt <- lapply(unique(.data[[.wid]]),
                     function(id) {
                       .data <- .data[.data[[.wid]] == id &
                                        .data[[.wevid]] == 0, ]
                       .ret <- .data[[.wlow]]
                       .env$mt <- max(c(.ret, .env$mt))
                       .ret
                     })
  }
  .whi <- which(.nd == timeHi)
  .maxxt <- NULL
  if (length(.whi) == 1L) {
    .maxxt <- lapply(unique(.data[[.wid]]),
                     function(id) {
                       .data <- .data[.data[[.wid]] == id &
                                        .data[[.wevid]] == 0, ]
                       .ret <- .data[[.whi]]
                       .env$mt <- max(c(.ret, .env$mt))
                       .ret
                     })
  }
  .designSpace <-
    PopED::create_design_space(.design,
                               maxni = maxni,
                               minni = minni,
                               maxtotni = maxtotni,
                               mintotni = mintotni,
                               maxgroupsize = maxgroupsize,
                               mingroupsize = mingroupsize,
                               maxtotgroupsize = maxtotgroupsize,
                               mintotgroupsize = mintotgroupsize,
                               maxxt=.maxxt,
                               minxt=.minxt,
                               xt_space = xt_space,
                               mina = mina,
                               maxa = maxa,
                               a_space = a_space,
                               x_space = x_space,
                               use_grouped_xt = use_grouped_xt,
                               grouped_xt = grouped_xt,
                               use_grouped_a = use_grouped_a,
                               grouped_a = grouped_a,
                               use_grouped_x = use_grouped_x,
                               grouped_x = grouped_x,
                               our_zero = our_zero)
  .rx <- .popedRxModel(ui, maxNumTime=max(.designSpace$design$ni))
  if (!attr(.rx, "mtime")) {
    stop("mtime models are not supported yet",
         call.=FALSE)
  }
  .rx <- eval(.rx)
  .popEd$model <- .rx
  .wamt <- which(.nd == "amt")
  .wrate <- which(.nd == "rate")
  .wdur <- which(.nd == "dur")
  .wss <- which(.nd == "ss")
  .wii <- which(.nd == "ii")
  .waddl <- which(.nd == "addl")
  # Create a dataset without these design points with one observation
  # 0.5 units after
  .dat <- do.call(rbind,
                  lapply(unique(.data[[.wid]]),
                         function(id) {
                           .data <- .data[.data[[.wid]] == id &
                                            .data[[.wevid]] != 0, ]
                           .len <- length(.data[[.wid]])
                           .data2 <- .data[.len, ]
                           .data2[[.wtime]] <- .env$mt + 0.5
                           .data2[[.wevid]] <- 0
                           if (length(.wamt) == 1L) .data2[[.wamt]] <- NA
                           if (length(.wrate) == 1L) .data2[[.wrate]] <- NA
                           if (length(.wdur) == 1L) .data2[[.wdur]] <- NA
                           if (length(.wss) == 1L) .data2[[.wss]] <- NA
                           if (length(.wii) == 1L) .data2[[.wii]] <- NA
                           if (length(.waddl) == 1L) .data2[[.waddl]] <- NA
                           rbind(.data, .data2)
                         }))
  .popEd$data <- .dat
  #nlmixr2est::.popedSetup(.popEd)
  .designSpace
}

# modelswitch --perhaps for defining endpoint in nlmixr2 model
# If so, defined by the CMT and could be inferred from the model

## PopED control needs to include:
## - optsw
# iFIMCalculationType
# iApproximationMethod
# prior_fim
# d_switch
# ofv_calc_type
# ds_index


## - groupsize -- size of the different groups
## - maxgroupsize -- maximum size of the different groups
## - mingroupsize -- minimum group size of the different groups
## - maxtotgroupsize -- maximum total group size
## - mintotgroupsize -- minimum total group size

## Should we have something line:
## ni = list(c(min, sample, max),
##           )
##
## - ni -- vector defining number of samples in each group
## - maxni
## - minni
## - maxtotni -- maximum number of samples in the experiment
## - mintotni -- minimum number of samples in the experiment





## Things that can come from data (like NONMEM's $DESIGN or from pop ed control:
## - xt -- Matrix defining the initial sampling schedule
## - minxt -- Matrix defining the minimum sampling time
## - maxxt -- Matrix defining the maximum sampling time
## - m  -- number of groups in the study
## - x  -- A matrix defining the initial discrete values for the model Each row is a group/individual.
## - nx -- number of discrete values for each row/column
## - a -- continuous covariate values; not sure if taken care of elsewhere...?
#' Control for a PopED design task
#'
#' @inheritParams nlmixr2est::foceiControl
#' @param ... other parameters for PopED design task
#' @return popedControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
popedControl <- function(stickyRecalcN=4,
                         maxOdeRecalc=5,
                         odeRecalcFactor=10^(0.5),
                         rxControl=NULL,
                         sigdig=4,
                         ...) {
  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% c("genRxControl"))]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste
    (paste0("'", .bad, "'", sep=""), collapse=", "),
    call.=FALSE)
  }

  checkmate::assertIntegerish(stickyRecalcN, any.missing=FALSE, lower=0, len=1)
  checkmate::assertIntegerish(maxOdeRecalc, any.missing=FALSE, len=1)
  checkmate::assertNumeric(odeRecalcFactor, len=1, lower=1, any.missing=FALSE)

  .genRxControl <- FALSE
  if (!is.null(.xtra$genRxControl)) {
    .genRxControl <- .xtra$genRxControl
  }
  if (is.null(rxControl)) {
    if (!is.null(sigdig)) {
      rxControl <- rxode2::rxControl(sigdig=sigdig)
    } else {
      rxControl <- rxode2::rxControl(atol=1e-4, rtol=1e-4)
    }
    .genRxControl <- TRUE
  } else if (inherits(rxControl, "rxControl")) {
  } else if (is.list(rxControl)) {
    rxControl <- do.call(rxode2::rxControl, rxControl)
  } else {
    stop("solving options 'rxControl' needs to be generated from 'rxode2::rxControl'", call=FALSE)
  }
  if (!is.null(sigdig)) {
    checkmate::assertNumeric(sigdig, lower=1, finite=TRUE, any.missing=TRUE, len=1)
  }
  .ret <- list(rxControl=rxControl,
               stickyRecalcN=as.integer(stickyRecalcN),
               maxOdeRecalc=as.integer(maxOdeRecalc),
               odeRecalcFactor=odeRecalcFactor,
               sigdig=sigdig,
               genRxControl=.genRxControl)
  class(.ret) <- "popedControl"
  .ret
}
