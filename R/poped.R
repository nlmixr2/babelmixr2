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
  .param <- str2lang(paste0("params(", paste(c(.split$pureMuRef, .split$taintMuRef, .covDef,
                                               paste0("rxXt_", seq_len(maxNumTime))),
                                             collapse=","), ")"))
  .poped$param <- c(.split$pureMuRef, .split$taintMuRef, .covDef)
  .poped$param <- setNames(rep(1, length(.poped$param)), .poped$param)
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
attr(rxUiGet.popedNotfixedBpop, "desc") <- "Get PopED's $notfixed_bpop"

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
rxUiGet.popedCovd <- function(x, ...) {
  .ui <- x[[1]]
  .omega <- .ui$omega
  .ret <- .omega[lower.tri(.omega)]
  names(.ret) <- NULL
  if (all(.ret == 0)) return(NULL)
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
#' @noRd
#' @author Matthew L. Fidler
.popedDataToDesignSpace <- function(ui, data, groupsize=NULL, time="time", timeLow="low", timeHi="high",
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
  groupsize <- rxode2::rxGetControl(ui, "groupsize", groupsize)
  if (is.null(groupsize)) {
    .minfo("groupsize should be specified; but for now assuming 20")
    groupsize <- 20
  }

  # Get from assigned control inside of ui object
  for (opt in c("m", "x"," ni", "model_switch", "maxni","minni", "maxtotni",
                "mintotni", "maxgroupsize","mingroupsize","maxtotgroupsize",  "mintotgroupsize",
                "xt_space", "maxa", "mina", "a_space", "x_space", "grouped_xt", "use_grouped_a",
                "grouped_a", "grouped_x", "our_zero", "use_grouped_xt",
                "use_grouped_a", "use_grouped_x")) {
    assign(opt, rxode2::rxGetControl(ui, opt, get(opt)))
  }
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
  .poped$model <- .rx
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
  .poped$data <- .dat
  #nlmixr2est::.popedSetup(.popEd)
  .designSpace
}

# modelswitch --perhaps for defining endpoint in nlmixr2 model
# If so, defined by the CMT and could be inferred from the model

## PopED control needs to include:
## Settings:
## - optsw

#' @export
rxUiGet.popedSettings <- function(x, ...) {
  ui <- x[[1]]
  rxReq("PopED")
  .line_opta <- rxode2::rxGetControl(ui, "line_opta", NULL)
  if (checkmate::testLogical(.line_opta, any.missing = FALSE, len=1)) {
    .line_opta <- .line_opta * 1
  }
  .line_optx <- rxode2::rxGetControl(ui, "line_optx", NULL)
  if (checkmate::testLogical(.line_optx, any.missing = FALSE, len=1)) {
    .line_optx <- .line_optx * 1
  }
  list(settings=list(
    optsw=c(rxode2::rxGetControl(ui, "optSample", FALSE) * 1,
            rxode2::rxGetControl(ui, "optDiscreteDesignVar", FALSE) * 1,
            rxode2::rxGetControl(ui, "optContDesignVar", FALSE) * 1,
            rxode2::rxGetControl(ui, "optIdPerGroup", FALSE) * 1),
    iFIMCalculationType=rxode2::rxGetControl(ui, "iFIMCalculationType", 1),
    iApproximationMethod=rxode2::rxGetControl(ui, "iApproximationMethod", 0),
    iFOCENumInd=rxode2::rxGetControl(ui, "iFOCENumInd", 1000),
    prior_fim=rxode2::rxGetControl(ui, "prior_fim", matrix(0, 0, 1)),
    ofv_calc_type=rxode2::rxGetControl(ui, "ofv_calc_type", 4),
    strEDPenaltyFile=rxode2::rxGetControl(ui, "strEDPenaltyFile", ""),
    ofv_fun=rxode2::rxGetControl(ui, "ofv_fun", NULL),
    iEDCalculationType=rxode2::rxGetControl(ui, "iEDCalculationType", 0),
    ED_samp_size=rxode2::rxGetControl(ui, "ED_samp_size", 45),
    bLHS=rxode2::rxGetControl(ui, "bLHS", 1),
    bUseRandomSearch=rxode2::rxGetControl(ui, "bUseRandomSearch", TRUE) * 1,
    bUseStochasticGradient=rxode2::rxGetControl(ui, "bUseStochasticGradient", TRUE) * 1,
    bUseLineSearch=rxode2::rxGetControl(ui, "bUseLineSearch", TRUE) * 1,
    bUseExchangeAlgorithm=rxode2::rxGetControl(ui, "bUseExchangeAlgorithm", FALSE) * 1,
    bUseBFGSMinimizer=rxode2::rxGetControl(ui, "bUseBFGSMinimizer", FALSE) * 1,
    EACriteria=rxode2::rxGetControl(ui, "EACriteria", 1),
    run_file_pointer=rxode2::rxGetControl(ui, "strRunFile", ""),
    poped_version=rxode2::rxGetControl(ui, "poped_version", packageVersion("PopED")),
    modtit=rxode2::rxGetControl(ui, "modtit", "PopED babelmixr2 model"),
    output_file=rxode2::rxGetControl(ui, "output_file", "PopED_output_summary"),
    output_function_file=rxode2::rxGetControl(ui, "output_function_file", "PopED_output_"),
    strIterationFileName=rxode2::rxGetControl(ui, "strIterationFileName", "PopED_current.R"),
    user_data=rxode2::rxGetControl(ui, "user_data", PopED::cell(0, 0)),
    ourzero=rxode2::rxGetControl(ui, "ourzero", 1e-5),
    dSeed=rxode2::rxGetControl(ui, "dSeed", NULL),
    line_opta=.line_opta,
    line_optx=.line_optx,
    bShowGraphs=rxode2::rxGetControl(ui, "bShowGraphs", FALSE),
    use_logfile=rxode2::rxGetControl(ui, "use_logfile", FALSE),
    m1_switch=rxode2::rxGetControl(ui, "m1_switch", 1),
    m2_switch=rxode2::rxGetControl(ui, "m2_switch", 2),
    hle_switch=rxode2::rxGetControl(ui, "hle_switch", 1),
    gradff_switch=rxode2::rxGetControl(ui, "gradff_switch", 1),
    gradfg_switch=rxode2::rxGetControl(ui, "gradfg_switch", 1),
    grad_all_switch=rxode2::rxGetControl(ui, "grad_all_switch", 1),
    rsit_output=rxode2:rxGetControl(ui, "rsit_output", 5),
    sgit_output=rxode2::rxGetControl(ui, "sgit_output", 1),
    hm1=rxode2::rxGetControl(ui, "hm1", 1e-5),
    hlf=rxode2::rxGetControl(ui, "hlf", 1e-5),
    hlg=rxode2::rxGetControl(ui, "hlg", 1e-5),
    hm2=rxode2::rxGetControl(ui, "hm2", 1e-5),
    hgd=rxode2::rxGetControl(ui, "hgd", 1e-5),
    hle=rxode2::rxGetControl(ui, "hle", 1e-5),
    AbsTol=Rxode2::RxGetControl(Ui, "AbsTol", 1e-5),
    RelTol=rxode2::rxGetControl(ui, "RelTol", 1e-5),
    iDiffSolverMethod=rxode2::rxGetControl(ui, "iDiffSolverMethod", NULL),
    bUseMemorySolver=rxode2::rxGetControl(ui, "bUseMemorySolver", FALSE),
    rsit=rxode2::rxGetControl(ui, "rsit", 300),
    sgit=rxode2::rxGetControl(ui, "sgit", 150),
    intrsit=rxode2::rxGetControl(ui, "trsit", 250),
    intsgit=rxode2::rxGetControl(ui, "tsgit", 50),
    maxrsnullit=rxode2::rxGetControl(ui, "ullit", 50),
    convergence_eps=rxode2::rxGetControl(ui, "ce_eps", 1e-08),
    rslxt=rxode2::rxGetControl(ui, "rslxt", 10),
    rsla=rxode2::rxGetControl(ui, "rsla", 10),
    cfaxt=rxode2::rxGetControl(ui, "cfaxt", 0.001),
    cfaa=rxode2::rxGetControl(ui, "cfaa", 0.001),
    bGreedyGroupOpt=rxode2::rxGetControl(ui, "bGreedyGroupOpt", FALSE),
    EAStepSize=rxode2::rxGetControl(ui, "EAStepSize", 0.01),
    EANumPoints=rxode2::rxGetControl(ui, "EANumPoints", FALSE),
    EAConvergenceCriteria=rxode2::rxGetControl(ui, "EAConvergenceCriteria", 1e-20),
    bEANoReplicates=rxode2::rxGetControl(ui, "oReplicates", FALSE),
    BFGSProjectedGradientTol=rxode2::rxGetControl(ui, "BFGSProjectedGradientTol", 1e-04),
    BFGSTolerancef=rxode2::rxGetControl(ui, "BFGSTolerancef", 0.001),
    BFGSToleranceg=rxode2::rxGetControl(ui, "BFGSToleranceg", 0.9),
    BFGSTolerancex=rxode2::rxGetControl(ui, "BFGSTolerancex", 0.1),
    ED_diff_it=rxode2::rxGetControl(ui, "ED_diff_it", 30),
    ED_diff_percent=rxode2::rxGetControl(ui, "ED_diff_percent", 10),
    line_search_it=rxode2::rxGetControl(ui, "line_search_it", 50),
    Doptim_iter=rxode2::rxGetControl(ui, "Doptim_iter", 1),
    parallel=list(
      iCompileOption = rxode2::rxGetControl(ui, "iCompileOption", -1),
      iUseParallelMethod=rxode2::rxGetControl(ui, "iUseParallelMethod", 1),
      strExecuteName=rxode2::rxGetControl(ui, "strExecuteName", "calc_fim.exe"),
      iNumProcesses=rxode2::rxGetControl(ui, "iNumProcesses", 2),
      iNumChunkDesignEvals=rxode2::rxGetControl(ui, "iNumChunkDesignEvals", -2),
      Mat_Out_Pre=rxode2::rxGetControl(ui, "Mat_Out_Pre", "parallel_output"),
      strExtraRunOptions=rxode2::rxGetControl(ui, "strExtraRunOptions", ""),
      dPollResultTime=rxode2::rxGetControl(ui, "dPollResultTime", 0.1),
      strFunctionInputName=rxode2::rxGetControl(ui, "strFunctionInputName", "function_input"),
      bParallelRS=rxode2::rxGetControl(ui, "bParallelRS", FALSE),
      bParallelSG=rxode2::rxGetControl(ui, "bParallelSG", FALSE),
      bParallelMFEA=rxode2::rxGetControl(ui, "bParallelMFEA", FALSE),
      bParallelLS=rxode2::rxGetControl(ui, "bParallelLS", FALSE)
      )
  ))
}


# iFIMCalculationType
# iApproximationMethod
# prior_fimr
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

################################################################################
# $parameter

#' @export
rxUiGet.popedParameters <- function(x, ...) {
  ui <- x[[1]]
  .ret <- list(parameters=list(
    bpop=rxUiGet.popedBpop(x, ...),
    notfixed_bpop=rxUiGet.popedNotfixedBpop(x, ...),

    d=rxUiGet.popedD(x, ...),
    notfixed_d=rxUiGet.popedNotfixedD(x, ...),

    covd=rxUiGet.popedCovd(x, ...),
    notfixed_covd=rxUiGet.popedNotfixedCovd(x, ...),

    sigma=rxUiGet.popedSigma(x, ...),
    notfixed_sigma=rxUiGet.popedNotfixedSigma(x, ...)
    # no diagonals from nlmixr2 here:
    ## covsigma=rxUiGet.poped(x, ...),
    ## notfixed_covsigma=rxUiGet.popedNotfixedSigma(x, ...),
  ))
  if (rxode2::rxGetControl(ui, "ofv_calc_type", 4) == 6) {
    .ret <- .popedImportant(ui, .ret)
  }
  .ret
}
attr(rxUiGet.popedParameters, "desc") <- "PopED input $parameters"
#' Update poped$parameters with ds_index
#'
#' @param ui rxode2 ui (to change importance of standard deviation parameters/lambda etc)
#' @param lst Curent popedInput list
#' @param important character vector of important parameters or NULL for default
#' @param unimportant character vector of unimportant parameters or NULL for default
#' @return updated $parameters with ds_index added
#' @noRd
#' @author Matthew L. Fidler
.popedImportant <- function(ui, lst, important=NULL, unimportant=NULL) {
  important <- rxode2::rxGetControl(ui, "important", important)
  important <- rxode2::rxGetControl(ui, "unimportant", unimportant)
  .par <- lst$parameters
  .err <- ui$iniDf$name[!is.na(ui$iniDf$err)]
  .interstingFixed <- names(.par$bpop[which(.par$notfixed_bpop==1)])
  .uninterstingRandom <- names(.par$d[which(.par$notfixed_d==1)])
  # defaults
  #  1 if a parameter is uninteresting, otherwise 0
  .ds <- setNames(c(rep(0, length(.interstingFixed)),
                    rep(1, length(.uninterstingRandom))),
                  c(.interstingFixed, .uninterstingRandom))
  # should make sd uninteresting
  .ds[.err] <- 1
  important <- important[important %in% names(.ds)]
  unimportant <- unimportant[unimportant %in% names(.ds)]
  if (!is.null(important)) {
    .ds[important] <- 0
  }
  if (!is.null(unimportant)) {
    .ds[unimportant] <- 1
  }
  .minfo(paste0("Ds-optimality important parameters: ", paste(names(.ds)[which(.ds == 0)], collapse=", ")))
  .minfo(paste0("Ds-optimality unimportant parameters: ", paste(names(.ds)[which(.ds == 1)], collapse=", ")))
  .ret <- lst
  .ret$parameters$ds_index <- .ds
  .ret
}

#' Control for a PopED design task
#'
#' @param important character vector of important parameters or NULL
#'   for default.  This is used with Ds-optimality
#'
#' @param unimportant character vector of unimportant parameters or
#'   NULL for default.  This is used with Ds-optimality
#'
#' @param iFIMCalculationType can be either an integer or a named
#'   value of the Fisher Information Matrix type:
#'
#' - 0/"full" = Full FIM
#'
#' - 1/"reduced" = Reduced FIM
#'
#' - 2/"weighted" = weighted models
#'
#' - 3/"loc" = Loc models
#'
#' - 4/"reducedPFIM" = reduced FIM with derivative of SD of sigma as in PFIM
#'
#' - 5/"fullABC" = FULL FIM parameterized with A,B,C matrices & derivative of variance
#'
#' - 6/"largeMat" = Calculate one model switch at a time, good for large matrices
#'
#' - 7/"reducedFIMABC" = =Reduced FIM parameterized with A,B,C matrices & derivative of variance
#'
#' @param optSample logical; Optimize Samples per subject as well as FIM with current design
#'
#' @param optDiscreteDesignVar logical; optimize discrete design variable
#'
#' @param optContDesignVar logical; optimize continuous design variable
#'
#' @param optIdPerGroup logical; optimize individual per group
#'
#' @param iFOCENumInd integer; number of individuals in focei solve
#'
#' @param prior_fim matrix; prior FIM
#'
#' @param d_switch integer or character option:
#'
#' - 0/"ed" = ED design
#'
#' - 1/"d" = D design
#'
#' @param ofv_calc_type objective calculation type:
#'
#' - 1/"d" = D-optimality". Determinant of the FIM: det(FIM)
#'
#' - 2/"a" =  "A-optimality". Inverse of the sum of the expected parameter variances: 1/trace_matrix(inv(FIM))
#'
#' - 4/"lnD" = "lnD-optimality". Natural logarithm of the determinant of the FIM: log(det(FIM))
#'
#' - 6/"Ds" = "Ds-optimality". Ratio of the Determinant of the FIM and the Determinant of the uninteresting rows and columns of the FIM: det(FIM)/det(FIM_u)
#'
#' - 7/"inverse" = Inverse of the sum of the expected parameter RSE: 1/sum(get_rse(FIM,poped.db,use_percent=FALSE))
#'
#' @param iEDCalculationType ED Integral Calculation type:
#'
#' - 0/"mc" = Monte-Carlo-Integration
#'
#' - 1/"laplace" = Laplace Approximation
#'
#' - 2/"bfgs-laplace" = BFGS Laplace Approximation
#'
#' @param EACriteria Exchange Algorithm Criteria:
#'
#' - 1/"modified" = Modified
#'
#' - 2/"fedorov"  = Fedorov
#'
#' @param m1_switch Method used to calculate M1:
#'
#' - 1/"central" = Central difference
#'
#' - 0/"complex" = Complex difference
#'
#' - 20/"analytic" = Analytic derivative
#'
#' - 30/"ad" = Automatic differentiation
#'
#' @param m2_switch  Method used to calculate M2:
#'
#' - 1/"central" = Central difference
#'
#' - 0/"complex" = Complex difference
#'
#' - 20/"analytic" = Analytic derivative
#'
#' - 30/"ad" = Automatic differentiation
#'
#' @param hle_switch Method used to calculate linearization of residual error:
#'
#' - 1/"central" = Central difference
#'
#' - 0/"complex" = Complex difference
#'
#' - 30/"ad" = Automatic differentiation
#'
#' @param gradff_switch Method used to calculate the gradient of the model:
#'
#' - 1/"central" = Central difference
#'
#' - 0/"complex" = Complex difference
#'
#' - 20/"analytic" = Analytic derivative
#'
#' - 30/"ad" = Automatic differentiation
#'
#' @param gradfg_switch Method used to calculate the gradient of the
#'   parameter vector g:
#'
#' - 1/"central" = Central difference
#'
#' - 0/"complex" = Complex difference
#'
#' - 20/"analytic" = Analytic derivative
#'
#' - 30/"ad" = Automatic differentiation
#'
#' @param grad_all_switch Method used to calculate all the gradients:
#'
#' - 1/"central" = Central difference
#'
#' - 0/"complex" = Complex difference
#'
#' @param iCompileOption Compile options for PopED
#'
#'  - "none"/-1 = No compilation
#'
#'  - "full/0 or 3 = Full compilation
#'
#'  - "mcc"/1 or 4 = Only using MCC (shared lib)
#'
#'  - "mpi"/2 or 5 = Only MPI,
#'
#' When using numbers, option 0,1,2 runs PopED and option 3,4,5 stops
#' after compilation.
#'
#' When using characters, the option `compileOnly` determines if the
#' model is only compiled (and PopED is not run).
#'
#' @param compileOnly logical; only compile the model, do not run
#'   PopED (in conjunction with `iCompileOption`)
#'
#' @param iUseParallelMethod Parallel method to use
#'
#' - 0/"matlab"= Matlab PCT
#'
#' - 1/"mpi" = MPI
#'
#' @param time string that represents the time in the dataset (ie xt)
#' @param timeLow string that represents the lower design time (ie minxt)
#' @param timeHi string that represents the upper design time (ie maxmt)
#' @param id The id variable
#' @inheritParams nlmixr2est::foceiControl
#' @inheritParams PopED::create.poped.database
#' @inheritParams PopED::create_design_space
#' @inheritParams PopED::create_design
#' @param ... other parameters for PopED control
#' @return popedControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
popedControl <- function(stickyRecalcN=4,
                         maxOdeRecalc=5,
                         odeRecalcFactor=10^(0.5),
                         rxControl=NULL,
                         sigdig=4,
                         important=NULL,
                         unimportant=NULL,
                         iFIMCalculationType=c("reduced", "full", "weighted",
                                               "loc", "reducedPFIM",
                                               "fullABC", "largeMat", "reducedFIMABC"),
                         iApproximationMethod=c("fo", "foce", "focei", "foi"),
                         optSample=FALSE,
                         optDiscreteDesignVar=FALSE,
                         optContDesignVar=FALSE,
                         optIdPerGroup=FALSE,
                         iFOCENumInd=1000,
                         prior_fim=matrix(0, 0, 1),
                         d_switch=c("d", "ed"),
                         ofv_calc_type=c("lnD", "d", "a", "Ds", "inverse"),
                         strEDPenaltyFile="",
                         ofv_fun=NULL,
                         iEDCalculationType=c("mc", "laplace", "bfgs-laplace"),
                         ED_samp_size=45,
                         bLHS=c("hypercube", "random"),
                         bUseRandomSearch=TRUE,
                         bUseStochasticGradient=TRUE,
                         bUseLineSearch=TRUE,
                         bUseExchangeAlgorithm=FALSE,
                         bUseBFGSMinimizer=FALSE,
                         EACriteria=c("modified", "fedorov"),
                         strRunFile="",
                         poped_version=NULL,
                         modtit="PopED babelmixr2 model",
                         output_file="PopED_output_summary",
                         output_function_file="PopED_output_",
                         strIterationFileName="PopED_current.R",
                         user_data=NULL,
                         ourzero=1e-5,
                         dSeed=NULL,
                         line_opta=NULL,
                         line_optx=NULL,
                         bShowGraphs=FALSE,
                         use_logfile=FALSE,
                         m1_switch=c("central", "complex", "analytic", "ad"),
                         m2_switch=c("central", "complex", "analytic", "ad"),
                         hle_switch=c("central", "complex", "ad"),
                         gradff_switch=c("central", "complex", "analytic", "ad"),
                         gradfg_switch=c("central", "complex", "analytic", "ad"),
                         grad_all_switch=c("central", "complex"),
                         rsit_output=5,
                         sgit_output=1,
                         hm1 = 1e-05,
                         hlf = 1e-05,
                         hlg = 1e-05,
                         hm2 = 1e-05,
                         hgd = 1e-05,
                         hle = 1e-05,
                         AbsTol = 1e-06,
                         RelTol = 1e-06,
                         iDiffSolverMethod = NULL,
                         bUseMemorySolver = FALSE,
                         rsit = 300,
                         sgit = 150,
                         intrsit = 250,
                         intsgit = 50,
                         maxrsnullit = 50,
                         convergence_eps = 1e-08,
                         rslxt = 10,
                         rsla = 10,
                         cfaxt = 0.001,
                         cfaa = 0.001,
                         bGreedyGroupOpt = FALSE,
                         EAStepSize = 0.01,
                         EANumPoints = FALSE,
                         EAConvergenceCriteria = 1e-20,
                         bEANoReplicates = FALSE,
                         BFGSProjectedGradientTol = 1e-04,
                         BFGSTolerancef = 0.001,
                         BFGSToleranceg = 0.9,
                         BFGSTolerancex = 0.1,
                         ED_diff_it = 30,
                         ED_diff_percent = 10,
                         line_search_it = 50,
                         Doptim_iter = 1,
                         iCompileOption = c("none", "full", "mcc", "mpi"),
                         compileOnly=FALSE,
                         iUseParallelMethod = c("mpi", "matlab"),
                         MCC_Dep = NULL,
                         strExecuteName = "calc_fim.exe",
                         iNumProcesses = 2,
                         iNumChunkDesignEvals = -2,
                         Mat_Out_Pre = "parallel_output",
                         strExtraRunOptions = "",
                         dPollResultTime = 0.1,
                         strFunctionInputName = "function_input",
                         bParallelRS = FALSE,
                         bParallelSG = FALSE,
                         bParallelMFEA = FALSE,
                         bParallelLS = FALSE,
                         # design options
                         groupsize=NULL, time="time", timeLow="low", timeHi="high",
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
                         our_zero = NULL,
                         ...) {
  rxode2::rxReq("PopED")
  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% c("genRxControl"))]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste
    (paste0("'", .bad, "'", sep=""), collapse=", "),
    call.=FALSE)
  }

  if (!checkmate::testIntegerish(iFIMCalculationType, len=1, lower=0, upper=7,
                                 any.missing=FALSE)) {
    iFIMCalculationType <- match.arg(iFIMCalculationType)
    iFIMCalculationType <- c("reduced"=1, "full"=0, "weighted"=2,
                             "loc"=3, "reducedPFIM" = 4,
                             "fullABC"=5, "largeMat"=6,
                             "reducedFIMABC"=7)[iFIMCalculationType]
  }
  if (!checkmate::testIntegerish(iApproximationMethod, len=1, lower=0, upper=3)) {
    iApproximationMethod <- match.arg(iApproximationMethod)
    iApproximationMethod <- c("fo"=0, "foce"=1, "focei"=2, "foi"=3)[iApproximationMethod]
  }
  if (!checkmate::testIntegerish(m1_switch, len=1, lower=0, upper=30)) {
    m1_switch <- match.arg(m1_switch)
    m1_switch <- c("central"=1, "complex"=0, "analytic"=20, "ad"=30)[m1_switch]
  }
  if (!checkmate::testIntegerish(m2_switch, len=1, lower=0, upper=30)) {
    m2_switch <- match.arg(m2_switch)
    m2_switch <- c("central"=1, "complex"=0, "analytic"=20, "ad"=30)[m2_switch]
  }
  if (!checkmate::testIntegerish(hle_switch, len=1, lower=0, upper=30)) {
    hle_switch <- match.arg(hle_switch)
    hle_switch <- c("central"=1, "complex"=0, "ad"=30)[hle_switch]
  }
  if (!checkmate::testIntegerish(gradff_switch, len=1, lower=0, upper=30)) {
    gradff_switch <- match.arg(gradff_switch)
    gradff_switch <- c("central"=1, "complex"=0, "analytic"=20, "ad"=30)[gradff_switch]
  }
  if (!checkmate::testIntegerish(gradfg_switch, len=1, lower=0, upper=30)) {
    gradfg_switch <- match.arg(gradfg_switch)
    gradfg_switch <- c("central"=1, "complex"=0, "analytic"=20, "ad"=30)[gradfg_switch]
  }
  if (!checkmate::testIntegerish(grad_all_switch, len=1, lower=0, upper=30)) {
    grad_all_switch <- match.arg(grad_all_switch)
    grad_all_switch <- c("central"=1, "complex"=0)[grad_all_switch]
  }
  if (!checkmate::testIntegerish(ofv_calc_type, len=1, lower=1, upper=7)) {
    ofv_calc_type <- match.arg(ofv_calc_type)
    ofv_calc_type <- c("lnD"=4, "d"=1, "a"=2, "Ds"=6, "inverse"=7)[ofv_calc_type]
  }

  if (!checkmate::testIntegerish(d_switch, len=1, lower=0, upper=1)) {
    d_switch <- match.arg(d_switch)
    d_switch <- c("d"=1, "ed"=0)[d_switch]
  }
  if (!checkmate::testIntegerish(iEDCalculationType, len=1, lower=0, upper=2)) {
    iEDCalculationType <- match.arg(iEDCalculationType)
    iEDCalculationType <- c("mc"=0, "laplace"=1, "bfgs-laplace"=2)[iEDCalculationType]
  }
  if (!checkmate::testIntegerish(bLHS, len=1, lower=0, upper=1)) {
    bLHS <- match.arg(bLHS)
    bLHS <- c("hypercube"=1, "random"=0)[bLHS]
  }
  if (!checkmate::testIntegerish(EACriteria, len=1, lower=1, upper=2)) {
    EACriteria <- match.arg(EACriteria)
    EACriteria <- c("modified"=1, "fedorov"=2)[EACriteria]
  }
  if (!checkmate::testIntegerish(iCompileOption, len=1, lower= -1, upper=5)) {
    iCompileOption <- match.arg(iCompileOption)
    iCompileOption <- c("none"= -1, "full"= 0, "mcc"=1, "mpi"=2)[iCompileOption]
    if (compileOnly && iCompileOption != -1) {
      iCompileOption <- iCompileOption + 3
    }
  }
  if (!checkmate::testIntegerish(iUseParallelMethod, len=1, lower=0, upper=1)) {
    iUseParallelMethod <- match.arg(iUseParallelMethod)
    iUseParallelMethod <- c("mpi"=1, "matlab"=0)[iUseParallelMethod]
  }
  checkmate::assertIntegerish(stickyRecalcN, any.missing=FALSE, lower=0, len=1)
  checkmate::assertIntegerish(maxOdeRecalc, any.missing=FALSE, len=1)
  checkmate::assertNumeric(odeRecalcFactor, len=1, lower=1, any.missing=FALSE)
  checkmate::assertIntegerish(iFOCENumInd, len=1, lower=10, any.missing=FALSE)
  checkmate::assertMatrix(prior_fim)
  checkmate::assertIntegerish(ED_samp_size, len=1, lower=1, any.missing=FALSE)
  checkmate::assertLogical(bUseRandomSearch, len=1, any.missing = FALSE)
  checkmate::assertLogical(bUseStochasticGradient, len=1, any.missing=FALSE)
  checkmate::assertLogical(bUseLineSearch, len=1, any.missing=FALSE)
  checkmate::assertLogical(bUseExchangeAlgorithm, len=1, any.missing=FALSE)
  checkmate::assertLogical(bUseBFGSMinimizer, len=1, any.missing=FALSE)
  if (is.null(poped_version)) {
    poped_version <- packageVersion("PopED")
  }
  checkmate::assertCharacter(modtit, len=1, any.missing=FALSE, min.chars=1)
  checkmate::assertCharacter(output_file, len=1, any.missing=FALSE, min.chars = 1)
  checkmate::assertCharacter(output_function_file, len=1, any.missing = FALSE, min.chars = 1)
  checkmate::assertCharacter(strIterationFileName, len=1, any.missing=FALSE, min.chars = 1)
  checkmate::assertNumeric(ourzero, lower=0, len=1, any.missing=FALSE, finite=TRUE)

  if (!checkmate::testCharacter(strEDPenaltyFile, len=1, max.chars=1, any.missing=FALSE)) {
    checkmate::assertFileExists(strEDPenaltyFile, access="r")
  }
  if (!checkmate::testCharacter(strRunFile, len=1, max.chars=1, any.missing=FALSE)) {
    checkmate::assertFileExists(strRunFile, access="r")
  }

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
  checkmate::assertFunction(ofv_fun, null.ok=TRUE)
  if (is.null(user_data)) user_data <- PopED::cell(0, 0)
  checkmate::assertIntegerish(dSeed, any.missing=FALSE, null.ok = TRUE, len=1)
  checkmate::assertLogical(line_opta, any.missing=FALSE, null.ok = TRUE, len=1)
  checkmate::assertLogical(line_optx, any.missing=FALSE, null.ok = TRUE, len=1)
  checkmate::assertLogical(bShowGraphs, any.missing = FALSE, len=1)
  checkmate::assertLogical(use_logfile, any.missing = FALSE, len=1)
  checkmate::assertIntegerish(rsit_output, any.missing=FALSE, len=1, lower=1)
  checkmate::assertIntegerish(sgit_output, any.missing=FALSE, len=1, lower=1)
  checkmate::assertNumeric(hm1, any.missing=FALSE, len=1, lower=0)
  checkmate::assertNumeric(hlf, any.missing=FALSE, len=1, lower=0)
  checkmate::assertNumeric(hlg, any.missing=FALSE, len=1, lower=0)
  checkmate::assertNumeric(hm2, any.missing=FALSE, len=1, lower=0)
  checkmate::assertNumeric(hgd, any.missing=FALSE, len=1, lower=0)
  checkmate::assertNumeric(hle, any.missing=FALSE, len=1, lower=0)
  checkmate::assertNumeric(AbsTol, any.missing=FALSE, len=1, lower=0)
  checkmate::assertNumeric(RelTol, any.missing=FALSE, len=1, lower=0)
  checkmate::assertNumeric(convergence_eps, any.missing=FALSE, len=1, lower=0, finite=TRUE)
  checkmate::assertNumeric(cfaxt, any.missing=FALSE, len=1, lower=0, finite=TRUE)
  checkmate::assertNumeric(cfaa, any.missing=FALSE, len=1, lower=0, finite=TRUE)
  checkmate::assertNumeric(EAStepSize, any.missing=FALSE, len=1, lower=0, finite=TRUE)
  checkmate::assertNumeric(EAConvergenceCriteria, any.missing=FALSE, len=1, lower=0, finite=TRUE)
  checkmate::assertNumeric(BFGSProjectedGradientTol, any.missing=FALSE, len=1, lower=0, finite=TRUE)
  checkmate::assertNumeric(BFGSTolerancef, any.missing=FALSE, len=1, lower=0, finite=TRUE)
  checkmate::assertNumeric(BFGSToleranceg, any.missing=FALSE, len=1, lower=0, finite=TRUE)
  checkmate::assertNumeric(BFGSTolerancex, any.missing=FALSE, len=1, lower=0, finite=TRUE)
  checkmate::assertNumeric(dPollResultTime, any.missing=FALSE, len=1, lower=0, finite=TRUE)

  checkmate::testCharacter(strExecuteName, len=1, min.chars=1, any.missing=FALSE)
  checkmate::testCharacter(Mat_Out_Pre, len=1, min.chars=1, any.missing=FALSE)
  checkmate::testCharacter(strExtraRunOptions, len=1, any.missing=FALSE)

  checkmate::assertIntegerish(rsit, any.missing=FALSE, len=1, lower=1)
  checkmate::assertIntegerish(sgit, any.missing=FALSE, len=1, lower=1)
  checkmate::assertIntegerish(intrsit, any.missing=FALSE, len=1, lower=1)
  checkmate::assertIntegerish(intsgit, any.missing=FALSE, len=1, lower=1)
  checkmate::assertIntegerish(maxrsnullit, any.missing=FALSE, len=1, lower=1)
  checkmate::assertIntegerish(rslxt, any.missing=FALSE, len=1, lower=1)
  checkmate::assertIntegerish(rsla, any.missing=FALSE, len=1, lower=1)
  checkmate::assertIntegerish(ED_diff_it, any.missing=FALSE, len=1, lower=1)
  checkmate::assertIntegerish(ED_diff_percent, any.missing=FALSE, len=1, lower=1)
  checkmate::assertIntegerish(line_search_it, any.missing=FALSE, len=1, lower=1)
  checkmate::assertIntegerish(Doptim_iter, any.missing=FALSE, len=1, lower=1)
  checkmate::assertIntegerish(iNumProcesses, any.missing=FALSE, len=1, lower=1)

  checkmate::assertLogical(bUseMemorySolver, any.missing=FALSE, len=1)
  checkmate::assertLogical(bGreedyGroupOpt, any.missing=FALSE, len=1)
  checkmate::assertLogical(EANumPoints, any.missing=FALSE, len=1)
  checkmate::assertLogical(bEANoReplicates, any.missing=FALSE, len=1)
  checkmate::assertLogical(bParallelRS, any.missing=FALSE, len=1)
  checkmate::assertLogical(bParallelSG, any.missing=FALSE, len=1)
  checkmate::assertLogical(bParallelMFEA, any.missing=FALSE, len=1)
  checkmate::assertLogical(bParallelLS, any.missing=FALSE, len=1)

  .ret <- list(rxControl=rxControl,
               stickyRecalcN=as.integer(stickyRecalcN),
               maxOdeRecalc=as.integer(maxOdeRecalc),
               odeRecalcFactor=odeRecalcFactor,
               sigdig=sigdig,
               genRxControl=.genRxControl,
               optSample=optSample,
               optDiscreteDesignVar=optDiscreteDesignVar,
               optContDesignVar=optContDesignVar,
               optIdPerGroup=optIdPerGroup,
               iFIMCalculationType=iFIMCalculationType,
               iApproximationMethod=iApproximationMethod,
               iFOCENumInd=iFOCENumInd,
               prior_fim=prior_fim,
               d_switch=d_switch,
               ofv_calc_type=ofv_calc_type,
               strEDPenaltyFile=strEDPenaltyFile,
               ofv_fun=ofv_fun,
               iEDCalculationType=iEDCalculationType,
               ED_samp_size=ED_samp_size,
               bLHS=bLHS,
               bUseRandomSearch=bUseRandomSearch,
               bUseLineSearch=bUseLineSearch,
               bUseExchangeAlgorithm=bUseExchangeAlgorithm,
               bUseBFGSMinimizer=bUseBFGSMinimizer,
               EACriteria=EACriteria,
               strRunFile=strRunFile,
               poped_version=poped_version,
               modtit=modtit,
               output_file=output_file,
               output_function_file=output_function_file,
               strIterationFileName=strIterationFileName,
               user_data=user_data,
               ourzero=ourzero,
               dSeed=dSeed,
               line_opta=line_opta,
               line_optx=line_optx,
               bShowGraphs=bShowGraphs,
               use_logfile=use_logfile,
               m1_switch=m1_switch,
               m2_switch=m2_switch,
               hle_switch=hle_switch,
               gradff_switch=gradff_switch,
               gradfg_switch=gradfg_switch,
               grad_all_switch=grad_all_switch,
               rsit_output=rsit_output,
               sgit_output=sgit_output,
               hm1=hm1,
               hlf=hlf,
               hlg=hlg,
               hm2=hm2,
               hgd=hgd,
               hle=hle,
               AbsTol=AbsTol,
               RelTol = RelTol,
               iDiffSolverMethod=iDiffSolverMethod,
               bUseMemorySolver=bUseMemorySolver,
               rsit=rsit,
               sgit=sgit,
               intrsit=intrsit,
               intsgit=intsgit,
               maxrsnullit=maxrsnullit,
               convergence_eps = 1e-08,
               rslxt=rslxt,
               rsla=rsla,
               cfaxt=cfaxt,
               cfaa=cfaa,
               bGreedyGroupOpt=bGreedyGroupOpt,
               EAStepSize=EAStepSize,
               EANumPoints = EANumPoints,
               EAConvergenceCriteria = EAConvergenceCriteria,
               bEANoReplicates = bEANoReplicates,
               BFGSProjectedGradientTol = BFGSProjectedGradientTol,
               BFGSTolerancef = BFGSTolerancef,
               BFGSToleranceg = BFGSToleranceg,
               BFGSTolerancex = BFGSTolerancex,
               ED_diff_it=ED_diff_it,
               ED_diff_percent = ED_diff_percent,
               line_search_it = line_search_it,
               Doptim_iter=Doptim_iter,
               iCompileOption = iCompileOption,
               iUseParallelMethod=iUseParallelMethod,
               MCC_Dep=MCC_Dep,
               strExecuteName = strExecuteName,
               iNumProcesses = iNumProcesses,
               iNumChunkDesignEvals = iNumChunkDesignEvals,
               Mat_Out_Pre=Mat_Out_Pre,
               strExtraRunOptions = strExtraRunOptions,
               dPollResultTime=dPollResultTime,
               strFunctionInputName = strFunctionInputName,
               bParallelRS=bParallelRS,
               bParallelSG=bParallelSG,
               bParallelMFEA=bParallelMFEA,
               bParallelLS=bParallelLS,
               m=m,
               model_switch=model_switch,
               maxni=maxni,
               minni=minni,
               maxtotni=maxtotni,
               mintotni=mintotni,
               maxgroupsize=maxgroupsize,
               mingroupsize=mingroupsize,
               maxtotgroupsize=maxtotgroupsize,
               mintotgroupsize=mintotgroupsize,
               xt_space=xt_space,
               maxa=maxa,
               mina=mina,
               a_space=a_space,
               x_space=x_space,
               use_grouped_xt=use_grouped_xt,
               grouped_xt=grouped_xt,
               use_grouped_a=use_grouped_a,
               grouped_a=grouped_a,
               use_grouped_x=use_grouped_x,
               grouped_x=grouped_x,
               our_zero=our_zero)
  class(.ret) <- "popedControl"
  .ret
}
