.lastNobs <- 0
#' Convert nlmixr compatible data to other formats (if possible)
#'
#' @param model rxode2 model for conversion
#'
#' @param data Input dataset.
#'
#' @param table is the table control; this is mostly to figure out if
#'   there are additional columns to keep.
#'
#' @param env When `NULL` (default) nothing is done.  When an
#'   environment, the function `nlmixr2est::.foceiPreProcessData(data,
#'   env, model)` is called on the provided environment.
#'
#' @return
#'
#' With the function `bblDatToMonolix()` return a list with:
#'  - Monolix compatible dataset ($monolix)
#'  - Monolix ADM information ($adm)
#'
#' With the function `nlmixrDataToNonmem()` return a dataset that is
#' compatible with NONMEM.
#'
#' With the function `nlmixrDataToMrgsolve()` return a dataset that is
#' compatible with `mrgsolve`.  Unlike NONMEM, it supports replacement
#' events with `evid=8` (note with `rxode2` replacement `evid` is `5`).
#'
#' With the function `nlmixrDataToRxode()` this will normalize the
#' dataset to use newer `evid` definitions that are closer to NONMEM
#' instead of any classic definitions that are used at a lower level
#'
#' @author Matthew L. Fidler
#' @export
#' @examples
#'
#' pk.turnover.emax3 <- function() {
#'   ini({
#'     tktr <- log(1)
#'     tka <- log(1)
#'     tcl <- log(0.1)
#'     tv <- log(10)
#'     ##
#'     eta.ktr ~ 1
#'     eta.ka ~ 1
#'     eta.cl ~ 2
#'     eta.v ~ 1
#'     prop.err <- 0.1
#'     pkadd.err <- 0.1
#'     ##
#'     temax <- logit(0.8)
#'     tec50 <- log(0.5)
#'     tkout <- log(0.05)
#'     te0 <- log(100)
#'     ##
#'     eta.emax ~ .5
#'     eta.ec50  ~ .5
#'     eta.kout ~ .5
#'     eta.e0 ~ .5
#'     ##
#'     pdadd.err <- 10
#'   })
#'   model({
#'     ktr <- exp(tktr + eta.ktr)
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv + eta.v)
#'     emax = expit(temax+eta.emax)
#'     ec50 =  exp(tec50 + eta.ec50)
#'     kout = exp(tkout + eta.kout)
#'     e0 = exp(te0 + eta.e0)
#'     ##
#'     DCP = center/v
#'     PD=1-emax*DCP/(ec50+DCP)
#'     ##
#'     effect(0) = e0
#'     kin = e0*kout
#'     ##
#'     d/dt(depot) = -ktr * depot
#'     d/dt(gut) =  ktr * depot -ka * gut
#'     d/dt(center) =  ka * gut - cl / v * center
#'     d/dt(effect) = kin*PD -kout*effect
#'     ##
#'     cp = center / v
#'     cp ~ prop(prop.err) + add(pkadd.err)
#'     effect ~ add(pdadd.err) | pca
#'   })
#' }
#'
#' bblDatToMonolix(pk.turnover.emax3, nlmixr2data::warfarin)
#'
#' bblDatToNonmem(pk.turnover.emax3, nlmixr2data::warfarin)
#'
#' bblDatToMrgsolve(pk.turnover.emax3, nlmixr2data::warfarin)
#'
#' bblDatToRxode(pk.turnover.emax3, nlmixr2data::warfarin)
#'
#' @useDynLib babelmixr2, .registration=TRUE
bblDatToMonolix <- function(model, data, table=nlmixr2est::tableControl(), env=NULL) {
  # https://dataset.lixoft.com/faq/translating-your-dataset-from-nonmem-format-to-the-monolix-suite-format/
  model <- rxode2::assertRxUi(model, extra=" to convert the data with 'bblDatToMonolix'")
  rxode2::assertRxUiPrediction(model, extra=" to convert the data with 'bblDatToMonolix'")
  if (is.environment(env)) {
    .env <- env
  } else {
    .env <- new.env(parent=emptyenv())
  }
  .env$table <- table
  nlmixr2est::.foceiPreProcessData(data, .env, model)
  .mv <- rxode2::rxModelVars(model)
  .flag <- .mv$flags
  .conv0 <- .Call(`_babelmixr2_convertDataBack`, .env$dataSav$ID, .env$dataSav$TIME, .env$dataSav$AMT,
                  .env$dataSav$II, .env$dataSav$EVID, .env$dataSav$CMT,
                  model$predDf$cmt, model$predDf$dvid, .flag["ncmt"], .flag["ka"], length(.mv$state),
                  replaceEvid=5L)
  assignInMyNamespace(".lastNobs", .conv0$nobs)
  if (.conv0$hasTinf && .conv0$hasRate) {
    stop("monolix does not support a fixed duration (`tinf`) and rate (`rate`) at the same time",
         call.=FALSE)
  }
  if (.conv0$hasTinf) {
    warning("monolix changes infusion times for `tinf` with bioavailability differently than `nlmixr2`, make sure there is no bioavailibilty changes for this infusion in the model",
            call.=FALSE)
  }
  if (.conv0$hasPhantom) {
    stop("transit compartment phantom events are not supported in babelmixr2 to monolix conversion",
         call.=FALSE)
  }
  if (.conv0$hasReplace) {
    stop("replacement events are not supported in monolix",
         call.=FALSE)
  }
  if (.conv0$hasMult) {
    stop("multiply events are not supported in monolix",
         call.=FALSE)
  }
  if (.conv0$hasSsRate) {
    stop("steady state infusions are not supported in monolix",
         call.=FALSE)
  }
  if (.conv0$hasSs2) {
    stop("complex steady state (ss=2) are not supported in monolix",
         call.=FALSE)
  }
  .df <- .conv0$df
  .new <- .env$dataSav
  .new$EVID <-.df$EVID
  .new$AMT <- .df$AMT
  .new$YTYPE <- .df$DVID
  .new$SS <- .df$SS
  .new$ADM <- .df$ADM
  .col0 <- c("ID", "TIME", "EVID", "AMT", "II", "DV", "ADM", "YTYPE", "SS")
  if (.conv0$hasRate) {
    .new$RATE <- .df$RATE
    .col0 <- c(.col0, "RATE")
  } else if (.conv0$hasTinf) {
    .new$TINF <- .df$TINF
    .col0 <- c(.col0, "TINF")
  }

  .censData <- NULL
  if (any(names(.new) == "CENS")) {
    .censData <- "CENS"
  }

  .limitData <- NULL
  if (any(names(.new) == "LIMIT")) {
    .limitData <- "LIMIT"
  }

  .col0 <- c(.col0, model$allCovs, .censData, .limitData, "nlmixrRowNums")
  list("monolix"=.new[.df$.nlmixrKeep, .col0],
       "adm"=.conv0$adm)
}


.bblGetLambda <- function(ui) {
  .env <- new.env(parent=emptyenv())
  .env$isFixed <- TRUE
  .predDf <- ui$predDf
  .iniDf <- ui$iniDf
  .ret <- vapply(seq_along(.predDf$cond),
                 function(i){
                   .cond <- .predDf$cond[i]
                   .df <- .iniDf[which(.iniDf$condition == .cond), ]
                   if (length(.df$cond) == 0) return(1.0)
                   .w <- which(.df$err == "yeoJohnson")
                   if (length(.w) == 1) {
                     if (!.df$fix[.w]) .env$isFixed <- FALSE
                     return(.df$est[.w])
                   }
                   .w <- which(.df$err == "boxCox")
                   if (length(.w) == 1) {
                     if (!.df$fix[.w]) .env$isFixed <- FALSE
                     return(.df$est[.w])
                   }
                   return(1.0)
                 }, numeric(1), USE.NAMES=FALSE)
  list(lambda=.ret, allFixed=.env$isFixed)
}

.bblTransform <- function(dvIn, cmtIn, ui, onlyFixed=TRUE) {
  .predDf <- ui$predDf
  .l <- .bblGetLambda(ui)
  if (onlyFixed & !.l$allFixed) {
    stop("all transformations need to be fixed (not estimated) for multiple endpoint models with babelmixr2 NONMEM",
         call.=FALSE)
  }
  .Call(`_babelmixr2_transDv`, dvIn, cmtIn, .predDf$cmt, .l$lambda,
        as.integer(.predDf$transform) - 1,
        .predDf$trLow, .predDf$trHi)
}

.bblDatToNonmem <- function(model, data, table=nlmixr2est::tableControl(),
                                fun="bblDatToNonmem", replaceEvid=5L,
                                replaceOK=FALSE, software="NONMEM", env=NULL) {
  .xtra <- paste0(" to convert the data with '", fun, "'")
  model <- rxode2::assertRxUi(model, extra=.xtra)
  rxode2::assertRxUiPrediction(model, extra=.xtra)
  if (is.environment(env)) {
    .env <- env
  } else {
    .env <- new.env(parent=emptyenv())
  }
  .env$table <- table
  nlmixr2est::.foceiPreProcessData(data, .env, model)
  .mv <- rxode2::rxModelVars(model)
  .flag <- .mv$flags
  .conv0 <- .Call(`_babelmixr2_convertDataBack`, .env$dataSav$ID, .env$dataSav$TIME, .env$dataSav$AMT,
                  .env$dataSav$II, .env$dataSav$EVID, .env$dataSav$CMT,
                  model$predDf$cmt, model$predDf$dvid, .flag["ncmt"], .flag["ka"], length(.mv$state),
                  replaceEvid=5L)
  assignInMyNamespace(".lastNobs", .conv0$nobs)
  if (!is.na(replaceOK)) {
    if (.conv0$hasPhantom) {
      stop("transit compartment phantom events are not supported in babelmixr2 to ", software, " conversion",
           call.=FALSE)
    }
    if (replaceOK) {
    } else if (.conv0$hasReplace) {
      stop("replacement events are not supported in ", software,
           call.=FALSE)
    }
    if (.conv0$hasMult) {
      stop("multiply events are not supported in ", software,
           call.=FALSE)
    }
    if (.conv0$hasTinf) {
      stop(software, "does not support a duration/tinf data item",
           call.=FALSE)
    }
  }
  .df <- .conv0$df
  .new <- .env$dataSav
  .new$EVID <-.df$EVID
  .new$AMT <- .df$AMT
  .new$DVID <- .df$DVID
  .new$SS <- .df$SS
  .new$CMT <- .df$CMT
  .col0 <- c("ID", "TIME", "EVID", "AMT", "II", "DV", "CMT", "DVID", "SS")
  if (.conv0$hasRate) {
    .new$RATE <- .df$RATE
    .col0 <- c(.col0, "RATE")
  } else if (.conv0$hasTinf) {
    .new$TINF <- .df$TINF
    .col0 <- c(.col0, "TINF")
  }

  .censData <- NULL
  .w <- which(toupper(names(.new)) == "CENS")
  if (length(.w) == 1) {
    names(.new)[.w] <- "CENS"
    .censData <- "CENS"
  }

  .limitData <- NULL
  .w <- which(toupper(names(.new)) == "LIMIT")
  if (length(.w) == 1) {
    names(.new)[.w] <- "LIMIT"
    .limitData <- "LIMIT"
  }

  .col0 <- c(.col0, model$allCovs, .censData, .limitData, "nlmixrRowNums")
  .new[.df$.nlmixrKeep, .col0]
}

#' @rdname bblDatToMonolix
#' @export
bblDatToNonmem <- function(model, data, table=nlmixr2est::tableControl(), env=NULL) {
  .xtra <- paste0(" to convert the data with 'bblDatToNonmem'")
  model <- rxode2::assertRxUi(model, extra=.xtra)
  .ret <- .bblDatToNonmem (model, data, table,
                           fun="bblDatToNonmem", replaceEvid=5L,
                           replaceOK=FALSE, software="NONMEM", env=env)

  .ret <- .ret[, names(.ret) != "DVID"]
  if (any(names(.ret) == "LIMIT")) {
    # This converts LIMIT to NONMEM's definition of infinity
    # (according to manual for $THETA)
    .ret$LIMIT <- ifelse(is.finite(.ret$LIMIT),
                         .ret$LIMIT,
                         ifelse(.ret$LIMIT < 0, -1000000, 1000000))
  }
  .ui <- model
  env$nobs <- .lastNobs
  env$nmLikAdj <- 0
  if (length(model$predDf$cond) > 1) {
    .dv2 <- .bblTransform(.ret$DV, .ret$CMT, model)
    .ret$DV <- .dv2$dv
    .ret$CMT <- .dv2$cmt
    .ret$DVID <- .dv2$dvid
    env$nmLikAdj <- .dv2$likAdj
  }
  .names <- c(
    "ID", "TIME", "EVID", "AMT",
    ifelse(rxode2::rxGetControl(.ui, ".hasIi", FALSE), "II", ""),
    "DV", "CMT",
    ifelse(length(.ui$predDf$cond) > 1, "DVID", ""),
    ifelse(rxode2::rxGetControl(.ui, ".hasSs", FALSE), "SS", ""),
    ifelse(rxode2::rxGetControl(.ui, ".hasRate", FALSE), "RATE", ""),
    vapply(.ui$allCovs, .nmGetVar, character(1), ui=.ui,
           USE.NAMES=FALSE),
    ifelse(rxode2::rxGetControl(.ui, ".hasCens", FALSE), "CENS", ""),
    ifelse(rxode2::rxGetControl(.ui, ".hasLimit", FALSE), "LIMIT", ""),
    "nlmixrRowNums")
  .names <- .names[.names != ""]
  .ret <- .ret[, .names]
  .ret
}

#' @rdname bblDatToMonolix
#' @export
bblDatToRxode <- function(model, data, table=nlmixr2est::tableControl(), env=NULL) {
  .bblDatToNonmem (model, data, table,
                       fun="bblDatToRxode", replaceEvid=5L,
                       replaceOK=NA, software="rxode2", env=env)
}


#' @rdname bblDatToMonolix
#' @export
bblDatToMrgsolve <- function(model, data, table=nlmixr2est::tableControl(), env=NULL) {
  .bblDatToNonmem (model, data, table,
                       fun="bblDatToMrgsolve", replaceEvid=8L,
                       replaceOK=TRUE, software="mrgsolve", env=env)
}
