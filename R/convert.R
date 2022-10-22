.lastNobs <- 0L
.lastCmtCnt <- 0L
#' Convert nlmixr2-compatible data to other formats (if possible)
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
  nlmixr2est::nmObjUiSetCompressed(FALSE)
  on.exit({nlmixr2est::nmObjUiSetCompressed(TRUE)})
  model <- rxode2::assertRxUi(model, extra=" to convert the data with 'bblDatToMonolix'")
  model <- rxode2::rxUiDecompress(model)
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
                  replaceEvid=5L, zeroDose2 = FALSE)
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
                 function(i) {
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
  nlmixr2est::nmObjUiSetCompressed(FALSE)
  on.exit({  nlmixr2est::nmObjUiSetCompressed(TRUE)})
  .xtra <- paste0(" to convert the data with '", fun, "'")
  model <- rxode2::assertRxUi(model, extra=.xtra)
  model <- rxode2::rxUiDecompress(model)
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
                  replaceEvid=replaceEvid, zeroDose2 = TRUE)
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
  nlmixr2est::nmObjUiSetCompressed(FALSE)
  on.exit({nlmixr2est::nmObjUiSetCompressed(TRUE)})
  .xtra <- paste0(" to convert the data with 'bblDatToNonmem'")
  model <- rxode2::assertRxUi(model, extra=.xtra)
  model <- rxode2::rxUiDecompress(model)
  .ret <- .bblDatToNonmem (model, data, table,
                           fun="bblDatToNonmem", replaceEvid=5L,
                           replaceOK=FALSE, software="NONMEM", env=env)
  nlmixr2est::nmObjUiSetCompressed(FALSE)
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
  env$nmNcmt <- .lastNobs
  if (length(model$predDf$cond) > 1) {
    .dv2 <- .bblTransform(.ret$DV, .ret$CMT, model)
    .ret$DV <- .dv2$dv
    .ret$CMT <- .dv2$cmt
    .ret$DVID <- .dv2$dvid
    env$nmLikAdj <- .dv2$likAdj
    env$nmNcmt <- .dv2$nCmt
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
  .bblDatToNonmem(model, data, table,
                  fun="bblDatToRxode", replaceEvid=5L,
                  replaceOK=NA, software="rxode2", env=env)
}


#' @rdname bblDatToMonolix
#' @export
bblDatToMrgsolve <- function(model, data, table=nlmixr2est::tableControl(), env=NULL) {
  .bblDatToNonmem(model, data, table,
                  fun="bblDatToMrgsolve", replaceEvid=8L,
                  replaceOK=TRUE, software="mrgsolve", env=env)
}

#' @rdname bblDatToMonolix
#' @export
bblDatToPknca <- function(model, data, table=nlmixr2est::tableControl(), env=NULL) {
  newData <-
    .bblDatToNonmem(
      model, data, table,
      fun="bblDatToPknca", replaceEvid=5L,
      replaceOK=TRUE, software="pknca", env=env
    )
  # When we have columns in both the old and new data, which should be
  # preferred?
  betterNewCols <- c("time", "amt", "rate", "dur", "evid", "ss", "ii", "addl", "dv", "mdv", "cens", "limit")
  betterOldCols <- c("id", "cmt", "dvid")

  oldData <- data
  oldData$nlmixrRowNums <- seq_len(nrow(oldData))
  oldStdNames <- na.omit(getStandardColNames(oldData))
  newStdNames <- na.omit(getStandardColNames(newData))

  # Standardized column names to drop from the old and new data
  prepDropFromOld <- betterNewCols[betterNewCols %in% names(newStdNames)]
  prepDropFromNew <- betterOldCols[betterOldCols %in% names(oldStdNames)]
  # Actual column names to drop from the old and new data
  dropFromOld <- oldStdNames[prepDropFromOld]
  dropFromNew <- newStdNames[prepDropFromNew]

  # Prepare for merging and merge
  oldDataPrep <- oldData[, setdiff(names(oldData), dropFromOld), drop=FALSE]
  newDataPrep <- newData[, setdiff(names(newData), dropFromNew), drop=FALSE]
  stopifnot(intersect(names(oldDataPrep), names(newDataPrep)) == "nlmixrRowNums")
  # Some data may be dropped by .bblDatToNonmem above, so only keep the rows
  # that are maintained for both datasets.
  mergedData <- merge(oldDataPrep, newDataPrep, by = "nlmixrRowNums", all = FALSE)
  cleanStdNames <- getStandardColNames(mergedData)

  # Extract out the observation and dosing data
  obsData <- mergedData[mergedData[[cleanStdNames[["evid"]]]] == 0, ]
  if (nrow(obsData) < 1) {
    cli::cli_abort("no observation rows (EVID = 0) detected")
  }
  doseData <- mergedData[mergedData[[cleanStdNames[["evid"]]]] %in% c(1, 4), ]
  if (nrow(doseData) < 1) {
    cli::cli_abort("no dosing rows (EVID = 1 or 4) detected")
  }
  obsCmt <- unique(obsData[[cleanStdNames[["cmt"]]]])
  doseCmt <- unique(doseData[[cleanStdNames[["cmt"]]]])
  stopifnot(length(obsCmt) == 1)
  stopifnot(length(doseCmt) == 1)

  # Drop subjects using ADDL for dosing
  idWithAddl <- unique(doseData[[cleanStdNames["id"]]][doseData[[cleanStdNames["addl"]]] > 0])
  dropDoseAddl <- doseData[[cleanStdNames["id"]]] %in% idWithAddl
  if (any(dropDoseAddl)) {
    cli::cli_alert_info(paste("ADDL dosing not supported with PKNCA estimation, dropping subjects using ADDL:", sum(dropDoseAddl), "rows"))
    doseData <- doseData[!dropDoseAddl, ]
  }

  if (!is.na(cleanStdNames["cens"])) {
    # Drop subjects using right-censored data or LIMIT > 0 with left-censored data
    rowsRightCens <- !is.na(obsData[[cleanStdNames["cens"]]]) & obsData[[cleanStdNames["cens"]]] == -1
    rowsLeftCens <- !is.na(obsData[[cleanStdNames["cens"]]]) & obsData[[cleanStdNames["cens"]]] == 1
    idWithRightCens <- unique(obsData[[cleanStdNames["id"]]][rowsRightCens])
    idWithLeftCensNonzero <- c()
    if (!is.na(cleanStdNames["limit"])) {
      rowsWithLeftCensNonzero <-
        rowsLeftCens &
        !is.na(obsData[[cleanStdNames["limit"]]]) &
        obsData[[cleanStdNames["limit"]]] > 0
      idWithLeftCensNonzero <- unique(obsData[[cleanStdNames["id"]]][rowsWithLeftCensNonzero])
    }
    dropIdCens <- c(idWithRightCens, idWithLeftCensNonzero)
    dropObsCens <- obsData[[cleanStdNames["id"]]] %in% dropIdCens
    if (any(dropObsCens)) {
      cli::cli_alert_info(paste("Right censoring and left censoring with a value above zero is not supported with PKNCA estimation, dropping subjects with those censoring types:", sum(dropObsCens), "rows"))
      obsData <- obsData[!dropObsCens, ]
      rowsLeftCens <- rowsLeftCens[!dropObsCens]
    }

    # Modify LIMIT <= 0 so that DV is 0
    if (!is.na(cleanStdNames["limit"])) {
      rowsWithLeftCensZero <-
        rowsLeftCens &
        !is.na(obsData[[cleanStdNames["limit"]]]) &
        obsData[[cleanStdNames["limit"]]] <= 0
    } else {
      rowsWithLeftCensZero <- rowsLeftCens
    }
    if (any(rowsWithLeftCensZero)) {
      cli::cli_alert_info(paste("Setting DV to zero for PKNCA estimation with left censoring:", sum(rowsWithLeftCensZero), "rows"))
      obsData[[cleanStdNames["dv"]]][rowsWithLeftCensZero] <- 0
    }
  }

  # Drop subjects in only one dataset
  dropDoseData <- !(doseData[[cleanStdNames["id"]]] %in% unique(obsData[[cleanStdNames["id"]]]))
  dropObsData <- !(obsData[[cleanStdNames["id"]]] %in% unique(doseData[[cleanStdNames["id"]]]))
  if (any(dropDoseData)) {
    cli::cli_alert_info(paste("Dropping", sum(dropDoseData), "dosing rows with no observations for the subject with PKNCA estimation"))
    doseData <- doseData[!dropDoseData, ]
  }
  if (any(dropObsData)) {
    # This is mostly handled with the original data mapping above with
    # .bblDatToNonmem, but in some cases, the subject may be dropped subsequent
    # to that.
    cli::cli_alert_info(paste("Dropping", sum(dropObsData), "observation rows with no doses for the subject with PKNCA estimation"))
    obsData <- obsData[!dropObsData, ]
  }

  if (nrow(obsData) < 1) {
    cli::cli_abort("No observation rows after filtering for analysis")
  } else if (nrow(doseData) < 1) {
    cli::cli_abort("No dose rows after filtering for analysis")
  }
  list(
    obs=obsData,
    dose=doseData
  )
}

#' Determine standardized rxode2 column names from data
#'
#' @param data A data.frame as the source for column names
#' @return A named character vector where the names are the standardized names
#'   and the values are either the name of the column from the data or \code{NA}
#'   if the column is not present in the data.
#' @examples
#' getStandardColNames(data.frame(ID=1, DV=2, Time=3, CmT=4))
#' @export
getStandardColNames <- function(data) {
  stdCols <-
    c(
      # rxode2 columns
      "id", "time", "amt", "rate", "dur", "evid", "cmt", "ss", "ii", "addl",
      # nlmixr2 columns
      "dv", "mdv", "dvid", "cens", "limit"
    )
  stdCols <- setNames(rep(NA_character_, length(stdCols)), nm = stdCols)
  lowerNames <- tolower(names(data))
  for (nm in names(stdCols)) {
    found <- which(lowerNames %in% nm)
    if (length(found) == 1) {
      stdCols[nm] <- names(data)[found]
    } else if (length(found) > 1) {
      cli::cli_abort(paste("Multiple data columns match", nm, "when converted to lower case"))
    }
  }
  stdCols
}
