
##
## observationTypes (list): A list giving the type of each observation present in the data file. If there is only one y-type, the corresponding observation name can be omitted.
## The possible observation types are "continuous", "discrete", and "event".
##
## nbSSDoses [optional](int): Number of doses (if there is a SS column).


##' Monolix Controller for nlmixr
##'
##' @param nbSSDoses Number of steady state doses (default 7)
##' @return
##' @author Matthew Fidler
##' @export
monolixControl <- function(nbSSDoses=7) {
  checkmate::assertIntegerish(nbSSDoses, lower=1, max.len=1)
  list(nbSSDoses=as.integer(nbSSDoses))
}

##' This constructs the data->monolix header mapping & regressors
##'
##' @param data Input dataset
##' @param regressors Character vector of model based regressors
##' @param covariates Character vector of covariates.  If the covariates are
##'   factors, assume they are categorical.  Otherwise assume continuous.
##' @return list with (header=monolix header specification, regressor=monolix regressor specification)
##' @author Matthew Fidler
monolixMapData <- function(data, uif) {
  ## Monolix header types "id", "time", "observation", "amount",
  ## "contcov", "catcov", "occ", "evid", "mdv", "obsid", "cens",
  ## "limit", "regressor","admid", "rate", "tinf", "ss", "ii", "addl",
  ## "date", "ignore".
  ##' Generate the monolix input line
  regressors <- RxODE::rxModelVars(uif$rxode)$params
  covariates <- uif$saem.all.covs
  .env <- environment()
  .env$.regressor <- c(paste0("input={", paste(regressors, collapse=", "), "}"))
  .headers <- sapply(names(data), function(x) {
    if (tolower(x) == "id") return("id")
    if (tolower(x) == "dv") return("observation")
    if (tolower(x) == "amt") return("amount")
    if (tolower(x) == "evid") return("evid")
    if (tolower(x) == "mdv") return("mdv")
    if (tolower(x) == "cmt") return("obsid")
    if (tolower(x) == "cens") return("cens")
    if (tolower(x) == "limit") return("limit")
    if (tolower(x) == "rate") return("rate")
    if (tolower(x) == "dur") return("tinf")
    if (tolower(x) == "ss") return("ss")
    if (tolower(x) == "ii") return("ii")
    if (tolower(x) == "addl") return("addl")
    if (any(x == regressors)) {
      .env$.regressor <- c(.env$.regressor, paste0(x, " = {use=regressor}"))
      return("regressor")
    }
    if (any(x == covariates)) {
      if (inherits(x, "factor") || inherits(x, "character")) {
        return("catcov")
      } else {
        return("contcov")
      }
    }
    return("ignore")
  })
  return(list(headerType=.headers, regressors=paste(.env$.regressor, collapse="\n")))
}

##' Monolix estimation routine
##'
##' @param env nlmixr environment
##'
##' @param ... Other parameters
##'
##' @export
nlmixrEst.monolix <- function(env, ...){
  with(env, {
    ## obj$env$.curTv needs to be reset
    if (length(uif$noMuEtas) > 0) {
      stop(sprintf("Cannot run Monolix since some of the parameters are not mu-referenced (%s)", paste(uif$noMuEtas, collapse = ", ")))
    }
  })
}


## lixoftConnectors::newProject(data = list(dataFile = "/path/to/data/file.txt",
##                                          headerTypes = c("IGNORE","OBSERVATION"),
##                                          observationTypes = "continuous"),
##                              modelFile = "/path/to/model/file.txt")

## Example with warfarin_data.txt from demos and oral1_1cpt_kaVCl.txt from libraries in the current directory
## data = list(dataFile= "./warfarin_data.txt",
##             headerTypes =c("id", "time", "amount", "observation", "obsid", "contcov", "catcov", "ignore"),
##             observationTypes = list(y1 = "continuous", y2 = "continuous" ),
##             mapping = list("1" = "y1", "2" = "y2"))
## modelFile <- './oral1_1cpt_kaVCl.txt'
## newProject(modelFile = modelFile, data = data)
