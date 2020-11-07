# Initial version for a model translation function
# The following considerations are in place and should be extended in future version:
# - The function only handle models with differential equations
#   although closed form should be quite easy to add
# - Meta data is not yet used, this can have a lot of added value as $DATA/$INPUT
#   can be populated when used as well as $ESTIMATION (currently set to FOCE by default)


#' nlmixr translation function
#' @examples
#'
# Test model for the function:
#' f <- function(){
#'  ini({
#'  lCl <- 1.6      #log Cl (L/hr)
#'    lVc <- log(90)  #log Vc (L)
#'    lKA <- 0.1      #log Ka (1/hr)
#'    prop.err <- c(0, 0.2, 1)
#'    add.err <- 1
#'    #eta.Cl ~ 0.1 ## BSV Cl
#'    #eta.Vc ~ 0.1 ## BSV Vc
#'    eta.Cl + eta.Vc ~ c(1,0.01, 1)
#'    eta.KA ~ fix(0.1) ## BSV Ka
#'  })
#'  model({
#'    ## First parameters are defined in terms of the initial estimates
#'    ## parameter names.
#'    Cl <- exp(lCl + eta.Cl)
#'    Vc <- exp(lVc + eta.Vc)
#'    KA <- exp(lKA + eta.KA)
#'    ## After the differential equations are defined
#'    kel <- Cl / Vc;
#'    d/dt(depot)  = -KA*depot;
#'    d/dt(centr)  =  KA*depot-kel*centr;
#'    ## And the concentration is then calculated
#'    cp = centr / Vc;
#'    ## Last, nlmixr is told that the plasma concentration follows
#'    ## a proportional error (estimated by the parameter prop.err)
#'    cp ~ prop(prop.err) + add(add.err)
#'  })
#' }
#' # Run the function
#' nlmixr_trans(f)
#' @author Richard Hooijmaijers with contributions from Matt Fidler
#' @export
nlmixr_trans <- function(func){
  parsef  <- nlmixr::nlmixr(func)

  # Get odes and additional formulas for des block
  desb    <- trimws(unlist(strsplit(parsef$nmodel$rxode,"\n")))
  desb    <- desb[!grepl("cmt\\(.*\\);",desb)] # seems like one line intended for rxode
  desboth <- desb[!grepl("d/dt\\(.*\\)",desb)]
  odespos <- grep("d/dt\\(.*\\)",desb)
  odes    <- desb[grepl("d/dt\\(.*\\)",desb)]
  odesnam <- gsub(".*\\(|\\).*","",gsub("=.*|<-.*","",odes))
  odesnum <- 1:length(odes)
  odesrhs <- sub(".*=|.*<-","",odes)
  odesnm  <- paste0("DADT(",odesnum,") =",odesrhs)
  desall           <- desb
  desall[odespos]  <- odesnm
  desall           <- stringi::stri_replace_all_fixed(desall,odesnam,paste0("A",odesnum),vectorize_all = FALSE)

  # Get information for theta block
  pars <- lapply(body(parsef$theta.pars),as.character)
  th   <- lapply(pars,function(x) if(any(grepl("THETA\\[.*\\]",x))) return(x))
  thb  <- unlist(lapply(th,function(x) paste(x[2],x[1],sub("\\[","(",sub("\\]",")",x[3])))))

  # Get information for other coding (PK block)
  et   <- lapply(pars,function(x) if(any(grepl("^ETA\\[.*\\]",x))) return(x))
  etn  <- unlist(lapply(et,"[[",2))
  etnm <- sub("\\[","(",sub("\\]",")",unlist(lapply(et,"[[",3))))
  oth  <- lapply(pars,function(x) if(all(!grepl("ETA\\[.*\\]|\\{",x))) return(x))
  oth  <- Filter(Negate(is.null), oth)
  othl <- unlist(lapply(oth,"[[",2))
  othr <- unlist(lapply(oth,"[[",3))
  othr <- stringi::stri_replace_all_fixed(othr,etn,etnm,vectorize_all = FALSE)
  oths <- paste(othl,"=",othr)

  # Get initial estimates (th/om/sig separately)
  inits     <- parsef$ini
  inits$est <- round(inits$est,4);inits$lower <- round(inits$lower,4);inits$upper <- round(inits$upper,4)
  thval     <- inits[na.omit(inits$ntheta),]
  thfix     <- ifelse(thval$fix==TRUE," FIX","")
  finit     <- parsef$focei.inits
  ombl      <- sapply(finit$OMGA,function(x){
    om   <- deparse(x)
    lom  <- length(gregexpr("ETA",om)[[1]])
    vom  <- eval(parse(text=sub("~","",deparse(x[-2]))))
    paste0("$OMEGA ",ifelse(lom>1,paste0("BLOCK(",lom,")"),""),"\n",paste(vom,collapse="\n"))
  })

  reserr  <- inits[!is.na(inits$err),]
  reserr  <- reserr[order(reserr$err,decreasing=TRUE),]
  reserrv <- paste(reserr$err,collapse="+")
  reserrg <- paste(inits$name[!is.na(inits$err)],collapse="|")
  if(reserrv=="prop"){
    reserrv <- c(thb[grepl(reserrg,thb)],paste0("W=SQRT(",reserr$name,"**2*F*F)"))
  }else if(reserrv=="add"){
    reserrv <- c(thb[grepl(reserrg,thb)],paste0("W=SQRT(",reserr$name,"**2)"))
  }else if(reserrv=="add+prop"|reserrv=="prop+add"){
    reserrv <- c(thb[grepl(reserrg,thb)],paste0("W=SQRT(",reserr$name[reserr$err=="add"],
                                                "**2+",reserr$name[reserr$err=="prop"],"**2*F*F)"))
  }

  # Create character vector with model code
  nmout <- c("$PROBLEM translated model from nlmixr","$DATA [place your dataset here]",
             "$INPUT [place your input here]","$SUBROUTINE ADVAN=13 TOL=5",
             "$MODEL",paste0("COMP=(",odesnam,")"),"$PK",thb[!grepl(reserrg,thb)],oths,"$DES",
             desall,"$ERROR (OBSERVATIONS ONLY)","IPRED=F","IRES  = DV - IPRED",reserrv,
             "IWRES = IRES/W","Y = IPRED + W*EPS(1)",
             "$THETA",paste0("(",paste(thval$lower,thval$est,thval$upper,sep=","),") ",thfix),
             ombl,"$SIGMA","1 FIX",
             "$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=5 NOABORT POSTHOC",
             "$COV COMP","$TABLE [add table variable] NOPRINT ONEHEADER FILE=partab")
  transl <- c("exp\\("="EXP(","Inf"="INF")
  for(i in 1:length(transl))  nmout <- gsub(names(transl)[i], transl[i], nmout, ignore.case = TRUE)
  cat(nmout,sep="\n")
}
