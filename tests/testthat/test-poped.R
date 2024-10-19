if (requireNamespace("PopED", quietly=TRUE)) {

  pheno <- function() {
    ini({
      tcl <- log(0.008) # typical value of clearance
      tv <-  log(0.6)   # typical value of volume
      ## var(eta.cl)
      eta.cl + eta.v ~ c(1,
                         0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
      # interindividual variability on clearance and volume
      add.err <- 0.1    # residual variability
    })
    model({
      cl <- exp(tcl + eta.cl) # individual value of clearance
      v <- exp(tv + eta.v)    # individual value of volume
      ke <- cl / v            # elimination rate constant
      d/dt(A1) = - ke * A1    # model differential equation
      cp = A1 / v             # concentration in plasma
      cp ~ add(add.err)       # define error model
    })
  }


  # Add modeled time for err in procedure
  phenoMT <- function() {
    ini({
      tcl <- log(0.008) # typical value of clearance
      tv <-  log(0.6)   # typical value of volume
      ## var(eta.cl)
      eta.cl + eta.v ~ c(1,
                         0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
      # interindividual variability on clearance and volume
      add.err <- 0.1    # residual variability
    })
    model({
      cl <- exp(tcl + eta.cl) # individual value of clearance
      v <- exp(tv + eta.v)    # individual value of volume
      ke <- cl / v            # elimination rate constant
      d/dt(A1) = - ke * A1    # model differential equation
      cp = A1 / v             # concentration in plasma
      mtime(v) <- 4
      cp ~ add(add.err)       # define error model
    })
  }

  phenoWt <- function() {
    ini({
      tcl <- log(0.008) # typical value of clearance
      tv <-  log(0.6)   # typical value of volume
      clw <- fix(0.75)
      ## var(eta.cl)
      eta.cl + eta.v ~ fix(1,
                           0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
      # interindividual variability on clearance and volume
      add.err <- 0.1    # residual variability
    })
    model({
      cl <- exp(tcl + eta.cl + clw * log(WT/70)) # individual value of clearance
      v <- exp(tv + eta.v)    # individual value of volume
      ke <- cl / v            # elimination rate constant
      d/dt(A1) = - ke * A1    # model differential equation
      cp = A1 / v             # concentration in plasma
      cp ~ add(add.err)       # define error model
    })
  }



  test_that("test poped interface functions", {

    # Test covariance of omega

    p <- pheno()

    expect_equal(p$popedD, c(eta.cl=1, eta.v=1))

    expect_equal(p$popedNotfixedD, NULL)

    expect_equal(p$popedCovd, 0.01)

    expect_equal(p$popedNotfixedCovd, NULL)

    phenoF <- function() {
      ini({
        tcl <- log(0.008) # typical value of clearance
        tv <-  log(0.6)   # typical value of volume
        ## var(eta.cl)
        eta.cl + eta.v ~ fix(1,
                             0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
        # interindividual variability on clearance and volume
        add.err <- 0.1    # residual variability
      })
      model({
        cl <- exp(tcl + eta.cl) # individual value of clearance
        v <- exp(tv + eta.v)    # individual value of volume
        ke <- cl / v            # elimination rate constant
        d/dt(A1) = - ke * A1    # model differential equation
        cp = A1 / v             # concentration in plasma
        cp ~ add(add.err)       # define error model
      })
    }

    p <- phenoF()

    expect_equal(p$popedD, c(eta.cl=1, eta.v=1))

    expect_equal(p$popedNotfixedD, c(0L, 0L))

    expect_equal(p$popedCovd, 0.01)

    expect_equal(p$popedNotfixedCovd, 0L)

  })

  test_that("modeled time behavior for models", {

    p <- phenoMT()

    f <- suppressMessages(.popedRxModel(p, 5))

    expect_equal(attr(f, "mtime"),
                 FALSE)

    expect_equal(attr(f, "maxNumTime"), 0L)

    expect_error(eval(f), NA)

    p <- pheno()

    f <- .popedRxModel(p, 5)

    expect_equal(attr(f, "mtime"),
                 TRUE)

    expect_equal(attr(f, "maxNumTime"), 5L)

    expect_error(eval(f), NA)

  })

  test_that("importance test; should only be needed with ofv_calc_type=6", {

    p <- pheno()

    lst <- p$popedParameters

    expect_equal(.popedImportant(p, lst, important=NULL, unimportant=NULL)$parameters$ds_index,
                 c(tcl = 0, tv = 0, eta.cl = 1, eta.v = 1, add.err = 1))

    expect_equal(.popedImportant(p, lst, important="eta.cl", unimportant=NULL)$parameters$ds_index,
                 c(tcl = 0, tv = 0, eta.cl = 0, eta.v = 1, add.err = 1))

    expect_equal(.popedImportant(p, lst, important=NULL, unimportant="tcl")$parameters$ds_index,
                 c(tcl = 1, tv = 0, eta.cl = 1, eta.v = 1, add.err = 1))

    expect_equal(.popedImportant(p, lst, important=NULL, unimportant="matt")$parameters$ds_index,
                 c(tcl = 0, tv = 0, eta.cl = 1, eta.v = 1, add.err = 1))

    expect_equal(.popedImportant(p, lst, important="matt", unimportant=NULL)$parameters$ds_index,
                 c(tcl = 0, tv = 0, eta.cl = 1, eta.v = 1, add.err = 1))

  })


  test_that("data to design", {

    withr::with_seed(42, {

      set.seed(42)
      phenoWt <- function() {
        ini({
          tcl <- log(0.008) # typical value of clearance
          tv <-  log(0.6)   # typical value of volume
          clw <- fix(0.75)
          ## var(eta.cl)
          eta.cl + eta.v ~ fix(1,
                               0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
          # interindividual variability on clearance and volume
          add.err <- 0.1    # residual variability
        })
        model({
          cl <- exp(tcl + eta.cl + clw * log(WT/70)) # individual value of clearance
          v <- exp(tv + eta.v)    # individual value of volume
          ke <- cl / v            # elimination rate constant
          d/dt(A1) = - ke * A1    # model differential equation
          cp = A1 / v             # concentration in plasma
          cp ~ add(add.err)       # define error model
        })
      }

      p <- phenoWt()

      e <- et(amt=25) %>%
        et(list(c(0, 0.5),
                c(0.5, 1.5),
                c(1.5, 2.5),
                c(2.75,6),
                c(8, 14),
                c(18, 26))) %>%
        et(id=1:2) %>%
        as.data.frame() %>%
        merge(data.frame(id=1:2, WT=c(1.4, 1.5)))

      expect_error(.popedDataToDesignSpace(p, e, groupsize=20), NA)

      db <- nlmixr2(p, e, "poped", popedControl(maxn = 15))

      library(PopED)

      expect_equal(
        list(ofv = 20.9455595871912,
             fim = lotri({
               tcl ~ 39.0235064308582
               tv ~ c(-0.385137678397696, 40.0037346438455)
               sig_var_add.err ~ 800120.483866746
             }),
             rse = c(tcl = 3.31559905061048, tv = 30.9526395967922,
                     sig_var_add.err = 11.1794980759702))
        evaluate_design(db), tolerance = 1e-4)

      ## v <- poped_optim(db, opt_xt=TRUE)

      ## expect_equal(v$ofv, 20.9503530468227, tolerance = 1e-4)


      ## skip_if_not_installed("vdiffr")
      ## vdiffr::expect_doppelganger("pheno_pred", plot_model_prediction(db, model_num_points = 300))

    })

  })


  # Test the poped intro example in nlmixr2 speak

  test_that("poped intro", {

    library(PopED)

    withr::with_seed(42, {

      set.seed(42)

      e <- et(amt=1, ii=24, until=250) %>%
        et(list(c(0, 10),
                c(0, 10),
                c(0, 10),
                c(240, 248),
                c(240, 248))) %>%
        dplyr::mutate(time =c(0, 1, 2, 8, 240, 245))

      # model
      f <- function() {
        ini({
          tV <- 72.8
          tKA <- 0.25
          tCL <- 3.75
          Favail <- fix(0.9)
          eta.ka ~ 0.09
          eta.cl ~ 0.25 ^ 2
          eta.v ~ 0.09
          prop.sd <- sqrt(0.04)
          add.sd <- fix(sqrt(5e-6))
        })
        model({
          ka <- tKA * exp(eta.ka)
          v <- tV * exp(eta.v)
          cl <- tCL * exp(eta.cl)
          d/dt(depot) <- -ka * depot
          d/dt(central) <- ka * depot - cl / v * central
          cp <- central / v
          f(depot) <- DOSE * Favail
          cp ~ add(add.sd) + prop(prop.sd)
        })
      }

      withr::with_seed(42, {

        set.seed(42)
        db <- nlmixr2(f, e, "poped",
                      popedControl(a=list(c(DOSE=20),
                                          c(DOSE=40)),
                                   maxa=c(DOSE=200),
                                   mina=c(DOSE=0)))

        vdiffr::expect_doppelganger("pred-example1",
                                    plot_model_prediction(db, model_num_points = 300))

        vdiffr::expect_doppelganger("pred-example2",
                                    plot_model_prediction(db,
                                                          PI=TRUE,
                                                          separate.groups=TRUE,
                                                          model_num_points = 300,
                                                          sample.times = FALSE))


        expect_equal(evaluate_design(db),
                     list(ofv = 39.3090057657525,
                          fim = lotri({
                            tV + tKA + tCL ~ c(0.0533669124790745, -8.68396393002656,
                                               2999.85286248872, -0.058634133188746, -14.4306063906335,
                                               37.1524325291608)
                            d_eta.ka + d_eta.cl + d_eta.v + sig_prop.var ~
                              c(439.413101383167,
                                2.2878439908508, 3412.00514559432, 312.24028527664, 3.20285380277308,
                                999.953201635217, 638.582017761863, 1182.32547832173,
                                575.347350448117, 33864.3205676804)}),
                          rse = c(tV = 8.21533739750444,
                                  tKA = 10.0909508715706, tCL = 4.40030409644082, d_eta.ka = 60.6550666795227,
                                  d_eta.cl = 27.562540815783, d_eta.v = 39.8447671522606, sig_prop.var = 13.8653576106817)), tolerance=1e-4)

      })

    })


  })

  test_that("shrinkage", {

    f <- function() {
      ini({
        tV <- 72.8
        tKa <- 0.25
        tCl <- 3.75
        tF <- fix(0.9)
        eta.v ~ 0.09
        eta.ka ~ 0.09
        eta.cl ~0.25^2
        prop.sd <- fix(sqrt(0.04))
        add.sd <- fix(sqrt(5e-6))
      })
      model({
        V<-tV*exp(eta.v)
        KA<-tKa*exp(eta.ka)
        CL<-tCl*exp(eta.cl)
        Favail <- tF
        N <-  floor(time/TAU)+1
        y <- (DOSE*Favail/V)*(KA/(KA - CL/V)) *
          (exp(-CL/V * (time - (N - 1) * TAU)) *
             (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) -
             exp(-KA * (time - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))
        y ~ prop(prop.sd) + add(add.sd)
      })
    }

    # minxt, maxxt
    e <- et(list(c(0, 10),
                 c(0, 10),
                 c(0, 10),
                 c(240, 248),
                 c(240, 248))) %>%
      as.data.frame()

    #xt
    e$time <-  c(1,2,8,240,245)


    babel.db <- nlmixr2(f, e, "poped",
                        popedControl(groupsize=20,
                                     bUseGrouped_xt=TRUE,
                                     a=list(c(DOSE=20,TAU=24),
                                            c(DOSE=40, TAU=24)),
                                     maxa=c(DOSE=200,TAU=24),
                                     mina=c(DOSE=0,TAU=24)))

    expect_error(shrinkage(babel.db), NA)

  })

  test_that("mixing 2 poped databases at the same time", {
    library(PopED)

    # Define the model
    f <- function() {
      ini({
        tV <- 72.8
        tKa <- 0.25
        tCl <- 3.75
        tF <- fix(0.9)
        pedCL <- 0.8
        eta.v ~ 0.09
        eta.ka ~ 0.09
        eta.cl ~0.25^2
        prop.sd <- fix(sqrt(0.04))
        add.sd <- fix(sqrt(5e-6))
      })
      model({
        V<-tV*exp(eta.v)
        KA<-tKa*exp(eta.ka)
        CL<-tCl*exp(eta.cl)  * (pedCL^isPediatric) # add covariate for pediatrics
        Favail <- tF
        N <-  floor(t/TAU)+1
        y <- (DOSE*Favail/V)*(KA/(KA - CL/V)) *
          (exp(-CL/V * (t - (N - 1) * TAU)) *
             (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) -
             exp(-KA * (t - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))
        y ~ prop(prop.sd) + add(add.sd)
      })
    }

    e <- et(c( 1,8,10,240,245))

    babel.db <- nlmixr2(f, e, "poped",
                        popedControl(m = 2,
                                     groupsize=20,
                                     bUseGrouped_xt=TRUE,
                                     a=list(c(DOSE=20,TAU=24,isPediatric = 0),
                                            c(DOSE=40, TAU=24,isPediatric = 0))))


    # Note, to be able to use the adults FIM to combine with the pediatrics,
    # both have to have the parameter "pedCL" defined and set notfixed_bpop to 1.

    ## Define pediatric model/design (isPediatric = 1)
    ## One arm, 4 time points only

    e.ped <- et(c( 1,2,6,240))

    babel.db.ped <- nlmixr2(f, e.ped, "poped",
                            popedControl(m = 1,
                                         groupsize=6,
                                         bUseGrouped_xt=TRUE,
                                         a=list(c(DOSE=40,TAU=24,isPediatric = 1))))


    ##  Create plot of model of adult data without variability
    vdiffr::expect_doppelganger("poped-adult-first",
                                plot_model_prediction(babel.db, model_num_points = 300))

    vdiffr::expect_doppelganger("poped-ped-next",
                                plot_model_prediction(babel.db.ped,
                                                      model_num_points = 300))

    expect_equal(babelBpopIdx(babel.db.ped, "pedCL"), 4L)

    expect_error(babelBpopIdx(babel.db, "matt"))

  })

  ##  The tests run interactively runs OK
  ## However, the capture output seems to be interfere with the tests (from PopED)
  # So... these are commented out for now.

  ## test_that("example 1", {

  ##   library(PopED)

  ##   f <- function() {
  ##     ini({
  ##       tV <- 72.8
  ##       tKa <- 0.25
  ##       tCl <- 3.75
  ##       tF <- fix(0.9)
  ##       eta.v ~ 0.09
  ##       eta.ka ~ 0.09
  ##       eta.cl ~0.25^2
  ##       prop.sd <- fix(sqrt(0.04))
  ##       add.sd <- fix(sqrt(5e-6))
  ##     })
  ##     model({
  ##       V<-tV*exp(eta.v)
  ##       KA<-tKa*exp(eta.ka)
  ##       CL<-tCl*exp(eta.cl)
  ##       Favail <- tF
  ##       N <-  floor(time/TAU)+1
  ##       y <- (DOSE*Favail/V)*(KA/(KA - CL/V)) *
  ##         (exp(-CL/V * (time - (N - 1) * TAU)) *
  ##            (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) -
  ##            exp(-KA * (time - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))
  ##       y ~ prop(prop.sd) + add(add.sd)
  ##     })
  ##   }

  ##   # minxt, maxxt
  ##   e <- et(list(c(0, 10),
  ##                c(0, 10),
  ##                c(0, 10),
  ##                c(240, 248),
  ##                c(240, 248))) %>%
  ##     as.data.frame()

  ##   #xt
  ##   e$time <-  c(1,2,8,240,245)


  ##   babel.db <- nlmixr2(f, e, "poped",
  ##                       popedControl(groupsize=20,
  ##                                    bUseGrouped_xt=TRUE,
  ##                                    a=list(c(DOSE=20,TAU=24),
  ##                                           c(DOSE=40, TAU=24)),
  ##                                    maxa=c(DOSE=200,TAU=24),
  ##                                    mina=c(DOSE=0,TAU=24)))

  ##   babelDesign <- evaluate_design(babel.db)

  ##   ##-- Model: One comp first order absorption
  ##   ## -- Analytic solution for both mutiple and single dosing
  ##   ff <- function(model_switch,xt,parameters,poped.db) {
  ##     with(as.list(parameters),{
  ##       y=xt
  ##       N = floor(xt/TAU)+1
  ##       y=(DOSE*Favail/V)*(KA/(KA - CL/V)) *
  ##         (exp(-CL/V * (xt - (N - 1) * TAU)) * (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) -
  ##            exp(-KA * (xt - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))
  ##       return(list( y=y,poped.db=poped.db))
  ##     })
  ##   }

  ##   ## -- parameter definition function
  ##   ## -- names match parameters in function ff
  ##   sfg <- function(x,a,bpop,b,bocc) {
  ##     parameters=c( V=bpop[1]*exp(b[1]),
  ##                  KA=bpop[2]*exp(b[2]),
  ##                  CL=bpop[3]*exp(b[3]),
  ##                  Favail=bpop[4],
  ##                  DOSE=a[1],
  ##                  TAU=a[2])
  ##     return( parameters )
  ##   }

  ##   ## -- Residual unexplained variablity (RUV) function
  ##   ## -- Additive + Proportional
  ##   feps <- function(model_switch,xt,parameters,epsi,poped.db){
  ##     returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))
  ##     y <- returnArgs[[1]]
  ##     poped.db <- returnArgs[[2]]
  ##     y = y*(1+epsi[,1])+epsi[,2]
  ##     return(list( y= y,poped.db =poped.db ))
  ##   }

  ##   ## -- Define design and design space
  ##   poped.db <- create.poped.database(ff_fun="ff",
  ##                                     fg_fun="sfg",
  ##                                     fError_fun="feps",
  ##                                     bpop=c(V=72.8,KA=0.25,CL=3.75,Favail=0.9),
  ##                                     notfixed_bpop=c(1,1,1,0),
  ##                                     d=c(V=0.09,KA=0.09,CL=0.25^2),
  ##                                     sigma=c(0.04,5e-6),
  ##                                     notfixed_sigma=c(0,0),
  ##                                     m=2,
  ##                                     groupsize=20,
  ##                                     xt=c( 1,2,8,240,245),
  ##                                     minxt=c(0,0,0,240,240),
  ##                                     maxxt=c(10,10,10,248,248),
  ##                                     bUseGrouped_xt=1,
  ##                                     a=list(c(DOSE=20,TAU=24),c(DOSE=40, TAU=24)),
  ##                                     maxa=c(DOSE=200,TAU=24),
  ##                                     mina=c(DOSE=0,TAU=24))

  ##   popedDesign <- evaluate_design(poped.db)

  ##   expect_equal(babelDesign$ofv, popedDesign$ofv)

  ## })

  ## test_that("example 2", {

  ##   library(PopED)

  ##   f <- function() {
  ##     ini({
  ##       tCl <- 0.15
  ##       tV <- 8
  ##       tKA <- 1.0
  ##       tFavail <- fix(1)
  ##       eta.cl ~ 0.07
  ##       eta.v ~ 0.02
  ##       eta.ka ~ 0.6
  ##       prop.sd <- sqrt(0.01)
  ##     })
  ##     model({
  ##       CL <- tCl*exp(eta.cl)
  ##       V <- tV*exp(eta.v)
  ##       KA <- tKA*exp(eta.ka)
  ##       Favail <- tFavail
  ##       y <- (DOSE*Favail*KA/(V*(KA-CL/V)))*(exp(-CL/V*time)-exp(-KA*time))
  ##       y ~ prop(prop.sd)
  ##     })
  ##   }

  ##   e <-  et(c(0.5, 1,2,6,24,36,72,120)) %>%
  ##     as.data.frame()

  ##   ## -- Define initial design  and design space
  ##   babel.db <- nlmixr2(f, e, "poped",
  ##                       control=popedControl(
  ##                         groupsize=32,
  ##                         minxt=0,
  ##                         maxxt=120,
  ##                         a=70))

  ##   babelPred <- model_prediction(babel.db)

  ##   babelED <- withr::with_seed(42, {
  ##     evaluate_design(babel.db,d_switch=FALSE,ED_samp_size=20)
  ##   })


  ##   # Example changes the model to add+prop error model

  ##   f <- function() {
  ##     ini({
  ##       tCl <- 0.15
  ##       tV <- 8
  ##       tKA <- 1.0
  ##       tFavail <- fix(1)
  ##       eta.cl ~ 0.07
  ##       eta.v ~ 0.02
  ##       eta.ka ~ 0.6
  ##       prop.sd <- sqrt(0.01)
  ##       add.sd <- sqrt(0.25)
  ##     })
  ##     model({
  ##       CL <- tCl*exp(eta.cl)
  ##       V <- tV*exp(eta.v)
  ##       KA <- tKA*exp(eta.ka)
  ##       Favail <- tFavail
  ##       y <- (DOSE*Favail*KA/(V*(KA-CL/V)))*(exp(-CL/V*time)-exp(-KA*time))
  ##       y ~ prop(prop.sd) + add(add.sd)
  ##     })
  ##   }

  ##   e <-  et(c(0.5, 1,2,6,24,36,72,120)) %>%
  ##     as.data.frame()

  ##   babel.db <- nlmixr2(f, e, "poped",
  ##                       popedControl(groupsize=32,
  ##                                    minxt=0,
  ##                                    maxxt=120,
  ##                                    a=70,
  ##                                    mina=0,
  ##                                    maxa=100,
  ##                                    # selecting important/unimportant
  ##                                    # parameters assumes Ds optimal design.
  ##                                    important=c("tCl", "tV", "tKa")))


  ##   babelDs <- withr::with_seed(42, evaluate_design(babel.db))

  ##   ff <- function(model_switch,xt,parameters,poped.db){
  ##     ##-- Model: One comp first order absorption
  ##     with(as.list(parameters),{
  ##       y=xt
  ##       y=(DOSE*Favail*KA/(V*(KA-CL/V)))*(exp(-CL/V*xt)-exp(-KA*xt))
  ##       return(list(y=y,poped.db=poped.db))
  ##     })
  ##   }

  ##   sfg <- function(x,a,bpop,b,bocc){
  ##     ## -- parameter definition function
  ##     parameters=c(CL=bpop[1]*exp(b[1]),
  ##                  V=bpop[2]*exp(b[2]),
  ##                  KA=bpop[3]*exp(b[3]),
  ##                  Favail=bpop[4],
  ##                  DOSE=a[1])
  ##     return(parameters)
  ##   }

  ##   feps <- function(model_switch,xt,parameters,epsi,poped.db){
  ##     ## -- Residual Error function
  ##     ## -- Proportional
  ##     returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))
  ##     y <- returnArgs[[1]]
  ##     poped.db <- returnArgs[[2]]
  ##     y = y*(1+epsi[,1])
  ##     return(list(y=y,poped.db=poped.db))
  ##   }

  ##   ## -- Define initial design  and design space
  ##   poped.db <- create.poped.database(ff_file="ff",
  ##                                     fg_file="sfg",
  ##                                     fError_file="feps",
  ##                                     bpop=c(CL=0.15, V=8, KA=1.0, Favail=1),
  ##                                     notfixed_bpop=c(1,1,1,0),
  ##                                     d=c(CL=0.07, V=0.02, KA=0.6),
  ##                                     sigma=0.01,
  ##                                     groupsize=32,
  ##                                     xt=c( 0.5,1,2,6,24,36,72,120),
  ##                                     minxt=0,
  ##                                     maxxt=120,
  ##                                     a=70)

  ##   popedPred <- model_prediction(poped.db)

  ##   popedED <- withr::with_seed(42, {
  ##     evaluate_design(poped.db,d_switch=FALSE,ED_samp_size=20)
  ##   })

  ##   expect_equal(popedPred$ofv, babelPred$ofv)

  ##   expect_equal(popedED$ofv, babelED$ofv)

  ##   poped.db <- create.poped.database(ff_fun=ff,
  ##                                     fg_fun=sfg,
  ##                                     fError_fun=feps.add.prop,
  ##                                     bpop=c(CL=0.15, V=8, KA=1.0, Favail=1),
  ##                                     notfixed_bpop=c(1,1,1,0),
  ##                                     d=c(CL=0.07, V=0.02, KA=0.6),
  ##                                     sigma=c(0.01,0.25),
  ##                                     groupsize=32,
  ##                                     xt=c(0.5,1,2,6,24,36,72,120),
  ##                                     minxt=0,
  ##                                     maxxt=120,
  ##                                     a=70,
  ##                                     mina=0,
  ##                                     maxa=100,
  ##                                     ds_index=c(0,0,0,1,1,1,1,1),
  ##                                     ofv_calc_type=6)
  ##   popedDs <- withr::with_seed(42, evaluate_design(poped.db))

  ##   expect_equal(popedDs$ofv, babelDs$ofv)

  ## })

  ## test_that("example 3", {

  ## })

  test_that("PopED lhs", {

    pheno <- function() {
      ini({
        tcl <- log(0.008) # typical value of clearance
        tv <-  log(0.6)   # typical value of volume
        ## var(eta.cl)
        eta.cl + eta.v ~ c(1,
                           0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
        # interindividual variability on clearance and volume
        add.err <- 0.1    # residual variability
      })
      model({
        cl <- exp(tcl + eta.cl) # individual value of clearance
        if (fed == 1)
          cl <- cl*1.1
        v <- exp(tv + eta.v)    # individual value of volume
        ke <- cl / v            # elimination rate constant
        d/dt(A1) = - ke * A1    # model differential equation
        cp = A1 / v             # concentration in plasma
        cp ~ add(add.err)       # define error model
      })
    }

    p <- pheno()

    expect_equal(eval(.popedRxModel(p))$lhs, c("rx_pred_", "rx_r_"))

  })
}
