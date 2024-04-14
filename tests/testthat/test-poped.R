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

      expect_equal(list(ofv = 20.9439384305323,
                        fim = structure(c(38.9594470337844,
                                          -0.384587494707983, 0, -0.384587494710757, 40.0037326388712,
                                          0, 0, 0, 800137.835488059),
                                        dim = c(3L, 3L),
                                        dimnames = list(
                                          c("tcl", "tv", "sig_add.err"),
                                          c("tcl", "tv", "sig_add.err"))),
                        rse = c(tcl = 3.31832359015444, tv = 30.9526385849823, sig_add.err = 11.1793768571875)),
                   evaluate_design(db), tolerance = 1e-4)

      v <- poped_optim(db, opt_xt=TRUE)

      expect_equal(v$ofv, 20.9503530468227, tolerance = 1e-4)


      skip_if_not_installed("vdiffr")
      vdiffr::expect_doppelganger("pheno_pred", plot_model_prediction(db, model_num_points = 300))

    })

  })


  # Test the poped intro example in nlmixr2 speak

  test_that("poped intro", {

    library(PopED)

    withr::with_seed(42, {

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


        dat <- model_prediction(db,DV=TRUE)

        expect_equal(head(dat, n=4),
                     data.frame(ID = factor(c(1L, 1L, 1L, 1L), levels=paste(1:40)),
                                Time = c(1, 2, 8, 240),
                                DV = c(0.0636114439502449, 0.128497965443666, 0.146365309173223, 0.165854838702936),
                                IPRED = c(0.0682068386181826, 0.112300266786103, 0.167870981706669, 0.153239620769789),
                                PRED = c(0.0532502332862765, 0.0920480197661157, 0.164096088998621, 0.126713764327394),
                                Group = factor(c(1L, 1L, 1L, 1L), levels = c("1", "2")),
                                Model = factor(c(1L, 1L, 1L, 1L), levels = "1"),
                                DOSE = c(20, 20, 20, 20)),
                     tolerance=1e-4)

        expect_equal(evaluate_design(db),
                     list(ofv = 39.3090057657525,
                          fim = lotri({
                            tV + tKA + tCL ~ c(0.0533669124790745, -8.68396393002656,
                                               2999.85286248872, -0.058634133188746, -14.4306063906335,
                                               37.1524325291608)
                            d_eta.ka + d_eta.cl + d_eta.v + sig_prop.sd ~
                              c(439.413101383167,
                                2.2878439908508, 3412.00514559432, 312.24028527664, 3.20285380277308,
                                999.953201635217, 638.582017761863, 1182.32547832173,
                                575.347350448117, 33864.3205676804)}),
                          rse = c(tV = 8.21533739750444,
                                  tKA = 10.0909508715706, tCL = 4.40030409644082, d_eta.ka = 60.6550666795227,
                                  d_eta.cl = 27.562540815783, d_eta.v = 39.8447671522606, sig_prop.sd = 13.8653576106817)), tolerance=1e-4)

      })

    })


  })

  test_that("single endpoint solve", {

    ff_analytic <- function(model_switch,xt,parameters,poped.db){
      with(as.list(parameters),{
        y=xt
        N = floor(xt/TAU)+1
        f=(DOSE/V)*(KA/(KA - CL/V)) *
          (exp(-CL/V * (xt - (N - 1) * TAU)) * (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) -
             exp(-KA * (xt - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))
        return(list( f=f,poped.db=poped.db))
      })
    }

    sfg <- function(x,a,bpop,b,bocc){
      parameters=c(
        KA=bpop[1]*exp(b[1]),
        CL=bpop[2]*exp(b[2]),
        V=bpop[3]*exp(b[3]),
        DOSE=a[1],
        TAU=a[2])
      return( parameters )
    }

    feps <- function(model_switch,xt,parameters,epsi,poped.db){
      f <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db))[[1]]
      y = f*(1+epsi[,1])+epsi[,2]
      return(list(y=y,poped.db=poped.db))
    }

    poped_db_analytic <- create.poped.database(
      ff_fun =ff_analytic,
      fg_fun =sfg,
      fError_fun=feps,
      bpop=c(KA=0.25,CL=3.75,V=72.8),
      d=c(KA=0.09,CL=0.25^2,V=0.09),
      sigma=c(prop=0.04,add=0.0025),
      m=2,
      groupsize=20,
      xt=c( 1,2,8,240,245),
      minxt=c(0,0,0,240,240),
      maxxt=c(10,10,10,248,248),
      bUseGrouped_xt=1,
      a=cbind(DOSE=c(20,40),TAU=c(24,24)),
      maxa=c(DOSE=200,TAU=24),
      mina=c(DOSE=0,TAU=24))


    e <- et(amt=1, ii=24, until=250) %>%
      et(list(c(0, 10),
              c(0, 10),
              c(0, 10),
              c(240, 248),
              c(240, 248))) %>%
      dplyr::mutate(time =c(0, 1,2,8,240,245))

    # model
    f <- function() {
      ini({
        tKA <- 0.25
        tCL <- 3.75
        tV <- 72.8
        Favail <- fix(0.9)
        eta.ka ~ 0.09
        eta.cl ~ 0.25 ^ 2
        eta.v ~ 0.09
        prop.sd <- sqrt(0.04)
        add.sd <- sqrt(5e-6)
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

    poped_db_ode_babelmixr2 <- nlmixr(f, e,
                                      popedControl(a=list(c(DOSE=20),
                                                          c(DOSE=40)),
                                                   maxa=c(DOSE=200),
                                                   mina=c(DOSE=0)))

    eb <- evaluate_design(poped_db_ode_babelmixr2)

    expect_equal(list(ofv = 58.595188415962,
                      fim=lotri({
                        tKA + tCL + tV ~ c(2999.85286248872,
                                           -14.4306063906335, 37.1524325291608,
                                           -8.68396393002655, -0.0586341331887462, 0.0533669124790745)
                        d_eta.ka + d_eta.cl + d_eta.v + sig_prop.sd + sig_add.sd ~
                          c(439.413101383167,
                            2.2878439908508, 3412.00514559432,
                            312.24028527664, 3.20285380277308, 999.953201635217,
                            638.582017761863, 1182.32547832173, 575.347350448117, 33864.3205676804,
                            82065.2676506916, 38107.6501681295, 21309.2735081853, 2327835.71273584, 404178779.571231)
                      }),
                      rse = c(tKA = 10.0909508715706, tCL = 4.40030409644082, tV = 8.21533739750444,
                              d_eta.ka = 61.371504736694, d_eta.cl = 27.592223215635,
                              d_eta.v = 40.0643564351632,
                              sig_prop.sd = 17.6873376767579, sig_add.sd = 1297.44406776262)),
                 eb)

  })
}
