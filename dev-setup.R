withr::with_libpaths("~/nlmixr2libs", {
  if (!library(nlmixr2, logical.return=TRUE)) {
    remotes::install_github("nlmixr2/nlmixr2data")
    remotes::install_github("nlmixr2/lotri")
    remotes::install_github("nlmixr2/rxode2")
    remotes::install_github("nlmixr2/nlmixr2est")
    remotes::install_github("nlmixr2/nlmixr2extra")
    remotes::install_github("nlmixr2/nlmixr2plot")
    remotes::install_github("nlmixr2/nlmixr2")
  }
  require(nlmixr2)
}, action="prefix")



devtools::load_all()
