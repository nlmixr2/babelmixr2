.onLoad <- function(libname, pkgname) {
  rxode2::.s3register("nlmixr2est::nlmixr2Est", "monolix")
  rxode2::.s3register("nlmixr2est::getValidNlmixrCtl", "monolix")
}
