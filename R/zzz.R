.onLoad <- function(libname, pkgname) {
  rxode2::.s3register("nlmixr2est::nlmixr2Est", "monolix")
  rxode2::.s3register("nlmixr2est::getValidNlmixrCtl", "monolix")
  rxode2::.s3register("nlmixr2est::nmObjGetFoceiControl", "monolix")
  rxode2::.s3register("nlmixr2est::nmObjHandleControlObject", "monolixControl")
}
