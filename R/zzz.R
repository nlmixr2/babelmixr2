# This will be saved when compiled
rxode2.api <- names(rxode2::.rxode2ptrs())

.iniRxode2Ptr <- function() {
  .ptr <- rxode2::.rxode2ptrs()
  .nptr <- names(.ptr)
  if (length(rxode2.api) > length(.nptr)) {
    stop("babelmixr2 requires a newer version of rxode2 api, cannot run nlmixr2est\ntry `install.packages(\"rxode2\")` to get a newer version of rxode2", call.=FALSE)
  } else {
    .nptr <- .nptr[seq_along(rxode2.api)]
    if (!identical(rxode2.api, .nptr)) {
      .bad <- TRUE
      stop("babelmixr2 needs a different version of rxode2 api, cannot run nlmixr2est\ntry `install.packages(\"rxode2\")` to get a newer version of rxode2, or update both packages", call.=FALSE)
    }
  }
  .Call(`_babelmixr2_iniRxodePtrs`, .ptr,
        PACKAGE = "babelmixr2")
}




.onLoad <- function(libname, pkgname) {

  rxode2::.s3register("nlmixr2est::nmObjHandleControlObject", "bayesianToolsControl")
  rxode2::.s3register("nlmixr2est::getValidNlmixrCtl", "bayesianTools")
  rxode2::.s3register("nlmixr2est::nmObjGetControl", "bayesianTools")
  rxode2::.s3register("nlmixr2est::nlmixr2Est", "bayesianTools")

  rxode2::.s3register("nlmixr2est::nlmixr2Est", "monolix")
  rxode2::.s3register("nlmixr2est::getValidNlmixrCtl", "monolix")
  rxode2::.s3register("nlmixr2est::nmObjGetFoceiControl", "monolix")
  rxode2::.s3register("nlmixr2est::nmObjHandleControlObject", "monolixControl")
  rxode2::.s3register("nlmixr2est::nlmixr2", "pkncaEst")
  rxode2::.s3register("rxode2::rxUiDeparse", "monolixControl")
  rxode2::.s3register("rxode2::rxUiDeparse", "nonmemControl")
  rxode2::.s3register("rxode2::rxUiDeparse", "pkncaControl")
  rxode2::.s3register("rxode2::rxUiDeparse", "popedControl")
  .iniRxode2Ptr()
  .poped$loadInfo <- popedGetLoadedInfo()
}
