#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

SEXP _babelmixr2_convertDataBack(SEXP, SEXP, SEXP, SEXP, SEXP,
                                 SEXP, SEXP, SEXP, SEXP, SEXP,
                                 SEXP, SEXP, SEXP);

SEXP _babelmixr2_transDv(SEXP, SEXP, SEXP, SEXP, SEXP,
                         SEXP, SEXP);

SEXP _babelmixr2_iniRxodePtrs(SEXP in);

SEXP _babelmixr2_popedFree(void);
SEXP _babelmixr2_popedSetup(SEXP eSEXP, SEXP fullSEXP);
SEXP _babelmixr2_popedSolveIdME(SEXP thetaSEXP, SEXP idSEXP);
SEXP _babelmixr2_popedSolveIdME2(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP _babelmixr2_popedGetMultipleEndpointModelingTimes(SEXP, SEXP, SEXP);
SEXP _babelmixr2_popedMultipleEndpointResetTimeIndex(void);
SEXP _babelmixr2_popedMultipleEndpointIndexDataFrame(SEXP);
SEXP _babelmixr2_popedMultipleEndpointParam(SEXP, SEXP, SEXP, SEXP, SEXP,
                                            SEXP );
SEXP _babelmixr2_popedPostSolveMat(SEXP, SEXP);
SEXP _babelmixr2_popedPostSolveFull(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_babelmixr2_popedPostSolveFull", (DL_FUNC) &_babelmixr2_popedPostSolveFull, 3},
  {"_babelmixr2_popedPostSolveMat", (DL_FUNC) &_babelmixr2_popedPostSolveMat, 2},
  {"_babelmixr2_popedMultipleEndpointParam",
   (DL_FUNC) &_babelmixr2_popedMultipleEndpointParam, 6},
  {"_babelmixr2_popedMultipleEndpointIndexDataFrame", (DL_FUNC) &_babelmixr2_popedMultipleEndpointIndexDataFrame, 1},
  {"_babelmixr2_popedMultipleEndpointResetTimeIndex", (DL_FUNC) &_babelmixr2_popedMultipleEndpointResetTimeIndex, 0},
  {"_babelmixr2_popedGetMultipleEndpointModelingTimes",
   (DL_FUNC) &_babelmixr2_popedGetMultipleEndpointModelingTimes, 3},
  {"_babelmixr2_popedFree", (DL_FUNC) &_babelmixr2_popedFree, 0},
  {"_babelmixr2_popedSetup", (DL_FUNC) &_babelmixr2_popedSetup, 2},
  {"_babelmixr2_popedSolveIdME", (DL_FUNC) &_babelmixr2_popedSolveIdME, 2},
  {"_babelmixr2_popedSolveIdME2", (DL_FUNC) &_babelmixr2_popedSolveIdME2, 5},
  {"_babelmixr2_iniRxodePtrs", (DL_FUNC) &_babelmixr2_iniRxodePtrs, 1},
  {"_babelmixr2_convertDataBack", (DL_FUNC) &_babelmixr2_convertDataBack, 13},
  {"_babelmixr2_transDv", (DL_FUNC) &_babelmixr2_transDv, 7},
  {NULL, NULL, 0}
};

void R_init_babelmixr2(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
  R_forceSymbols(dll,FALSE);
}

void R_unload_babelmixr2(DllInfo *info){
}
