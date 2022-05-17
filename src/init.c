#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

SEXP _babelmixr2_convertDataBack(SEXP, SEXP, SEXP, SEXP, SEXP,
                                 SEXP, SEXP, SEXP, SEXP, SEXP,
                                 SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_babelmixr2_convertDataBack", (DL_FUNC) &_babelmixr2_convertDataBack, 12},
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

