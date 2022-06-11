#define ARMA_DONT_PRINT_ERRORS
#define STRICT_R_HEADER
#define ARMA_DONT_PRINT_ERRORS
#define ARMA_DONT_USE_OPENMP // Known to cause speed problems
// #ifdef _OPENMP
// #include <omp.h>
// #endif
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <rxode2.h>

using namespace arma;
using namespace Rcpp;

#define _safe_sqrt(a) ((a) <= 0 ? sqrt(DBL_EPSILON) : sqrt(a))

static inline void getDv(double dv, int cmt,
                         IntegerVector &cmtTrans, NumericVector &lambda,
                         IntegerVector& yj, NumericVector& low, NumericVector& high,
                         double& llikAdj,
                         double &out, int &dvid, int &cmtOut) {
  for (unsigned int i = cmtTrans.size(); i--;) {
    if (cmt == cmtTrans[i]) {
      llikAdj += _powerL(dv, lambda[i], yj[i], low[i], high[i]);
      out = _powerD(dv, lambda[i], yj[i], low[i], high[i]);
      dvid = i+1;
      cmtOut = NA_INTEGER;
      return;
    }
  }
  out = NA_REAL;
  dvid = NA_INTEGER;
  cmtOut = cmt;
}

//[[Rcpp::export]]
List transDv(NumericVector &inDv, IntegerVector &inCmt,
             IntegerVector &cmtTrans, NumericVector &lambda,
             IntegerVector& yj, NumericVector& low, NumericVector& high) {
  NumericVector out(inDv.size());
  IntegerVector dvid(inDv.size());
  IntegerVector newCmt(inDv.size());
  double llikAdj = 0.0;
  double dvOut = 0.0;
  int dvidOut = 0;
  int cmtOut = 0;
  for (unsigned int i = inDv.size(); i--;) {
    getDv(inDv[i], inCmt[i], cmtTrans, lambda, yj, low, high, llikAdj,
          dvOut, dvidOut, cmtOut);
    out[i] = dvOut;
    dvid[i] = dvidOut;
    newCmt[i] = cmtOut;
  }
  
  return List::create(_["dv"]=out,
                      _["dvid"]=dvid,
                      _["cmt"]=newCmt,
                      _["likAdj"]=llikAdj);
}


