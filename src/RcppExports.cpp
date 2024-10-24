// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// convertDataBack
List convertDataBack(IntegerVector id, NumericVector time, NumericVector amt, NumericVector ii, IntegerVector evid, IntegerVector cmt, IntegerVector cmtDvid, IntegerVector dvidDvid, int linNcmt, int linKa, int neq, int replaceEvid, bool zeroDose2);
RcppExport SEXP _babelmixr2_convertDataBack(SEXP idSEXP, SEXP timeSEXP, SEXP amtSEXP, SEXP iiSEXP, SEXP evidSEXP, SEXP cmtSEXP, SEXP cmtDvidSEXP, SEXP dvidDvidSEXP, SEXP linNcmtSEXP, SEXP linKaSEXP, SEXP neqSEXP, SEXP replaceEvidSEXP, SEXP zeroDose2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type id(idSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type time(timeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type amt(amtSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ii(iiSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type evid(evidSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cmt(cmtSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cmtDvid(cmtDvidSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dvidDvid(dvidDvidSEXP);
    Rcpp::traits::input_parameter< int >::type linNcmt(linNcmtSEXP);
    Rcpp::traits::input_parameter< int >::type linKa(linKaSEXP);
    Rcpp::traits::input_parameter< int >::type neq(neqSEXP);
    Rcpp::traits::input_parameter< int >::type replaceEvid(replaceEvidSEXP);
    Rcpp::traits::input_parameter< bool >::type zeroDose2(zeroDose2SEXP);
    rcpp_result_gen = Rcpp::wrap(convertDataBack(id, time, amt, ii, evid, cmt, cmtDvid, dvidDvid, linNcmt, linKa, neq, replaceEvid, zeroDose2));
    return rcpp_result_gen;
END_RCPP
}
// popedGetMultipleEndpointModelingTimes
Rcpp::NumericVector popedGetMultipleEndpointModelingTimes(Rcpp::NumericVector times, Rcpp::IntegerVector modelSwitch, bool sorted);
RcppExport SEXP _babelmixr2_popedGetMultipleEndpointModelingTimes(SEXP timesSEXP, SEXP modelSwitchSEXP, SEXP sortedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type times(timesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type modelSwitch(modelSwitchSEXP);
    Rcpp::traits::input_parameter< bool >::type sorted(sortedSEXP);
    rcpp_result_gen = Rcpp::wrap(popedGetMultipleEndpointModelingTimes(times, modelSwitch, sorted));
    return rcpp_result_gen;
END_RCPP
}
// popedMultipleEndpointResetTimeIndex
Rcpp::RObject popedMultipleEndpointResetTimeIndex();
RcppExport SEXP _babelmixr2_popedMultipleEndpointResetTimeIndex() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(popedMultipleEndpointResetTimeIndex());
    return rcpp_result_gen;
END_RCPP
}
// popedMultipleEndpointIndexDataFrame
Rcpp::List popedMultipleEndpointIndexDataFrame(bool print);
RcppExport SEXP _babelmixr2_popedMultipleEndpointIndexDataFrame(SEXP printSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type print(printSEXP);
    rcpp_result_gen = Rcpp::wrap(popedMultipleEndpointIndexDataFrame(print));
    return rcpp_result_gen;
END_RCPP
}
// popedMultipleEndpointParam
Rcpp::NumericVector popedMultipleEndpointParam(Rcpp::NumericVector p, Rcpp::NumericVector times, Rcpp::IntegerVector modelSwitch, int maxMT, bool optTime);
RcppExport SEXP _babelmixr2_popedMultipleEndpointParam(SEXP pSEXP, SEXP timesSEXP, SEXP modelSwitchSEXP, SEXP maxMTSEXP, SEXP optTimeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type times(timesSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type modelSwitch(modelSwitchSEXP);
    Rcpp::traits::input_parameter< int >::type maxMT(maxMTSEXP);
    Rcpp::traits::input_parameter< bool >::type optTime(optTimeSEXP);
    rcpp_result_gen = Rcpp::wrap(popedMultipleEndpointParam(p, times, modelSwitch, maxMT, optTime));
    return rcpp_result_gen;
END_RCPP
}
// popedFree
RObject popedFree();
RcppExport SEXP _babelmixr2_popedFree() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(popedFree());
    return rcpp_result_gen;
END_RCPP
}
// popedSetup
RObject popedSetup(Environment e, Environment eglobal, bool full);
RcppExport SEXP _babelmixr2_popedSetup(SEXP eSEXP, SEXP eglobalSEXP, SEXP fullSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type e(eSEXP);
    Rcpp::traits::input_parameter< Environment >::type eglobal(eglobalSEXP);
    Rcpp::traits::input_parameter< bool >::type full(fullSEXP);
    rcpp_result_gen = Rcpp::wrap(popedSetup(e, eglobal, full));
    return rcpp_result_gen;
END_RCPP
}
// popedSolveIdME
Rcpp::DataFrame popedSolveIdME(NumericVector& theta, NumericVector& umt, NumericVector& mt, IntegerVector& ms, int nend, int id, int totn, int mn);
RcppExport SEXP _babelmixr2_popedSolveIdME(SEXP thetaSEXP, SEXP umtSEXP, SEXP mtSEXP, SEXP msSEXP, SEXP nendSEXP, SEXP idSEXP, SEXP totnSEXP, SEXP mnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type umt(umtSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type mt(mtSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type ms(msSEXP);
    Rcpp::traits::input_parameter< int >::type nend(nendSEXP);
    Rcpp::traits::input_parameter< int >::type id(idSEXP);
    Rcpp::traits::input_parameter< int >::type totn(totnSEXP);
    Rcpp::traits::input_parameter< int >::type mn(mnSEXP);
    rcpp_result_gen = Rcpp::wrap(popedSolveIdME(theta, umt, mt, ms, nend, id, totn, mn));
    return rcpp_result_gen;
END_RCPP
}
// popedSolveIdME2
Rcpp::DataFrame popedSolveIdME2(NumericVector& theta, NumericVector& umt, NumericVector& mt, IntegerVector& ms, int nend, int id, int totn, int mn);
RcppExport SEXP _babelmixr2_popedSolveIdME2(SEXP thetaSEXP, SEXP umtSEXP, SEXP mtSEXP, SEXP msSEXP, SEXP nendSEXP, SEXP idSEXP, SEXP totnSEXP, SEXP mnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type umt(umtSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type mt(mtSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type ms(msSEXP);
    Rcpp::traits::input_parameter< int >::type nend(nendSEXP);
    Rcpp::traits::input_parameter< int >::type id(idSEXP);
    Rcpp::traits::input_parameter< int >::type totn(totnSEXP);
    Rcpp::traits::input_parameter< int >::type mn(mnSEXP);
    rcpp_result_gen = Rcpp::wrap(popedSolveIdME2(theta, umt, mt, ms, nend, id, totn, mn));
    return rcpp_result_gen;
END_RCPP
}
// transDv
List transDv(NumericVector& inDv, IntegerVector& inCmt, IntegerVector& cmtTrans, NumericVector& lambda, IntegerVector& yj, NumericVector& low, NumericVector& high);
RcppExport SEXP _babelmixr2_transDv(SEXP inDvSEXP, SEXP inCmtSEXP, SEXP cmtTransSEXP, SEXP lambdaSEXP, SEXP yjSEXP, SEXP lowSEXP, SEXP highSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type inDv(inDvSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type inCmt(inCmtSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type cmtTrans(cmtTransSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type yj(yjSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type low(lowSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type high(highSEXP);
    rcpp_result_gen = Rcpp::wrap(transDv(inDv, inCmt, cmtTrans, lambda, yj, low, high));
    return rcpp_result_gen;
END_RCPP
}
