#define ARMA_WARN_LEVEL 1
#define STRICT_R_HEADER
#define ARMA_WARN_LEVEL 1
#define ARMA_DONT_USE_OPENMP // Known to cause speed problems
// #ifdef _OPENMP
// #include <omp.h>
// #endif
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <Rmath.h>
#include <RcppArmadillo.h>
#include <rxode2ptr.h>

#include "timsort.h"
#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <set>

struct timeInfo {
  int id;
  int firstIndex;
  int count;
};


class timeIndexer {
public:

  timeIndexer() : initialized(false), sorted(false), nIds(0) {}

  timeIndexer(const std::vector<int>& ids,
              const std::vector<std::vector<double>>& times) {
    initialize(ids, times);
  }

  void initialize(const std::vector<int>& ids, const std::vector<std::vector<double>>& times) {
    reset();
    nIds = ids.size();
    size_t k = 0;
    std::unordered_map<double, std::unordered_map<int, int>> timeIdCount;
    std::unordered_map<double, std::unordered_map<int, int>> timeIdFirstIndex;

    for (size_t i = ids.size(); i--;) {
      int id = ids[i];
      for (size_t j = 0; j < times[i].size(); ++j) {
        double time = times[i][j];
        if (timeIdCount[time].find(id) == timeIdCount[time].end()) {
          timeIdCount[time][id] = 0;
          timeIdFirstIndex[time][id] = k;
        }
        timeIdCount[time][id]++;
        uniqueTimes.insert(time);
        k++;
      }
    }

    for (const auto& timeEntry : timeIdCount) {
      double time = timeEntry.first;
      for (const auto& idEntry : timeEntry.second) {
        int id = idEntry.first;
        int count = idEntry.second;
        int firstIndex = timeIdFirstIndex[time][id];
        timeToInfo[time].emplace_back(timeInfo{id, firstIndex, count});
      }
    }

    initialized = true;
    sorted = false;
  }

  void reset() {
    timeToInfo.clear();
    uniqueTimes.clear();
    sortedTimes.clear();
    initialized = false;
    sorted = false;
    nIds=0;
  }

  size_t getNid() {
    return nIds;
  }

  bool isInitialized() const {
    return initialized;
  }

  std::vector<double> getSortedUniqueTimes() {
    if (!initialized) {
      throw std::runtime_error("timeIndexer has not been initialized");
    }
    if (!sorted) {
      sortedTimes.assign(uniqueTimes.begin(), uniqueTimes.end());
      gfx::timsort(sortedTimes.begin(), sortedTimes.end());
      sorted = true;
    }
    return sortedTimes;
  }

  const std::vector<timeInfo>& getTimeInfo(double time) const {
    if (!initialized) {
      throw std::runtime_error("timeIndexer has not been initialized");
    }
    return timeToInfo.at(time);
  }

private:
  std::unordered_map<double, std::vector<timeInfo>> timeToInfo;
  std::set<double> uniqueTimes;
  std::vector<double> sortedTimes;
  bool initialized;
  bool sorted;
  size_t nIds;
};


timeIndexer globalTimeIndexer;

void convertToTimeIndexerStructure(const std::vector<double>& times,
                                   const std::vector<int>& ids,
                                   std::vector<int>& outIds,
                                   std::vector<std::vector<double>>& outTimes) {

  std::unordered_map<int, std::vector<double>> idToTimesMap;

  for (size_t i = 0; i < times.size(); ++i) {
    idToTimesMap[ids[i]].push_back(times[i]);
  }

  for (const auto& pair : idToTimesMap) {
    outIds.push_back(pair.first);
    outTimes.push_back(pair.second);
  }
}


//' @title Get Multiple Endpoint Modeling Times
//'
//' @description
//'
//' This function takes a vector of times and a corresponding vector
//'   of IDs, groups the times by their IDs, initializes an internal
//'   C++ global TimeIndexer, that is used to efficiently lookup the
//'   final output from the rxode2 solve and then returns the sorted
//'   unique times
//'
//' @param times A numeric vector of times.
//'
//' @param modelSwitch An integer vector of model switch indicator
//'   corresponding to the times
//'
//' @return A numeric vector of sorted unique times.
//'
//' @examples
//'
//' \donttest{
//'
//' times <- c(1.1, 1.2, 1.3, 2.1, 2.2, 3.1)
//' modelSwitch <- c(1, 1, 1, 2, 2, 3)
//' sortedTimes <- popedGetMultipleEndpointModelingTimes(times, modelSwitch)
//' print(sortedTimes)
//'
//' # now show the output of the data frame representing the model
//' # switch to endpoint index
//'
//' popedMultipleEndpointIndexDataFrame()
//'
//' # now show a more complex example with overlaps etc.
//'
//' times <- c(1.1, 1.2, 1.3, 0.5, 2.2, 1.1, 0.75,0.75)
//' modelSwitch <- c(1, 1, 1, 2, 2, 2, 3, 3)
//' sortedTimes <- popedGetMultipleEndpointModelingTimes(times, modelSwitch)
//' print(sortedTimes)
//'
//' popedMultipleEndpointIndexDataFrame()
//'
//' }
// [[Rcpp::export]]
Rcpp::NumericVector popedGetMultipleEndpointModelingTimes(Rcpp::NumericVector times,
                                                          Rcpp::IntegerVector modelSwitch) {
  std::vector<int> idsVec = Rcpp::as<std::vector<int>>(modelSwitch);
  std::vector<double> timesVec = Rcpp::as<std::vector<double>>(times);

  std::vector<int> outIds;
  std::vector<std::vector<double>> outTimes;

  convertToTimeIndexerStructure(timesVec, idsVec, outIds, outTimes);

  // this will reset if needed
  globalTimeIndexer.initialize(outIds, outTimes);

  return Rcpp::wrap(globalTimeIndexer.getSortedUniqueTimes());
}

//' @title Reset the Global Time Indexer for Multiple Endpoint Modeling
//'
//' @description
//'
//' This clears the memory and resets the global time indexer used for
//'   multiple endpoint modeling.
//'
//' @return NULL, called for side effects
//'
//' @examples
//'
//' \donttest{
//'
//' popedMultipleEndpointResetTimeIndex()
//'
//' }
// [[Rcpp::export]]
Rcpp::RObject popedMultipleEndpointResetTimeIndex() {
  globalTimeIndexer.reset();
  return R_NilValue;
}

//' @rdname popedGetMultipleEndpointModelingTimes
//' @export
//[[Rcpp::export]]
Rcpp::List popedMultipleEndpointIndexDataFrame() {
  if (!globalTimeIndexer.isInitialized()) {
    Rcpp::stop("Time indexer has not been initialized");
  }
  Rcpp::NumericVector times = Rcpp::wrap(globalTimeIndexer.getSortedUniqueTimes());
  size_t nId = globalTimeIndexer.getNid();
  Rcpp::List ret(1+nId*2);
  ret[0] = times;
  // initialize the rest with NA_integer_
  for (size_t i = 0; i < nId*2; ++i) {
    Rcpp::IntegerVector cur(times.size());
    std::fill_n(cur.begin(), cur.size(), NA_INTEGER);
    ret[i+1] = cur;
  }
  size_t timei = 0;
  for (const auto& time : times) {
    const auto& infos= globalTimeIndexer.getTimeInfo(time);
    for (const auto& info : infos) {
      if (info.id > (int)nId ||
          info.id <= 0) {
        Rcpp::stop("IDs need to be sequential 1, 2, 3, ..., n");
      }
      INTEGER(ret[2*(info.id-1)+1])[timei] = info.firstIndex + 1;
      INTEGER(ret[2*(info.id-1)+2])[timei] = info.count;
    }
    timei++;
  }
  Rcpp::CharacterVector names(nId*2+1);
  names[0] = "time";
  for (size_t i = 0; i < nId*2; i+=2) {
    names[i+1] = "MS:" + std::to_string(i+1);
    names[i+2] = "N:" + std::to_string(i+1);
  }
  ret.names() = names;
  ret.attr("class") = "data.frame";
  ret.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -times.size());
  return ret;
}

#define max2( a , b )  ( (a) > (b) ? (a) : (b) )
#define isSameTime(xout, xp) (fabs((xout)-(xp))  <= DBL_EPSILON*max2(fabs(xout),fabs(xp)))

extern "C" {
#define iniRxodePtrs _babelmixr2_iniRxodePtrs
  iniRxode2ptr
}

using namespace arma;
using namespace Rcpp;

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("babelmixr2", String)
/* replace pkg as appropriate */
#else
#define _(String) (String)
#endif

#define popedOde(id) ind_solve(rx, id, rxInner.dydt_liblsoda, rxInner.dydt_lsoda_dum, rxInner.jdum_lsoda, rxInner.dydt, rxInner.update_inis, rxInner.global_jt)

struct rxSolveF {
  //
  // std::string estStr;
  // std::string gradStr;
  // std::string obfStr;
  //
  t_dydt dydt = NULL;
  t_calc_jac calc_jac = NULL;
  t_calc_lhs calc_lhs = NULL;
  t_update_inis update_inis = NULL;
  t_dydt_lsoda_dum dydt_lsoda_dum = NULL;
  t_dydt_liblsoda dydt_liblsoda = NULL;
  t_jdum_lsoda jdum_lsoda = NULL;
  t_set_solve set_solve = NULL;
  t_get_solve get_solve = NULL;
  int global_jt = 2;
  int global_mf = 22;
  int global_debug = 0;
  int neq = NA_INTEGER;
};

rxSolveF rxInner;
void rxUpdateFuns(SEXP trans, rxSolveF *inner){
  const char *lib, *s_dydt, *s_calc_jac, *s_calc_lhs, *s_inis, *s_dydt_lsoda_dum, *s_dydt_jdum_lsoda,
    *s_ode_solver_solvedata, *s_ode_solver_get_solvedata, *s_dydt_liblsoda;
  lib = CHAR(STRING_ELT(trans, 0));
  s_dydt = CHAR(STRING_ELT(trans, 3));
  s_calc_jac = CHAR(STRING_ELT(trans, 4));
  s_calc_lhs = CHAR(STRING_ELT(trans, 5));
  s_inis = CHAR(STRING_ELT(trans, 8));
  s_dydt_lsoda_dum = CHAR(STRING_ELT(trans, 9));
  s_dydt_jdum_lsoda = CHAR(STRING_ELT(trans, 10));
  s_ode_solver_solvedata = CHAR(STRING_ELT(trans, 11));
  s_ode_solver_get_solvedata = CHAR(STRING_ELT(trans, 12));
  s_dydt_liblsoda = CHAR(STRING_ELT(trans, 13));
  inner->global_jt = 2;
  inner->global_mf = 22;
  inner->global_debug = 0;
  if (strcmp(CHAR(STRING_ELT(trans, 1)),"fulluser") == 0){
    inner->global_jt = 1;
    inner->global_mf = 21;
  } else {
    inner->global_jt = 2;
    inner->global_mf = 22;
  }
  inner->calc_lhs =(t_calc_lhs) R_GetCCallable(lib, s_calc_lhs);
  inner->dydt =(t_dydt) R_GetCCallable(lib, s_dydt);
  inner->calc_jac =(t_calc_jac) R_GetCCallable(lib, s_calc_jac);
  inner->update_inis =(t_update_inis) R_GetCCallable(lib, s_inis);
  inner->dydt_lsoda_dum =(t_dydt_lsoda_dum) R_GetCCallable(lib, s_dydt_lsoda_dum);
  inner->jdum_lsoda =(t_jdum_lsoda) R_GetCCallable(lib, s_dydt_jdum_lsoda);
  inner->set_solve = (t_set_solve)R_GetCCallable(lib, s_ode_solver_solvedata);
  inner->get_solve = (t_get_solve)R_GetCCallable(lib, s_ode_solver_get_solvedata);
  inner->dydt_liblsoda = (t_dydt_liblsoda)R_GetCCallable(lib, s_dydt_liblsoda);
}

void rxClearFuns(rxSolveF *inner){
  inner->calc_lhs              = NULL;
  inner->dydt                  = NULL;
  inner->calc_jac              = NULL;
  inner->update_inis           = NULL;
  inner->dydt_lsoda_dum        = NULL;
  inner->jdum_lsoda            = NULL;
  inner->set_solve             = NULL;
  inner->get_solve             = NULL;
  inner->dydt_liblsoda         = NULL;
}

rx_solve *rx;

struct popedOptions {
  int ntheta=0;
  int stickyTol=0;
  int stickyRecalcN=1;
  int stickyRecalcN2=0;
  int stickyRecalcN1=0;
  int maxOdeRecalc;
  int reducedTol;
  int reducedTol2;
  int naZero;
  double odeRecalcFactor;
  bool loaded=false;
};

popedOptions popedOp;

//[[Rcpp::export]]
RObject popedFree() {
  return R_NilValue;
}

Environment _popedE;

//[[Rcpp::export]]
RObject popedSetup(Environment e, bool full) {
  popedFree();
  _popedE=e;
  List control = e["control"];
  List rxControl = as<List>(e["rxControl"]);

  RObject model;
  NumericVector p;
  RObject data;
  if (full) {
    model = e["modelF"];
    p = as<NumericVector>(e["paramF"]);
    data = e["dataF"]; //const RObject &events =
  } else {
    model = e["modelMT"];
    p = as<NumericVector>(e["paramMT"]);
    data = e["dataMT"]; //const RObject &events =
  }
  NumericVector p2 = p;
  std::fill_n(p2.begin(), p2.size(), NA_REAL);
  e["paramCache"]=p2;
  e["lid"] = NA_INTEGER;
  List mvp = rxode2::rxModelVars_(model);
  rxUpdateFuns(as<SEXP>(mvp["trans"]), &rxInner);

  // initial value of parameters
  CharacterVector pars = mvp[RxMv_params];
  popedOp.ntheta = pars.size();
  if (p.size() != popedOp.ntheta) {
    Rprintf("pars\n");
    print(pars);
    Rprintf("p\n");
    print(p);
    Rcpp::stop("size mismatch");
  }
  popedOp.stickyRecalcN=as<int>(control["stickyRecalcN"]);
  popedOp.stickyTol=0;
  popedOp.stickyRecalcN2=0;
  popedOp.stickyRecalcN1=0;
  popedOp.reducedTol = 0;
  popedOp.reducedTol2 = 0;
  popedOp.naZero=0;
  popedOp.maxOdeRecalc = as<int>(control["maxOdeRecalc"]);
  popedOp.odeRecalcFactor = as<double>(control["odeRecalcFactor"]);
  rxode2::rxSolve_(model, rxControl,
                   R_NilValue,//const Nullable<CharacterVector> &specParams =
                   R_NilValue,//const Nullable<List> &extraArgs =
                   p,//const RObject &params =
                   data,//const RObject &events =
                   R_NilValue, // inits
                   1);//const int setupOnly = 0
  rx = getRxSolve_();
  return R_NilValue;
}

void popedSolve(int &id) {
  rx_solving_options *op = getSolvingOptions(rx);
  rx_solving_options_ind *ind =  getSolvingOptionsInd(rx, id);
  popedOde(id);
  int j=0;
  while (popedOp.stickyRecalcN2 <= popedOp.stickyRecalcN &&
         hasOpBadSolve(op) && j < popedOp.maxOdeRecalc) {
    popedOp.stickyRecalcN2++;
    popedOp.reducedTol2 = 1;
    // Not thread safe
    rxode2::atolRtolFactor_(popedOp.odeRecalcFactor);
    setIndSolve(ind, -1);
    popedOde(id);
    j++;
  }
  if (j != 0) {
    if (popedOp.stickyRecalcN2 <= popedOp.stickyRecalcN){
      // Not thread safe
      rxode2::atolRtolFactor_(pow(popedOp.odeRecalcFactor, -j));
    } else {
      popedOp.stickyTol=1;
    }
  }
}

static inline rx_solving_options_ind* updateParamRetInd(NumericVector &theta, int &id) {
  rx = getRxSolve_();
  rx_solving_options_ind *ind = getSolvingOptionsInd(rx, id);
  for (int i = popedOp.ntheta; i--;) {
    setIndParPtr(ind, i, theta[i]);
  }
  return ind;
}

static inline int getSafeId(int id) {
  if (id < 0) return 0;
  rx = getRxSolve_();
  if (id >= getRxNsub(rx)) {
    return getRxNsub(rx)-1;
  }
  return id;
}

// Solve prediction and saved based on modeling time
void popedSolveFid(double *f, double *w, double *t, NumericVector &theta, int id, int totn) {
  // arma::vec ret(retD, nobs, false, true);
  rx_solving_options_ind *ind =  updateParamRetInd(theta, id);
  rx_solving_options *op = getSolvingOptions(rx);
  iniSubjectE(id, 1, ind, op, rx, rxInner.update_inis);
  popedSolve(id);
  int kk, k=0;
  double curT;
  for (int j = 0; j < getIndNallTimes(ind); ++j) {
    setIndIdx(ind, j);
    kk = getIndIx(ind, j);
    curT = getTime(kk, ind);
    double *lhs = getIndLhs(ind);
    if (isDose(getIndEvid(ind, kk))) {
      rxInner.calc_lhs(id, curT, getOpIndSolve(op, ind, j), lhs);
      continue;
    } else if (getIndEvid(ind, kk) == 0) {
      rxInner.calc_lhs(id, curT, getOpIndSolve(op, ind, j), lhs);
      if (ISNA(lhs[0])) {
        popedOp.naZero=1;
        lhs[0] = 0.0;
      }
      // ret(k) = lhs[0];
      // k++;
    } else if (getIndEvid(ind, kk) >= 10 && getIndEvid(ind, kk) <= 99) {
      // mtimes to calculate information
      rxInner.calc_lhs(id, curT, getOpIndSolve(op, ind, j), lhs);
      f[k] = lhs[0];
      w[k] = sqrt(lhs[1]);
      t[k] = curT;
      k++;
      if (k >= totn) return; // vector has been created, break
    }
  }
}

void popedSolveFid2(double *f, double *w, double *t, NumericVector &theta, int id, int totn) {
  // arma::vec ret(retD, nobs, false, true);
  rx_solving_options_ind *ind =  updateParamRetInd(theta, id);
  rx_solving_options *op = getSolvingOptions(rx);
  iniSubjectE(id, 1, ind, op, rx, rxInner.update_inis);
  popedSolve(id);
  int kk, k=0;
  double curT;
  for (int j = 0; j < getIndNallTimes(ind); ++j) {
    setIndIdx(ind, j);
    kk = getIndIx(ind, j);
    curT = getTime(kk, ind);
    double *lhs = getIndLhs(ind);
    if (isDose(getIndEvid(ind, kk))) {
      rxInner.calc_lhs(id, curT, getOpIndSolve(op, ind, j), lhs);
      continue;
    } else if (getIndEvid(ind, kk) == 0) {
      rxInner.calc_lhs(id, curT, getOpIndSolve(op, ind, j), lhs);
      if (ISNA(lhs[0])) {
        popedOp.naZero=1;
        lhs[0] = 0.0;
      }
      // ret(k) = lhs[0];
      // k++;
      f[k] = lhs[0];
      w[k] = sqrt(lhs[1]);
      t[k] = curT;
      k++;
      if (k >= totn) return; // vector has been created, break
    } else if (getIndEvid(ind, kk) >= 10 && getIndEvid(ind, kk) <= 99) {
      // mtimes to calculate information
      rxInner.calc_lhs(id, curT, getOpIndSolve(op, ind, j), lhs);
    }
  }
}

static inline bool solveCached(NumericVector &theta, int &id) {
  int lid = as<int>(_popedE["lid"]);
  if (lid != id) return false;
  NumericVector last = as<NumericVector>(_popedE["paramCache"]);
  return as<bool>(all(last == theta));
}

//[[Rcpp::export]]
Rcpp::DataFrame popedSolveIdN2(NumericVector &theta, NumericVector &mt, int iid, int totn) {
  int id = getSafeId(iid);
  if (solveCached(theta, id)) return(as<Rcpp::DataFrame>(_popedE["s"]));
  NumericVector t(totn);
  arma::vec f(totn);
  arma::vec w(totn);
  popedSolveFid2(&f[0], &w[0], &t[0], theta, id, totn);
  DataFrame ret = DataFrame::create(_["t"]=t,
                                    _["rx_pred_"]=f, // match rxode2/nlmixr2 to simplify code of mtime models
                                    _["w"]=w); // w = sqrt(rx_r_)
  _popedE["s"] = ret;

  return ret;
}

//[[Rcpp::export]]
Rcpp::DataFrame popedSolveIdN(NumericVector &theta, NumericVector &mt, int iid, int totn) {
  int id = getSafeId(iid);
  if (solveCached(theta, id)) return(as<Rcpp::DataFrame>(_popedE["s"]));
  NumericVector t(totn);
  arma::vec f(totn);
  arma::vec w(totn);
  popedSolveFid(&f[0], &w[0], &t[0], theta, id, totn);
  // arma::uvec m = as<arma::uvec>(match(mt, t))-1;
  // f = f(m);
  // w = w(m);
  DataFrame ret = DataFrame::create(_["t"]=t,
                                    _["rx_pred_"]=f, // match rxode2/nlmixr2 to simplify code of mtime models
                                    _["w"]=w); // w = sqrt(rx_r_)
  _popedE["s"] = ret;
  return ret;
}

void popedSolveFidMat(arma::mat &matMT, NumericVector &theta, int id, int nrow, int nend) {
  // arma::vec ret(retD, nobs, false, true);
  rx_solving_options_ind *ind =  updateParamRetInd(theta, id);
  rx_solving_options *op = getSolvingOptions(rx);
  iniSubjectE(id, 1, ind, op, rx, rxInner.update_inis);
  popedSolve(id);
  int kk, k=0;
  double curT;
  bool isMT = false;
  for (int j = 0; j < getIndNallTimes(ind); ++j) {
    setIndIdx(ind, j);
    kk = getIndIx(ind, j);
    curT = getTime(kk, ind);
    isMT = getIndEvid(ind, kk) >= 10 && getIndEvid(ind, kk) <= 99;
    double *lhs = getIndLhs(ind);
    if (isDose(getIndEvid(ind, kk))) {
      rxInner.calc_lhs(id, curT, getOpIndSolve(op, ind, j), lhs);
      continue;
    } else if (isMT) {
      // mtimes to calculate information
      rxInner.calc_lhs(id, curT, getOpIndSolve(op, ind, j), lhs);
      if (ISNA(lhs[0])) {
        popedOp.naZero=1;
        lhs[0] = 0.0;
      }
      matMT(k, 0) = curT;
      for (int i = 0; i < nend; ++i) {
        matMT(k, i*2+1) = lhs[i*2];
        matMT(k, i*2+2) = lhs[i*2+1];
      }
      k++;
      if (k >= nrow) {
        return; // vector has been created, break
      }
    } else if (getIndEvid(ind, kk) == 0) {
      rxInner.calc_lhs(id, curT, getOpIndSolve(op, ind, j), lhs);
      if (ISNA(lhs[0])) {
        popedOp.naZero=1;
        lhs[0] = 0.0;
      }
    }
  }
}

//[[Rcpp::export]]
Rcpp::DataFrame popedSolveIdME(NumericVector &theta,
                               NumericVector &umt,
                               NumericVector &mt, IntegerVector &ms,
                               int nend, int id, int totn) {
  if (solveCached(theta, id)) return(as<Rcpp::DataFrame>(_popedE["s"]));
  NumericVector t(totn);
  arma::vec f(totn);
  arma::vec w(totn);
  int nrow = umt.size();
  arma::mat matMT(nrow, nend*2+1);
  List we(nend);
  for (int i = 0; i < nend; i++) {
    we[i] = LogicalVector(totn);
  }

  popedSolveFidMat(matMT, theta, id, nrow, nend);
  // this gets the information from:
  // - the model time
  // - the model_switch
  // It will match the model switch being supplied
  // and the time being supplied
  //
  // Hence for models sent to this routine, it should match what seems
  // to be PopED's prefered way of handling the information
  // that is: ordering by model_switch, then model time, for example:
  // model_switch = 1, 1, 1, 2, 2, 2
  // time         = 0, 1, 2, 0, 1, 2
  //
  // This is not how rxode2/nlmixr2 handles the information, but this
  // routine should put it in whatever order is supplied to
  // model_switch and time
  for (int i = 0; i < totn; ++i) {
    double curT = mt[i];
    int curMS = ms[i];
    // Create a logical vector for which endpoint (used in error per endpoint identification)
    for (int j = 0; j < nend; j++) {
      LogicalVector cur = we[j];
      cur[i] = (curMS-1 == j);
      we[j] = cur;
    }
    for (int j = 0; j < nrow; ++j) {
      if (curT == matMT(j, 0)) {
        f[i] = matMT(j, (curMS-1)*2+1);
        w[i] = matMT(j, (curMS-1)*2+2);
        break;
      }
      if (j == nrow-1) {
        f[i] = NA_REAL;
        w[i] = NA_REAL;
      }
    }
  }
  DataFrame ret = DataFrame::create(_["t"]=mt,
                                    _["ms"]=ms,
                                    _["rx_pred_"]=f, // match rxode2/nlmixr2 to simplify code of mtime models
                                    _["w"]=w); // w = sqrt(rx_r_)
  _popedE["s"] = ret;
  _popedE["we"] = we;
  return ret;
}


void popedSolveFidMat2(arma::mat &matMT, NumericVector &theta, int id, int nrow, int nend) {
  // arma::vec ret(retD, nobs, false, true);
  rx_solving_options_ind *ind =  updateParamRetInd(theta, id);
  rx_solving_options *op = getSolvingOptions(rx);
  iniSubjectE(id, 1, ind, op, rx, rxInner.update_inis);
  popedSolve(id);
  int kk, k=0;
  double curT;
  for (int j = 0; j < getIndNallTimes(ind); ++j) {
    setIndIdx(ind, j);
    kk = getIndIx(ind, j);
    curT = getTime(kk, ind);
    double *lhs = getIndLhs(ind);
    if (isDose(getIndEvid(ind, kk))) {
      rxInner.calc_lhs(id, curT, getOpIndSolve(op, ind, j), lhs);
      continue;
    } else if (getIndEvid(ind, kk) == 0) {
      rxInner.calc_lhs(id, curT, getOpIndSolve(op, ind, j), lhs);
      if (ISNA(lhs[0])) {
        popedOp.naZero=1;
        lhs[0] = 0.0;
      }
      matMT(k, 0) = curT;
      for (int i = 0; i < nend; ++i) {
        matMT(k, i*2+1) = lhs[i*2];
        matMT(k, i*2+2) = lhs[i*2+1];
      }
      k++;
      if (k >= nrow) {
        return; // vector has been created, break
      }
    }
  }
}

//[[Rcpp::export]]
Rcpp::DataFrame popedSolveIdME2(NumericVector &theta,
                                NumericVector &umt,
                                NumericVector &mt, IntegerVector &ms,
                                int nend, int id, int totn) {
  if (solveCached(theta, id)) return(as<Rcpp::DataFrame>(_popedE["s"]));
  NumericVector t(totn);
  arma::vec f(totn);
  arma::vec w(totn);
  int nrow = umt.size();
  arma::mat matMT(nrow, nend*2+1);
  List we(nend);
  for (int i = 0; i < nend; i++) {
    we[i] = LogicalVector(totn);
  }

  popedSolveFidMat2(matMT, theta, id, nrow, nend);
  // arma::uvec m = as<arma::uvec>(match(mt, t))-1;
  // f = f(m);
  // w = w(m);
  for (int i = 0; i < totn; ++i) {
    double curT = mt[i];
    int curMS = ms[i];
    // Create a logical vector for which endpoint (used in error per endpoint identification)
    for (int j = 0; j < nend; j++) {
      LogicalVector cur = we[j];
      cur[i] = (curMS-1 == j);
      we[j] = cur;
    }
    for (int j = 0; j < nrow; ++j) {
      if (curT == matMT(j, 0)) {
        f[i] = matMT(j, (curMS-1)*2+1);
        w[i] = matMT(j, (curMS-1)*2+2);
        break;
      }
      if (j == nrow-1) {
        f[i] = NA_REAL;
        w[i] = NA_REAL;
      }
    }
  }
  DataFrame ret = DataFrame::create(_["t"]=mt,
                                    _["ms"]=ms,
                                    _["rx_pred_"]=f, // match rxode2/nlmixr2 to simplify code of mtime models
                                    _["w"]=w); // w = sqrt(rx_r_)
  _popedE["s"] = ret;
  _popedE["we"] = we;
  return ret;
}
