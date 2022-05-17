#define STRICT_R_HEADER
#include <Rcpp.h>

// EVID = 0; Observations
// EVID = 1; is illegal, but converted from NONMEM
// EVID = 2; Non-observation, possibly covariate
// EVID = 3; Reset ODE states to zero; Non-observation event
// EVID = 4; Reset and then dose event;  Illegal
// EVID = 9; Non-observation event to ini system at time zero; This is to set the INIs at the correct place.
// EVID = 10-99; mtime events (from ODE system)
// When EVID > 100
// EVID: ## # ## ##
//       c2 I c1 xx
// c2 = Compartment numbers over 100
//  I = Infusion Flag/ Special event flag
#define EVIDF_NORMAL 0

#define EVIDF_INF_RATE 1
#define EVIDF_INF_DUR  2

#define EVIDF_REPLACE  4
#define EVIDF_MULT     5

#define EVIDF_MODEL_DUR_ON   8
#define EVIDF_MODEL_DUR_OFF  6

#define EVIDF_MODEL_RATE_ON  9
#define EVIDF_MODEL_RATE_OFF 7
//      0 = no Infusion
//      1 = Infusion, AMT=rate (mg/hr for instance)
//      2 = Infusion, duration is fixed
//      4 = Replacement event
//      5 = Multiplication event
//      6 = Turn off modeled duration
//      7 = Turn off modeled rate compartment
//      8 = Duration is modeled, AMT=dose; Rate = AMT/(Modeled Duration) NONMEM RATE=-2
//      9 = Rate is modeled, AMT=dose; Duration = AMT/(Modeled Rate) NONMEM RATE=-1
// c1 = Compartment numbers below 99
// xx = 1, regular event
// xx = 10, steady state event SS=1
// xx = 20, steady state event + last observed info.
// xx = 30, Turn off compartment
// xx = 40, Steady state constant infusion
// xx = 50, Phantom event, used for transit compartments
// Steady state events need a II data item > 0
#define EVID0_REGULAR 1
#define EVID0_SS 10
#define EVID0_SS2 20
#define EVID0_OFF 30
#define EVID0_SSINF 40 
#define EVID0_PHANTOM 50

static inline void getWh(int evid, int *wh, int *cmt, int *wh100, int *whI, int *wh0,
                         int linNcmt, int linKa, int neq){
  *wh = evid;
  *cmt = 0;
  *wh100 = std::floor(*wh/1e5L);
  *whI   = std::floor(*wh/1e4L-*wh100*10);
  *wh    = *wh - *wh100*1e5 - (*whI-1)*1e4;
  *wh0 = std::floor((*wh%10000)/100);
  *cmt = *wh0 - 1 + *wh100*100;
  *wh0 = evid - *wh100*1e5 - *whI*1e4 - *wh0*100;
  if (linNcmt != 0) {
    if (linKa) {
      switch (*cmt) {
      case 0:
				*cmt = neq;
				break;
      case 1:
				*cmt = neq+1;
				break;
      case 2:
				*cmt -= 2;
				break;
      }
    } else {
      if (*cmt == 0) {
				*cmt = neq;
      } else {
				*cmt -= 1;
      }
    }
  }
}

static inline int getSs(int wh0, bool &hasSs, bool &hasSs2, bool &hasSsRate) {
  if (wh0 == EVID0_SS2) {
    hasSs=true;
    return  2;
  } else if (wh0 == EVID0_SS) {
    hasSs=true;
    return 1;
  } else if (wh0 == EVID0_SSINF) {
    // corresponds to ss=1, amt=0, ii=0
    // both amt and ii should be 0, and rate should contain the infusion rate
    hasSsRate=false;
    hasSs=true;
    return 1;
  }
  return 0;
}

using namespace Rcpp;

static inline int getDvid(int &cmt, IntegerVector &dvidDvid, IntegerVector &cmtDvid) {
  for (int j=cmtDvid.size(); j--;) {
    if (cmtDvid[j] == cmt) {
      return(dvidDvid[j]);
    }
  }
  return 0;
}


#define MONOLIX_EMPTY_TYPE 1
#define MONOLIX_MODEL_RATE 2
#define MONOLIX_MODEL_DUR  3
#define MONOLIX_INFUSION   4
#define MONOLIX_BOLUS      5

static inline int getAdm(int cmt, int type, std::vector<int> &admIds) {
  int id = cmt*10+type;
  for (int i = 0; i < admIds.size(); ++i) {
    if (admIds[i] == id) {
      return i+1;
    }
  }
  admIds.push_back(id);
  return admIds.size();
}

static inline DataFrame createAdm(std::vector<int> &admIds) {
  IntegerVector adm(admIds.size());
  IntegerVector cmt(admIds.size());
  IntegerVector type(admIds.size());
  int c=0, t=0;
  for (int i = 0; i < admIds.size(); ++i) {
    adm[i] = i+1;
    c =admIds[i]/10;
    t = admIds[i] - c*10;
    cmt[i] = c;
    type[i] = t;
  }
  type.attr("levels") = CharacterVector::create("empty", "modelRate", "modelDur", "infusion", "bolus");
  type.attr("class") = "factor";
  return DataFrame::create(_["adm"]=adm,
                           _["cmt"]=cmt,
                           _["type"]=type);
}

//[[Rcpp::export]]
List convertDataBack(IntegerVector id, NumericVector time, NumericVector amt, NumericVector ii,
                     IntegerVector evid, IntegerVector cmt,
                     IntegerVector cmtDvid, IntegerVector dvidDvid,
                     int linNcmt=0, int linKa=0, int neq=0,
                     int replaceEvid=5) {
  IntegerVector newEvid(evid.size());
  IntegerVector newSs(evid.size());
  IntegerVector newDvid(evid.size());
  LogicalVector keepItem(evid.size());
  NumericVector newAmt(evid.size());
  NumericVector newRate(evid.size());
  NumericVector newTinf(evid.size());
  IntegerVector newCmt(evid.size());
  IntegerVector newAdm(evid.size());
  std::vector<int> admIds;
  int wh=0, cmt0=0, wh100=0, whI=0, wh0=0;
  bool turnOffCmt=false;
  bool hasTinf = false;
  bool hasRate = false;
  bool hasPhantom = false;
  bool hasReplace = false;
  bool hasMult = false;
  bool hasSs=false;
  bool hasSs2=false;
  bool hasSsRate=false;
  double curAmt=0.0;
  for (unsigned int i = 0; i < evid.size(); ++i) {
    int curEvid = evid[i];
    // put in defaults
    newEvid[i] = evid[i];
    newSs[i] = 0;
    newDvid[i] = 0;
    newRate[i] = 0;
    newTinf[i] = 0;
    newAmt[i] = amt[i];
    newCmt[i] = cmt[i];
    newAdm[i] = 0;
    if (curEvid == 0 || curEvid==2 || curEvid == 3) {
      // 0, 2 and 3 are preserved
      newDvid[i] = curEvid == 3 ? 0 : getDvid(cmt[i], dvidDvid, cmtDvid);
      keepItem[i] = true;
    } else if ((curEvid) >= 9 && (curEvid) <= 99) {
      // mtimes and zero events are ignored
      keepItem[i] = false;
    } else {
      // these are doses
      getWh(curEvid, &wh, &cmt0, &wh100, &whI, &wh0,
            linNcmt, linKa, neq);
      if (wh0 ==EVID0_OFF) {
        // turn off a compartment; supported in nonmem
        // In monoix "Turning off compartments should instead be defined in the model file"
        // This is done by the macro empty(adm=X, target=CmtName)
        // MONOLIX_EMPTY_TYPE
        newCmt[i] = -cmt[i];
        newAdm[i] = getAdm(cmt[i], MONOLIX_EMPTY_TYPE, admIds);
        turnOffCmt=true;
        keepItem[i] = true;
      } else if (wh0 == EVID0_PHANTOM) {
        // In monolix this can be defined by:
        // depot(type=x, target=cmtName, ka, Ktr, Mtt)
        
        // But the ODE structure must be completely change so we will
        // not support this for now
        keepItem[i] = false;
        newEvid[i] = 7;
        hasPhantom = true;
      } else {
        switch (whI) {
        case EVIDF_MODEL_RATE_ON: // modeled rate.
          // In monolix this is also done by the macro; rates cannot be -1 
          newEvid[i] = 1;
          newSs[i] =getSs(wh0, hasSs, hasSs2, hasSsRate);
          newRate[i] = -1;
          newAmt[i] = amt[i];
          newAdm[i] = getAdm(cmt[i], MONOLIX_MODEL_RATE, admIds);
          keepItem[i] = true;
          break;
        case EVIDF_MODEL_DUR_ON: // modeled duration.
          // In monolix this is done by macro; rates cannot be -2
          // Rather you define
          // depot(type=ADMId, target=Cmt, ModeledDuration)
          newEvid[i] = 1;
          newSs[i] = getSs(wh0, hasSs, hasSs2, hasSsRate);
          newRate[i] = -2;
          newAmt[i] = amt[i];
          newAdm[i] = getAdm(cmt[i], MONOLIX_MODEL_DUR, admIds);
          keepItem[i] = true;
          break;
        case EVIDF_MODEL_RATE_OFF: // End modeled rate
        case EVIDF_MODEL_DUR_OFF: // end modeled duration
          keepItem[i] = false;
          break;
        case EVIDF_INF_DUR:
        case EVIDF_INF_RATE:
          // classic rxode2 uses rate instead of amount here
          // in monolix this can be denfined by iv() macro
          curAmt = amt[i];
          if (amt[i] > 0) {
            bool foundOff=false;
            double t2=0;
            int curId = id[i];
            for (int j=i; j < evid.size(); j++){
              if (id[j] != curId) {
                break;
              }
              if (evid[j] == curEvid && curAmt == -amt[j]) {
                t2 = time[j];
                foundOff = true;
                break;
              }
            }
            if (foundOff) {
              double dur = t2-time[i];
              double rate = amt[i];
              curAmt = rate*dur;
              if (whI==EVIDF_INF_DUR) {
                newTinf[i] = dur;
                hasTinf = true;
              } else {
                newRate[i] = rate;
                hasRate = true;
              }
              newEvid[i] = 1;
              newSs[i] = getSs(wh0, hasSs, hasSs2, hasSsRate);
              newAmt[i] = curAmt;
              newAdm[i] = getAdm(cmt[i], MONOLIX_INFUSION, admIds);
              keepItem[i] = true;
            } else {
              keepItem[i] = false;
            }
          } else {
            keepItem[i] = false;
          }
        case EVIDF_REPLACE:
          newEvid[i] = replaceEvid;
          keepItem[i] = false;
          hasReplace=true;
        case EVIDF_MULT:
          newEvid[i] = 6;
          newSs[i] =getSs(wh0, hasSs, hasSs2, hasSsRate);
          keepItem[i] = false;
          hasMult=true;
        case EVIDF_NORMAL:
          // defined by 
          newEvid[i] = 1;
          newSs[i] =getSs(wh0, hasSs, hasSs2, hasSsRate);
          newAdm[i] = getAdm(cmt[i], MONOLIX_BOLUS, admIds);
          keepItem[i] = true;          
        }
      }
    }
  }
  // now return the dataset
  return List::create(_["df"]=DataFrame::create(_["EVID"]=newEvid,
                                                _["SS"]=newSs,
                                                _["DVID"]=newDvid,
                                                _["AMT"]=newAmt,
                                                _["RATE"]=newRate,
                                                _["TINF"]=newTinf,
                                                _["CMT"]=newCmt,
                                                _["ADM"]=newAdm,
                                                _[".nlmixrKeep"]=keepItem
                                                ),
                      _["adm"]=createAdm(admIds),
                      _["turnOffCmt"]=turnOffCmt,
                      _["hasTinf"]=hasTinf,
                      _["hasRate"]=hasRate,
                      _["hasPhantom"]=hasPhantom,
                      _["hasReplace"]=hasReplace,
                      _["hasMult"]=hasMult,
                      _["hasSs"]=hasSs,
                      _["hasSs2"]=hasSs2,
                      _["hasSsRate"]=hasSsRate);
}
