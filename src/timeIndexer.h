#ifndef __TIMEINDEXER_H__
#define __TIMEINDEXER_H__


struct timeInfo {
  int id;
  std::vector<int> indices;
};

class timeIndexer {
 public:
  timeIndexer() : initialized(false), sorted(false), nIds(0) {}

  timeIndexer(const std::vector<int>& ids,
              const std::vector<std::vector<double>>& times) {
    initialize(ids, times);
  }

  timeIndexer(const std::vector<int>& ids,
              const std::vector<double>& times) {
    initialize(ids, times);
  }

  timeIndexer(Rcpp::IntegerVector& ids,
              Rcpp::NumericVector& times) {
    initialize(ids, times);
  }

  void initialize(const std::vector<int>& ids,
                  const std::vector<std::vector<double>>& times) {
    nIds = ids.size();
    size_t k = 0;
    std::unordered_map<double, std::unordered_map<int, std::vector<int>>> timeIdIndices;

    maxTime = -INFINITY;
    for (size_t i = ids.size(); i--;) {
      int id = ids[i];
      for (size_t j = 0; j < times[i].size(); ++j) {
        double time = times[i][j];
        if (maxTime < time) {
          maxTime = time;
        }
        timeIdIndices[time][id].push_back(k);
        uniqueTimes.insert(time);
        k++;
      }
    }

    for (const auto& timeEntry : timeIdIndices) {
      double time = timeEntry.first;
      for (const auto& idEntry : timeEntry.second) {
        int id = idEntry.first;
        const std::vector<int>& indices = idEntry.second;
        timeToInfo[time].emplace_back(timeInfo{id, indices});
      }
    }
    initialized = true;
    sorted = false;
  }

  void initialize(const std::vector<int>& ids, const std::vector<double>& times,
                  const bool optTime) {
    if (optTime && initialized) {
      return;
    }
    if (isSame(ids, times)) {
      return;
    }
    std::unordered_map<int, std::vector<double>> idToTimesMap;
    reset();

    this->times = times;
    this->modelSwitch = ids;

    std::vector<int> outIds;
    std::vector<std::vector<double>> outTimes;

    for (size_t i = 0; i < times.size(); ++i) {
      idToTimesMap[ids[i]].push_back(times[i]);
    }

    for (const auto& pair : idToTimesMap) {
      outIds.push_back(pair.first);
      outTimes.push_back(pair.second);
    }

    initialize(outIds, outTimes);
  }

  void initialize(const std::vector<int>& ids, const std::vector<double>& times) {
    initialize(ids, times, false);
  }

  void initialize(Rcpp::IntegerVector& ids, Rcpp::NumericVector& time,
                  const bool optTime) {
    std::vector<int> idsVec = Rcpp::as<std::vector<int>>(ids);
    std::vector<double> timeVec = Rcpp::as<std::vector<double>>(time);
    initialize(idsVec, timeVec, optTime);
  }

  void initialize(Rcpp::IntegerVector& ids, Rcpp::NumericVector& time) {
    initialize(ids, time, false);
  }

  bool isSame(const std::vector<int>& ids, const std::vector<double>& times) {
    if (ids.size() != this->modelSwitch.size()) {
      return false;
    }
    if (times.size() != this->times.size()) {
      return false;
    }
    if (!std::equal(ids.begin(), ids.end(), this->modelSwitch.begin())) {
      return false;
    }
    if (!std::equal(times.begin(), times.end(), this->times.begin())) {
      return false;
    }
    return true;
  }

  void reset() {
    timeToInfo.clear();
    uniqueTimes.clear();
    sortedTimes.clear();
    times.clear();
    modelSwitch.clear();

    initialized = false;
    sorted = false;
    nIds=0;
  }

  size_t getNid() {
    return nIds;
  }

  double getMaxTime() {
    return maxTime;
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

  std::vector<double> getUniqueTimes() {
    if (!initialized) {
      throw std::runtime_error("timeIndexer has not been initialized");
    }
    if (!sorted) {
      return std::vector<double>(uniqueTimes.begin(), uniqueTimes.end());
    } else {
      return sortedTimes;
    }
  }

  const std::vector<timeInfo>& getTimeInfo(double time) const {
    if (!initialized) {
      throw std::runtime_error("timeIndexer has not been initialized");
    }
    try {
      return timeToInfo.at(time);
    } catch (const std::out_of_range& e) {
      REprintf("times:\n");
      Rcpp::print(Rcpp::wrap(this->times));
      REprintf("modelSwitch:");
      Rcpp::print(Rcpp::wrap(this->modelSwitch));
      Rcpp::stop("time '%f' not found in timeIndexer",
                 time);
    }
  }

  void setTimes(const std::vector<double>& times) {
    this->times = times;
  }

  void setModelSwitch(const std::vector<int>& modelSwitch) {
    this->modelSwitch = modelSwitch;
  }

  std::vector<int> getModelSwitch() {
    return this->modelSwitch;
  }

  std::vector<double> getTimes() {
    return this->times;
  }

 private:
  std::unordered_map<double, std::vector<timeInfo>> timeToInfo;
  std::set<double> uniqueTimes;
  std::vector<double> sortedTimes;
  bool initialized;
  bool sorted;
  size_t nIds;
  double maxTime;
  std::vector<double> times;
  std::vector<int> modelSwitch;
};

extern timeIndexer globalTimeIndexer;

static inline void convertToTimeIndexerStructure(const std::vector<double>& times,
                                                 const std::vector<int>& ids,
                                                 std::vector<int>& outIds,
                                                 std::vector<std::vector<double>>& outTimes) {
  if (globalTimeIndexer.isInitialized()) return;
  std::unordered_map<int, std::vector<double>> idToTimesMap;

  globalTimeIndexer.setTimes(times);
  globalTimeIndexer.setModelSwitch(ids);

  for (size_t i = 0; i < times.size(); ++i) {
    idToTimesMap[ids[i]].push_back(times[i]);
  }

  for (const auto& pair : idToTimesMap) {
    outIds.push_back(pair.first);
    outTimes.push_back(pair.second);
  }
}


#endif
