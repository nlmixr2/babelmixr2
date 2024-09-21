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

  void initialize(const std::vector<int>& ids, const std::vector<std::vector<double>>& times) {
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

  void reset() {
    timeToInfo.clear();
    uniqueTimes.clear();
    sortedTimes.clear();
    times.clear();

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
    return timeToInfo.at(time);
  }

  void setTimes(const std::vector<double>& times) {
    this->times = times;
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
};

extern timeIndexer globalTimeIndexer;

static inline void convertToTimeIndexerStructure(const std::vector<double>& times,
                                                 const std::vector<int>& ids,
                                                 std::vector<int>& outIds,
                                                 std::vector<std::vector<double>>& outTimes) {
  if (globalTimeIndexer.isInitialized()) return;
  std::unordered_map<int, std::vector<double>> idToTimesMap;

  for (size_t i = 0; i < times.size(); ++i) {
    idToTimesMap[ids[i]].push_back(times[i]);
  }

  for (const auto& pair : idToTimesMap) {
    outIds.push_back(pair.first);
    outTimes.push_back(pair.second);
  }
}
#endif
