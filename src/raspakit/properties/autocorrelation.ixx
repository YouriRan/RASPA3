module;

#ifdef USE_LEGACY_HEADERS
#include <iostream>
#include <string>
#include <vector>
#endif

export module property_autocorrelation;

#ifndef USE_LEGACY_HEADERS
import <vector>;
import <string>;
#endif

import archive;

export struct PropertyAutocorrelation
{
  PropertyAutocorrelation() {};

  PropertyAutocorrelation(size_t numberOfBuffersACF, size_t bufferLengthACF, size_t sampleEvery, size_t writeEvery,
                          std::string variableName);

  uint64_t versionNumber{1};

  size_t numberOfBuffersACF;  // more buffers for better sampling
  size_t bufferLengthACF;     // tauMax
  size_t sampleEvery;
  size_t writeEvery;
  std::string variableName;

  std::vector<std::vector<double>> buffer;
  std::vector<double> accumulatedAcf;
  std::vector<int64_t> counts;  // M for every buffer
  size_t countAccumulatedACF;
  double sumValue;

  void addSample(size_t currentCycle, double value);
  void writeOutput(size_t systemId, size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyAutocorrelation &acf);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyAutocorrelation &acf);
};
