module;

#ifdef USE_LEGACY_HEADERS
#include <iostream>
#include <string>
#include <vector>
#endif

export module property_autocorrelation;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import archive;

export struct PropertyAutocorrelation
{
  PropertyAutocorrelation() {};

  PropertyAutocorrelation(std::size_t numberOfBuffersACF, std::size_t bufferLengthACF, std::size_t sampleEvery, std::size_t writeEvery,
                          std::string variableName);

  std::uint64_t versionNumber{1};

  std::size_t numberOfBuffersACF;  // more buffers for better sampling
  std::size_t bufferLengthACF;     // tauMax
  std::size_t sampleEvery;
  std::size_t writeEvery;
  std::string variableName;

  std::vector<std::vector<double>> buffer;
  std::vector<double> accumulatedAcf;
  std::vector<std::int64_t> counts;  // M for every buffer
  std::size_t countAccumulatedACF;
  double sumValue;

  void addSample(std::size_t currentCycle, double value);
  void writeOutput(std::size_t systemId, std::size_t currentCycle);

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const PropertyAutocorrelation &acf);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, PropertyAutocorrelation &acf);
};
