module;

#ifdef USE_LEGACY_HEADERS
#include <iostream>
#include <vector>
#endif

#include "featomic.hpp"
export module property_soap;

#ifndef USE_LEGACY_HEADERS
import <iostream>;
import <vector>;
#endif

import atom;
import simulationbox;
import featomic;

export struct PropertySoap
{
  PropertySoap(size_t sampleEvery, size_t writeOutputEvery, size_t numberOfBuffers, double cutOff, double smoothingWidth, double gaussianWidth,
               size_t numberOfRadialBasisFunctions, size_t numberOfAngularBasisFunctions)
      : sampleEvery(sampleEvery),
      writeOutputEvery(writeOutputEvery),
      numberOfBuffers(numberOfBuffers),
      calculator(getSoapPowerSpectrumCalculator(cutOff, smoothingWidth, gaussianWidth, numberOfRadialBasisFunctions, numberOfAngularBasisFunctions
      )),
      options(),
      featomicSystems(numberOfBuffers)
  {
    options.use_native_system = true;
    // options.gradients.push_back("positions");
    std::cout << "Creating property " << std::endl;
  }

  uint64_t versionNumber{0};
  size_t sampleEvery;
  size_t writeOutputEvery;
  size_t numberOfBuffers;
  FeatomicCalculator calculator;
  featomic::CalculationOptions options;

  std::vector<FeatomicSystem> featomicSystems;

  void sample(size_t currentCycle, std::span<Atom> atoms, SimulationBox box);
  void writeOutput(size_t currentCycle);
};