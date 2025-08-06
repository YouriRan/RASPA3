module;

#ifdef USE_LEGACY_HEADERS
#include <iostream>
#include <vector>
#endif

#include "featomic.hpp"
export module property_soap;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import atom;
import simulationbox;
import featomic;

export struct PropertySoap
{
  PropertySoap(std::size_t sampleEvery, std::size_t writeOutputEvery, std::size_t numberOfBuffers, double cutOff, double smoothingWidth, double gaussianWidth,
               std::size_t numberOfRadialBasisFunctions, std::size_t numberOfAngularBasisFunctions)
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

  std::uint64_t versionNumber{0};
  std::size_t sampleEvery;
  std::size_t writeOutputEvery;
  std::size_t numberOfBuffers;
  FeatomicCalculator calculator;
  featomic::CalculationOptions options;

  std::vector<FeatomicSystem> featomicSystems;

  void sample(std::size_t currentCycle, std::span<Atom> atoms, SimulationBox box);
  void writeOutput(std::size_t currentCycle);
};