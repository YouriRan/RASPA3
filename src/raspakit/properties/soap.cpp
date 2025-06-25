module;

#ifdef USE_LEGACY_HEADERS
#include <iostream>
#include <vector>
#endif

module property_soap;

#ifndef USE_LEGACY_HEADERS
import <iostream>;
import <vector>;
#endif

import atom;
import simulationbox;
import featomic;

void PropertySoap::sample(size_t currentCycle, std::span<Atom> atoms, SimulationBox box) 
{
    if (currentCycle % sampleEvery != 0uz) return;

    std::cout << "Sampling ..." << std::endl;
    featomicSystems.push_back(FeatomicSystem(atoms, box));
}

void PropertySoap::writeOutput(size_t currentCycle) 
{
    if (currentCycle % writeOutputEvery != 0uz) return;
    std::cout << featomicSystems.size() << std::endl;
    auto tsm = calculator.compute(featomicSystems, options);
    tsm.save(std::format("out_{}.mts", currentCycle));

}