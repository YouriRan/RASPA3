module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <optional>
#include <span>
#include <tuple>
#include <type_traits>
#include <vector>
#endif

module mc_moves_swap_ncmc;

#ifndef USE_LEGACY_HEADERS
import <complex>;
import <vector>;
import <array>;
import <tuple>;
import <optional>;
import <span>;
import <optional>;
import <tuple>;
import <algorithm>;
import <chrono>;
import <cmath>;
import <iostream>;
import <iomanip>;
import <type_traits>;
#endif

import component;
import molecule;
import atom;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import randomnumbers;
import system;
import energy_factor;
import energy_status;
import energy_status_inter;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import running_energy;
import forcefield;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import mc_moves_move_types;

std::span<Atom> spanOfMolecule(size_t selectedComponent, size_t selectedMolecule, std::vector<Atom> atomPositions,
                               const System& system)
{
  return std::span(&atomPositions[system.atomsOffset(selectedComponent, selectedMolecule)],
                   system.components[selectedComponent].atoms.size());
}

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::swapNCMC(RandomNumber& random, System& system,
                                                                    size_t selectedComponent, size_t selectedMolecule,
                                                                    bool insertionDisabled, bool deletionDisabled)
{
}
/*
{
std::chrono::system_clock::time_point time_begin, time_end;
MoveTypes move = MoveTypes::SwapNCMC;
Component& component = system.components[selectedComponent];

// Retrieve lambda parameters and select a new lambda bin for the move
PropertyLambdaProbabilityHistogram& lambda = component.lambdaGC;
size_t oldBin = lambda.currentBin;
double deltaLambda = lambda.delta;
double maxChange = component.mc_moves_statistics.getMaxChange(move, 2);
std::make_signed_t<std::size_t> selectedNewBin = lambda.selectNewBin(random, maxChange);
size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];

size_t indexFractionalMolecule = system.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);

// all copied data: moleculePositions, moleculeAtomPositions, thermostat, dt
// all const data: components, forcefield, simulationbox, numberofmoleculespercomponents, fixedFrameworkStoredEik
// all scratch data: eik_x, eik_y, eik_z, eik_xy
std::span<Atom> atomPositions = system.spanOfMoleculeAtoms();
std::vector<Atom> moleculeAtomPositions(atomPositions.size());
std::copy(atomPositions.begin(), atomPositions.end(), moleculeAtomPositions.begin());

std::vector<Molecule> moleculePositions(system.moleculePositions);
std::optional<Thermostat> thermostat(system.thermostat);

std::vector<size_t> numberOfMoleculesPerComponent(system.numberOfMoleculesPerComponent);

RunningEnergy referenceEnergy = system.runningEnergies;
referenceEnergy.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(moleculePositions);
referenceEnergy.rotationalKineticEnergy =
Integrators::computeRotationalKineticEnergy(moleculePositions, system.components);
RunningEnergy currentEnergy = referenceEnergy;

// get Timestep from the max change
double dt = system.timeStep;
size_t numberOfPerturbations = lambda.numberOfSamplePoints - 1;
size_t numberOfMDStepsPerIntegration = 10;

bool insert = (random.uniform() < 0.5);

if (insert)  // Insertion move
{
for (size_t i = 0; i < numberOfPerturbations; i++)
{
size_t newBin = (oldBin + i + 1) % numberOfPerturbations;
size_t newLambda = deltaLambda * static_cast<double>(newBin);

std::span<Atom> fractionalMolecule =
spanOfMolecule(selectedComponent, indexFractionalMolecule, atomPositions, system);
std::vector<Atom> oldFractionalMolecule(fractionalMolecule.begin(), fractionalMolecule.end());

if (newBin == 0)
{
for (Atom& atom : fractionalMolecule)
{
atom.setScalingToInteger();
}
}
else if (newBin == 1)
{
// insert lambda=(1/bins) particle
size_t newMolecule = system.numberOfMoleculesPerComponent[selectedComponent];

time_begin = std::chrono::system_clock::now();
std::optional<ChainData> growData = CBMC::growMoleculeSwapInsertion(
random, system.frameworkComponents, component, system.hasExternalField, system.components,
system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), moleculeAtomPositions, system.beta,
growType, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent, newMolecule, 0.0,
static_cast<size_t>(oldFractionalMolecule.front().groupId), system.numberOfTrialDirections);
time_end = std::chrono::system_clock::now();
component.mc_moves_cputime[move]["Insertion-NonEwald"] += (time_end - time_begin);
system.mc_moves_cputime[move]["Insertion-NonEwald"] += (time_end - time_begin);

std::vector<Atom>::const_iterator iterator = system.iteratorForMolecule(
selectedComponent, selectedMolecule, components, atomPositions, numberOfFrameworkAtoms);
moleculeAtomPositions.insert(iterator, growData->atom.begin(), growData->atom.end());

std::vector<Molecule>::iterator moleculeIterator = system.indexForMolecule(
selectedComponent, selectedMolecule, components, atomPositions, numberOfFrameworkAtoms);
moleculePositions.insert(moleculeIterator, growData->molecule);

numberOfMoleculesPerComponent[selectedComponent]++;

// size_t lastMoleculeId = system.numberOfMoleculesPerComponent[selectedComponent] - 1;
// std::span<Atom> lastMolecule = system.spanOfMolecule(selectedComponent, lastMoleculeId);
// fractionalMolecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);
// std::swap_ranges(fractionalMolecule.begin(), fractionalMolecule.end(), lastMolecule.begin());
// std::swap(system.moleculePositions[system.moleculeIndexOfComponent(selectedComponent, indexFractionalMolecule)],
//           system.moleculePositions[system.moleculeIndexOfComponent(selectedComponent, lastMoleculeId)]);

}
else
{
for (Atom& atom : fractionalMolecule)
{
atom.setScaling(newLambda);
}
}

currentEnergy = Integrators::updateGradients(
moleculeAtomPositions, system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox,
system.components, system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik,
system.fixedFrameworkStoredEik, numberOfMoleculesPerComponent);

for (size_t step = 0; step < numberOfMDStepsPerIntegration; step++)
{
currentEnergy = Integrators::velocityVerlet(
moleculePositions, moleculeAtomPositions, system.components, dt, thermostat, system.spanOfFrameworkAtoms(),
system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
system.totalEik, system.fixedFrameworkStoredEik,  system.interpolationGrids, system.numberOfMoleculesPerComponent);
}
}

if ((system.insideBlockedPockets(component, fractionalMolecule)))
{
// Reject, set fractional molecule back to old state
std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
return {std::nullopt, double3(0.0, 1.0, 0.0)};
}

double drift = currentEnergy.conservedEnergy() - referenceEnergy.conservedEnergy();

// Compute acceptance probability
double fugacity = component.fugacityCoefficient.value_or(1.0) * system.pressure;
double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
double preFactor = system.beta * component.molFraction * fugacity * system.simulationBox.volume /
static_cast<double>(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);
double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
double Pacc =
preFactor * (growData->RosenbluthWeight / idealGasRosenbluthWeight) * exp(-system.beta * drift + biasTerm);

// Retrieve bias from transition matrix
double biasTransitionMatrix = system.tmmc.biasFactor(oldN + 1, oldN);

if (system.tmmc.doTMMC)
{
size_t newN = oldN + 1;
if (newN > system.tmmc.maxMacrostate)
{
return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
}
}

if (random.uniform() < biasTransitionMatrix * Pacc)
{
Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
system.insertMolecule(selectedComponent, growData->)
}


}

else  // Deletion move
{
if (system.numberOfIntegerMoleculesPerComponent[selectedComponent] == 0)
{
return {std::nullopt, double3(0.0, 0.0, 0.0)};
}


for (size_t i = 0; i < numberOfPerturbations; i++)
{
size_t newBin = (oldBin - i - 1) % numberOfPerturbations;
size_t newLambda = deltaLambda * static_cast<double>(newBin);

std::span<Atom> fractionalMolecule =
spanOfMolecule(selectedComponent, indexFractionalMolecule, atomPositions, system);
std::vector<Atom> oldFractionalMolecule(fractionalMolecule.begin(), fractionalMolecule.end());

if (newBin == 0)
{
for (Atom& atom : fractionalMolecule)
{
atom.setScalingToInteger();
}
}
else if (newBin == 1)
{
// insert lambda=(1/bins) particle
size_t newMolecule = system.numberOfMoleculesPerComponent[selectedComponent];

time_begin = std::chrono::system_clock::now();
std::optional<ChainData> growData = CBMC::growMoleculeSwapInsertion(
random, system.frameworkComponents, component, system.hasExternalField, system.components,
system.forceField, system.simulationBox, system.spanOfFrameworkAtoms(), moleculeAtomPositions, system.beta,
growType, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent, newMolecule, 0.0,
static_cast<size_t>(oldFractionalMolecule.front().groupId), system.numberOfTrialDirections);
time_end = std::chrono::system_clock::now();
component.mc_moves_cputime[move]["Insertion-NonEwald"] += (time_end - time_begin);
system.mc_moves_cputime[move]["Insertion-NonEwald"] += (time_end - time_begin);

std::vector<Atom>::const_iterator iterator =
iteratorForMolecule(selectedComponent, selectedMolecule, components, atomPositions, numberOfFrameworkAtoms);
moleculeAtomPositions.insert(iterator, growData->atom.begin(), growData->atom.end());

std::vector<Molecule>::iterator moleculeIterator =
indexForMolecule(selectedComponent, selectedMolecule, components, atomPositions, numberOfFrameworkAtoms);
moleculePositions.insert(moleculeIterator, growData->molecule);

numberOfMoleculesPerComponent[selectedComponent]++;

// size_t lastMoleculeId = system.numberOfMoleculesPerComponent[selectedComponent] - 1;
// std::span<Atom> lastMolecule = system.spanOfMolecule(selectedComponent, lastMoleculeId);
// fractionalMolecule = system.spanOfMolecule(selectedComponent, indexFractionalMolecule);
// std::swap_ranges(fractionalMolecule.begin(), fractionalMolecule.end(), lastMolecule.begin());
// std::swap(system.moleculePositions[system.moleculeIndexOfComponent(selectedComponent, indexFractionalMolecule)],
//           system.moleculePositions[system.moleculeIndexOfComponent(selectedComponent, lastMoleculeId)]);

}
else
{
for (Atom& atom : fractionalMolecule)
{
atom.setScaling(newLambda);
}
}

currentEnergy = Integrators::updateGradients(
moleculeAtomPositions, system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox,
system.components, system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik,
system.fixedFrameworkStoredEik, numberOfMoleculesPerComponent);

for (size_t step = 0; step < numberOfMDStepsPerIntegration; step++)
{
currentEnergy = Integrators::velocityVerlet(
moleculePositions, moleculeAtomPositions, system.components, dt, thermostat, system.spanOfFrameworkAtoms(),
system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
system.totalEik, system.fixedFrameworkStoredEik, system.interpolationGrids, system.numberOfMoleculesPerComponent);
}
}

if ((system.insideBlockedPockets(component, fractionalMolecule)))
{
// Reject, set fractional molecule back to old state
std::copy(oldFractionalMolecule.begin(), oldFractionalMolecule.end(), fractionalMolecule.begin());
return {std::nullopt, double3(0.0, 1.0, 0.0)};
}

double drift = currentEnergy.conservedEnergy() - referenceEnergy.conservedEnergy();

// Compute acceptance probability
double fugacity = component.fugacityCoefficient.value_or(1.0) * system.pressure;
double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
double preFactor = system.beta * component.molFraction * fugacity * system.simulationBox.volume /
static_cast<double>(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);
double biasTerm = lambda.biasFactor[newBin] - lambda.biasFactor[oldBin];
double Pacc =
preFactor * (growData->RosenbluthWeight / idealGasRosenbluthWeight) * exp(-system.beta * drift + biasTerm);

// Retrieve bias from transition matrix
double biasTransitionMatrix = system.tmmc.biasFactor(oldN + 1, oldN);

if (system.tmmc.doTMMC)
{
size_t newN = oldN + 1;
if (newN > system.tmmc.maxMacrostate)
{
return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
}
}

if (random.uniform() < biasTransitionMatrix * Pacc)
{
Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
system.insertMolecule(selectedComponent, growData->)
}




}



return {std::nullopt, double3(0.0, 1.0, 0.0)};
}
*/