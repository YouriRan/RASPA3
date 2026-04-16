module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <optional>
#include <print>
#include <span>
#include <tuple>
#include <vector>
#endif

module mc_moves_noneq_cbmc;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import component;
import molecule;
import atom;
import double3;
import double3x3;
import simd_quatd;
import simulationbox;
import cbmc;
import cbmc_chain_data;
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
import transition_matrix;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import interactions_external_field;
import mc_moves_move_types;
import thermostat;
import integrators;
import integrators_update;
import integrators_compute;
import units;

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::NonEqCBMC(RandomNumber& random, System& system,
                                                                     std::size_t selectedComponent)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::SwapNonEqCBMC;
  Component& component = system.components[selectedComponent];

  std::vector<std::size_t> numberOfMoleculesPerComponent(system.numberOfMoleculesPerComponent);
  std::size_t oldN = numberOfMoleculesPerComponent[selectedComponent];

  // all copied data: moleculeData, moleculeAtomPositions, thermostat, dt
  // all const data: components, forcefield, simulationbox, numberofmoleculespercomponents, fixedFrameworkStoredEik
  // all scratch data: eik_x, eik_y, eik_z, eik_xy
  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::vector<Atom> moleculeAtomPositions(atomData.size());
  std::copy(atomData.begin(), atomData.end(), moleculeAtomPositions.begin());

  std::vector<Molecule> moleculeData(system.moleculeData);
  std::optional<Thermostat> thermostat(system.thermostat);

  // get Timestep from the max change
  double dt = system.timeStep;

  bool insert = (random.uniform() < 0.5);

  double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
  double cutOffCoulomb = system.forceField.cutOffCoulomb;
  Component::GrowType growType = component.growType;

  auto computeTemperature = [](const System& system, double translationalKineticEnergy) {
      return 2.0 * translationalKineticEnergy / (Units::KB * static_cast<double>(
          system.translationalDegreesOfFreedom - system.translationalCenterOfMassConstraint
      ));
  };

  std::ofstream stream("ncmc.txt", std::ios::app);

  // insertion / deletion without acceptance
  if (insert)  // Insertion
  {
    component.mc_moves_statistics.addTrial(move, 0);

    // Attempt to grow a new molecule using CBMC
    time_begin = std::chrono::system_clock::now();
    std::pair<Molecule, std::vector<Atom>> trialMolecule =
        component.equilibratedMoleculeRandomInBox(random, system.simulationBox);
    time_end = std::chrono::system_clock::now();

    // Check if the new molecule is inside blocked pockets
    if (system.insideBlockedPockets(system.components[selectedComponent], trialMolecule.second))
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    std::for_each(std::begin(trialMolecule.second), std::end(trialMolecule.second),
                  [selectedComponent, oldN](Atom& atom)
                  {
                    atom.moleculeId = static_cast<std::uint32_t>(oldN);
                    atom.componentId = static_cast<std::uint8_t>(selectedComponent);
                    atom.groupId = static_cast<std::uint8_t>(0);
                    atom.setScaling(1.0);
                  });

    // compute framework-molecule energy contribution
    std::optional<RunningEnergy> frameworkMolecule;
      frameworkMolecule = Interactions::computeFrameworkMoleculeEnergyDifference(
          system.forceField, system.simulationBox, system.interpolationGrids, system.framework,
          system.spanOfFrameworkAtoms(), trialMolecule.second, {});
    if (!frameworkMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // compute molecule-molecule energy contribution
    std::optional<RunningEnergy> interMolecule = Interactions::computeInterMolecularEnergyDifference(
        system.forceField, system.simulationBox, system.spanOfMoleculeAtoms(), trialMolecule.second, {});
    if (!interMolecule.has_value()) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    // insert molecule into copied atom positions
    std::size_t index = 0;
    for (std::size_t comp = 0; comp < selectedComponent + 1; ++comp)
    {
      index += system.components[comp].atoms.size() * system.numberOfMoleculesPerComponent[comp];
    }

    std::vector<Atom>::const_iterator iterator =
        moleculeAtomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(index);
    moleculeAtomPositions.insert(iterator, trialMolecule.second.begin(), trialMolecule.second.end());

    index = std::accumulate(system.numberOfMoleculesPerComponent.begin(),
                            system.numberOfMoleculesPerComponent.begin() + selectedComponent + 1, 0uz);

    std::vector<Molecule>::iterator moleculeIterator =
        moleculeData.begin() + static_cast<std::vector<Atom>::difference_type>(index);
    moleculeData.insert(moleculeIterator, trialMolecule.first);
    
    numberOfMoleculesPerComponent[selectedComponent] += 1;

    // MD INTEGRATION
    // initialize the velocities according to Boltzmann distribution
    // NOTE: it is important that the reference energy has the initial kinetic energies
    Integrators::initializeVelocities(random, moleculeData, system.components, system.temperature);

    if (system.numberOfFrameworkAtoms == 0)
    {
      Integrators::removeCenterOfMassVelocityDrift(moleculeData);
    }

    // before getting energy, recompute current energy
    RunningEnergy intermediateEnergy = Integrators::updateGradients(
        moleculeAtomPositions, system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox,
        system.components, system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik,
        system.fixedFrameworkStoredEik, system.interpolationGrids, numberOfMoleculesPerComponent);
    
    if (intermediateEnergy.potentialEnergy() > system.forceField.energyOverlapCriteria) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    intermediateEnergy.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(moleculeData);
    intermediateEnergy.rotationalKineticEnergy =
        Integrators::computeRotationalKineticEnergy(moleculeData, system.components);
    RunningEnergy currentEnergy = intermediateEnergy;

    // integrate for N steps
    time_begin = std::chrono::system_clock::now();
    // std::size_t idx = 0;
    // std::print(stream, "\ninsertion:\n");
    // std::print(stream, "{}: {} {} {}\n", ++idx, system.runningEnergies.potentialEnergy(), system.runningEnergies.translationalKineticEnergy, computeTemperature(system, system.runningEnergies.translationalKineticEnergy));
    // std::print(stream, "{}: {} {} {}\n", ++idx, currentEnergy.potentialEnergy(), currentEnergy.translationalKineticEnergy, computeTemperature(system, currentEnergy.translationalKineticEnergy));
    for (std::size_t step = 0; step < system.numberOfHybridMCSteps; ++step)
    {
      currentEnergy = Integrators::velocityVerlet(moleculeData, moleculeAtomPositions, system.components, dt,
                                                  thermostat, system.spanOfFrameworkAtoms(), system.forceField,
                                                  system.simulationBox, system.eik_x, system.eik_y, system.eik_z,
                                                  system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
                                                  system.interpolationGrids, numberOfMoleculesPerComponent);

    // std::print(stream, "{}: {} {} {}\n", ++idx, currentEnergy.potentialEnergy(), currentEnergy.translationalKineticEnergy, computeTemperature(system, currentEnergy.translationalKineticEnergy));
    }
    time_end = std::chrono::system_clock::now();

    if (currentEnergy.potentialEnergy() > system.forceField.energyOverlapCriteria) return {std::nullopt, double3(0.0, 1.0, 0.0)};
    if (std::isnan(currentEnergy.potentialEnergy()))  return {std::nullopt, double3(0.0, 1.0, 0.0)};

    system.mc_moves_cputime[move]["Integration"] += (time_end - time_begin);

    // Calculate correction factor for Ewald energy difference
    double correctionFactorEwald = 1.0;
    // std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy()));

    // Compute the acceptance probability pre-factor
    double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
    double preFactor = correctionFactorEwald * system.beta * fugacity * system.simulationBox.volume /
                       double(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);

    double work = currentEnergy.potentialEnergy() - system.runningEnergies.potentialEnergy();
    double drift = std::abs(currentEnergy.conservedEnergy() - intermediateEnergy.conservedEnergy());
    double Pacc = preFactor * std::exp(-system.beta * (work + drift));
    
    std::cout << "insertion: " << work << " " << drift << " " << Pacc << std::endl;
    if (std::isnan(work) || std::isinf(work))  return {std::nullopt, double3(0.0, 1.0, 0.0)};
    if (std::isnan(drift) || std::isinf(drift))  return {std::nullopt, double3(0.0, 1.0, 0.0)};

    component.mc_moves_statistics.addConstructed(move, 0);

    if (random.uniform() < Pacc)
    {
      component.mc_moves_statistics.addAccepted(move, 0);

      system.insertMolecule(selectedComponent, trialMolecule.first, trialMolecule.second);
      std::copy(moleculeData.begin(), moleculeData.end(), system.moleculeData.begin());

      system.thermostat = thermostat;
      system.timeStep = dt;

      std::span<Atom> newAtomData = system.spanOfMoleculeAtoms();
      std::copy(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), newAtomData.begin());
      // system.spanOfMoleculeAtoms() = moleculeAtomPositions;

      system.updateMoleculeAtomInformation();
      Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components);
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      return {currentEnergy, double3(0.0, 1.0 - Pacc, Pacc)};
    }
    return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
  }
  else
  {
    component.mc_moves_statistics.addTrial(move, 1);
    if (oldN == 0)
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Get a reference to the molecule being deleted
    std::size_t selectedMolecule = system.randomIntegerMoleculeOfComponent(random, selectedComponent);
    std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

    // delete molecule from copied atom positions
    std::size_t index = 0;
    for (std::size_t comp = 0; comp < selectedComponent; ++comp)
    {
      index += system.components[comp].atoms.size() * system.numberOfMoleculesPerComponent[comp];
    }
    index += selectedMolecule * system.components[selectedComponent].atoms.size();

    std::vector<Atom>::const_iterator iterator =
        moleculeAtomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(index);
    moleculeAtomPositions.erase(iterator, iterator + static_cast<std::vector<Atom>::difference_type>(molecule.size()));

    index = std::accumulate(system.numberOfMoleculesPerComponent.begin(),
                            system.numberOfMoleculesPerComponent.begin() + selectedComponent, 0uz);
    index += selectedMolecule;

    std::vector<Molecule>::iterator moleculeIterator =
        moleculeData.begin() + static_cast<std::vector<Atom>::difference_type>(index);
    moleculeData.erase(moleculeIterator, moleculeIterator + 1);
    
    numberOfMoleculesPerComponent[selectedComponent] -= 1;

    // MD INTEGRATION

    // initialize the velocities according to Boltzmann distribution
    // NOTE: it is important that the reference energy has the initial kinetic energies
    Integrators::initializeVelocities(random, moleculeData, system.components, system.temperature);

    if (system.numberOfFrameworkAtoms == 0)
    {
      Integrators::removeCenterOfMassVelocityDrift(moleculeData);
    }

    RunningEnergy intermediateEnergy = Integrators::updateGradients(
        moleculeAtomPositions, system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox,
        system.components, system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik,
        system.fixedFrameworkStoredEik, system.interpolationGrids, numberOfMoleculesPerComponent);
        
        if (intermediateEnergy.potentialEnergy() > system.forceField.energyOverlapCriteria) return {std::nullopt, double3(0.0, 1.0, 0.0)};

    intermediateEnergy.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(moleculeData);
    intermediateEnergy.rotationalKineticEnergy =
        Integrators::computeRotationalKineticEnergy(moleculeData, system.components);
    RunningEnergy currentEnergy = intermediateEnergy;

    // integrate for N steps
    time_begin = std::chrono::system_clock::now();
    std::size_t idx = 0;
    std::print(stream, "\ndeletion:\n");
    std::print(stream, "{}: {} {} {}\n", ++idx, system.runningEnergies.potentialEnergy(), system.runningEnergies.translationalKineticEnergy, computeTemperature(system, system.runningEnergies.translationalKineticEnergy));
    std::print(stream, "{}: {} {} {}\n", ++idx, currentEnergy.potentialEnergy(), currentEnergy.translationalKineticEnergy, computeTemperature(system, currentEnergy.translationalKineticEnergy));
    for (std::size_t step = 0; step < system.numberOfHybridMCSteps; ++step)
    {
      currentEnergy = Integrators::velocityVerlet(moleculeData, moleculeAtomPositions, system.components, dt,
                                                  thermostat, system.spanOfFrameworkAtoms(), system.forceField,
                                                  system.simulationBox, system.eik_x, system.eik_y, system.eik_z,
                                                  system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
                                                  system.interpolationGrids, numberOfMoleculesPerComponent);
    std::print(stream, "{}: {} {} {}\n", ++idx, currentEnergy.potentialEnergy(), currentEnergy.translationalKineticEnergy, computeTemperature(system, currentEnergy.translationalKineticEnergy));
    }
    time_end = std::chrono::system_clock::now();

    if (currentEnergy.potentialEnergy() > system.forceField.energyOverlapCriteria) return {std::nullopt, double3(0.0, 1.0, 0.0)};
    if (std::isnan(currentEnergy.potentialEnergy()))  return {std::nullopt, double3(0.0, 1.0, 0.0)};
    // Calculate the correction factor for Ewald summation
    double correctionFactorEwald = 1.0;
    // std::exp(-system.beta * (energyFourierDifference.potentialEnergy) + tailEnergyDifference.potentialEnergy()));

    // Compute acceptance probability factors
    double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
    double preFactor = correctionFactorEwald * double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                       (system.beta * component.molFraction * fugacity * system.simulationBox.volume);

    system.mc_moves_cputime[move]["Integration"] += (time_end - time_begin);
    component.mc_moves_statistics.addConstructed(move, 1);

    double work = currentEnergy.potentialEnergy() - system.runningEnergies.potentialEnergy();
    double drift = std::abs(currentEnergy.conservedEnergy() - intermediateEnergy.conservedEnergy());

    double Pacc = preFactor * std::exp(-system.beta * (work + drift));
    std::cout << "deletion: " << work << " " << drift << " " << Pacc <<  std::endl;

    if (std::isnan(work) || std::isinf(work))  return {std::nullopt, double3(0.0, 1.0, 0.0)};
    if (std::isnan(drift) || std::isinf(drift))  return {std::nullopt, double3(0.0, 1.0, 0.0)};

    if (random.uniform() < Pacc)
    {
      component.mc_moves_statistics.addAccepted(move, 1);

      system.deleteMolecule(selectedComponent, selectedMolecule, molecule);
      std::copy(moleculeData.begin(), moleculeData.end(), system.moleculeData.begin());
      system.thermostat = thermostat;
      system.timeStep = dt;

      std::span<Atom> newAtomData = system.spanOfMoleculeAtoms();
      std::copy(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), newAtomData.begin());
      // system.spanOfMoleculeAtoms() = moleculeAtomPositions;

      system.updateMoleculeAtomInformation();
      Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components);
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      return {currentEnergy, double3(Pacc, 1.0 - Pacc, 0.0)};
    }

    return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
  }
}
