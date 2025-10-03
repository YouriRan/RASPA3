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

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::NonEqCBMC(RandomNumber& random, System& system,
                                                                     std::size_t selectedComponent)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::SwapNonEqCBMC;
  Component& component = system.components[selectedComponent];
  std::size_t oldN = system.numberOfMoleculesPerComponent[selectedComponent];

  // Update move counts statistics for swap insertion move

  // all copied data: moleculeData, moleculeAtomPositions, thermostat, dt
  // all const data: components, forcefield, simulationbox, numberofmoleculespercomponents, fixedFrameworkStoredEik
  // all scratch data: eik_x, eik_y, eik_z, eik_xy
  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::vector<Atom> moleculeAtomPositions(atomData.size());
  std::copy(atomData.begin(), atomData.end(), moleculeAtomPositions.begin());

  std::vector<Molecule> moleculeData(system.moleculeData);
  std::optional<Thermostat> thermostat(system.thermostat);

  // get Timestep from the max change
  double dt = component.mc_moves_statistics.getMaxChange(move, 0);

  bool insert = (random.uniform() < 0.5);

  double cutOffFrameworkVDW = system.forceField.cutOffFrameworkVDW;
  double cutOffMoleculeVDW = system.forceField.cutOffMoleculeVDW;
  double cutOffCoulomb = system.forceField.cutOffCoulomb;
  Component::GrowType growType = component.growType;

  // insertion / deletion without acceptance
  if (insert)  // Insertion
  {
    component.mc_moves_statistics.addTrial(move, 0);
    // Attempt to grow a new molecule using CBMC
    time_begin = std::chrono::system_clock::now();
    std::optional<ChainGrowData> growData = CBMC::growMoleculeSwapInsertion(
        random, component, system.hasExternalField, system.forceField, system.simulationBox,
        system.interpolationGrids, system.framework, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
        system.beta, growType, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, oldN, 1.0, false, false);
    time_end = std::chrono::system_clock::now();

    // Update CPU time statistics for the non-Ewald part of the move
    component.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);

    // If growth failed, reject the move
    if (!growData)
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    std::span<const Atom> newMolecule = std::span(growData->atom.begin(), growData->atom.end());

    // Check if the new molecule is inside blocked pockets
    if (system.insideBlockedPockets(system.components[selectedComponent], newMolecule))
    {
      std::print("Blocka\n");
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }

    // Compute energy difference due to Ewald Fourier components
    time_begin = std::chrono::system_clock::now();
    RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, newMolecule, {});

    time_end = std::chrono::system_clock::now();

    // Update CPU time statistics for the Ewald part of the move
    component.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

    // Compute tail energy difference due to long-range corrections
    time_begin = std::chrono::system_clock::now();
    RunningEnergy tailEnergyDifference =
        Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                                system.spanOfMoleculeAtoms(), newMolecule, {}) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                   system.spanOfFrameworkAtoms(), newMolecule, {});
    time_end = std::chrono::system_clock::now();

    // Update CPU time statistics for the tail corrections
    component.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);

    // insert molecule into copied atom positions
    std::size_t index = 0;
    for (std::size_t comp = 0; comp < selectedComponent + 1; ++comp)
    {
      index += system.components[comp].atoms.size() * system.numberOfMoleculesPerComponent[comp];
    }

    std::vector<Atom>::const_iterator iterator =
        moleculeAtomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(index);
    moleculeAtomPositions.insert(iterator, growData->atom.begin(), growData->atom.end());

    index = std::accumulate(system.numberOfMoleculesPerComponent.begin(),
                            system.numberOfMoleculesPerComponent.begin() + selectedComponent + 1, 0uz);

    std::vector<Molecule>::iterator moleculeIterator =
        moleculeData.begin() + static_cast<std::vector<Atom>::difference_type>(index);
    moleculeData.insert(moleculeIterator, growData->molecule);

    // MD INTEGRATION
    // initialize the velocities according to Boltzmann distribution
    // NOTE: it is important that the reference energy has the initial kinetic energies
    Integrators::initializeVelocities(random, moleculeData, system.components, system.temperature);

    if (system.numberOfFrameworkAtoms == 0)
    {
      Integrators::removeCenterOfMassVelocityDrift(moleculeData);
    }

    // before getting energy, recompute current energy
    system.precomputeTotalGradients();

    RunningEnergy referenceEnergy = system.runningEnergies;

    referenceEnergy.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(moleculeData);
    referenceEnergy.rotationalKineticEnergy =
        Integrators::computeRotationalKineticEnergy(moleculeData, system.components);
    RunningEnergy currentEnergy = referenceEnergy;

    // integrate for N steps
    time_begin = std::chrono::system_clock::now();
    for (std::size_t step = 0; step < system.numberOfHybridMCSteps; ++step)
    {
      currentEnergy = Integrators::velocityVerlet(moleculeData, moleculeAtomPositions, system.components, dt,
                                                  thermostat, system.spanOfFrameworkAtoms(), system.forceField,
                                                  system.simulationBox, system.eik_x, system.eik_y, system.eik_z,
                                                  system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
                                                  system.interpolationGrids, system.numberOfMoleculesPerComponent);
    }
    time_end = std::chrono::system_clock::now();

    system.mc_moves_cputime[move]["Integration"] += (time_end - time_begin);

    // Calculate correction factor for Ewald energy difference
    double correctionFactorEwald =
        std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy()));

    // Compute the acceptance probability pre-factor
    double fugacity = component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
    double preFactor = correctionFactorEwald * system.beta * component.molFraction * fugacity *
                       system.simulationBox.volume /
                       double(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);

    // Calculate the acceptance probability Pacc
    double drift = currentEnergy.potentialEnergy() - referenceEnergy.potentialEnergy();
    double Pacc = preFactor;  // * (growData->RosenbluthWeight / idealGasRosenbluthWeight);
    double biasTransitionMatrix = system.tmmc.biasFactor(oldN + 1, oldN);

    component.mc_moves_statistics.addConstructed(move, 0);
    // Check if TMMC is enabled and macrostate limit is not exceeded
    if (system.tmmc.doTMMC)
    {
      std::size_t newN = oldN + 1;
      if (newN > system.tmmc.maxMacrostate)
      {
        return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
      }
    }

    if (random.uniform() < biasTransitionMatrix * Pacc * std::exp(-system.beta * drift))
    {
      component.mc_moves_statistics.addAccepted(move, 0);

      system.insertMolecule(selectedComponent, growData->molecule, growData->atom);

      system.moleculeData = moleculeData;
      system.thermostat = thermostat;
      system.timeStep = dt;

      std::copy(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), atomData.begin());
      system.spanOfMoleculeAtoms() = moleculeAtomPositions;

      Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components);
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      return {currentEnergy, double3(0.0, 1.0 - Pacc, Pacc)};
    }
    return {std::nullopt, double3(0.0, 1.0 - Pacc, Pacc)};
  }
  else
  {
    if (oldN == 0)
    {
      return {std::nullopt, double3(0.0, 1.0, 0.0)};
    }
    component.mc_moves_statistics.addTrial(move, 1);

    // Get a reference to the molecule being deleted
    std::size_t selectedMolecule = system.randomIntegerMoleculeOfComponent(random, selectedComponent);
    std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);

    // MD INTEGRATION

    // initialize the velocities according to Boltzmann distribution
    // NOTE: it is important that the reference energy has the initial kinetic energies
    Integrators::initializeVelocities(random, moleculeData, system.components, system.temperature);

    if (system.numberOfFrameworkAtoms == 0)
    {
      Integrators::removeCenterOfMassVelocityDrift(moleculeData);
    }

    // before getting energy, recompute current energy
    system.precomputeTotalGradients();

    RunningEnergy referenceEnergy = system.runningEnergies;
    referenceEnergy.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(moleculeData);
    referenceEnergy.rotationalKineticEnergy =
        Integrators::computeRotationalKineticEnergy(moleculeData, system.components);
    RunningEnergy currentEnergy = referenceEnergy;

    // integrate for N steps
    time_begin = std::chrono::system_clock::now();
    for (std::size_t step = 0; step < system.numberOfHybridMCSteps; ++step)
    {
      currentEnergy = Integrators::velocityVerlet(moleculeData, moleculeAtomPositions, system.components, dt,
                                                  thermostat, system.spanOfFrameworkAtoms(), system.forceField,
                                                  system.simulationBox, system.eik_x, system.eik_y, system.eik_z,
                                                  system.eik_xy, system.totalEik, system.fixedFrameworkStoredEik,
                                                  system.interpolationGrids, system.numberOfMoleculesPerComponent);
    }
    time_end = std::chrono::system_clock::now();

    // Retrace the molecule for the swap deletion using CBMC algorithm
    time_begin = std::chrono::system_clock::now();
    ChainRetraceData retraceData = CBMC::retraceMoleculeSwapDeletion(
        random, component, system.hasExternalField, system.forceField, system.simulationBox,
        system.interpolationGrids, system.framework, system.spanOfFrameworkAtoms(), system.spanOfMoleculeAtoms(),
        system.beta, growType, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, molecule);
    time_end = std::chrono::system_clock::now();

    // Update the CPU time statistics for the non-Ewald part of the move
    component.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["NonEwald"] += (time_end - time_begin);

    // Compute the energy difference in Fourier space due to the deletion
    time_begin = std::chrono::system_clock::now();
    RunningEnergy energyFourierDifference = Interactions::energyDifferenceEwaldFourier(
        system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.storedEik, system.totalEik, system.forceField,
        system.simulationBox, {}, molecule);
    time_end = std::chrono::system_clock::now();
    // Update the CPU time statistics for the Ewald part of the move
    component.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Ewald"] += (time_end - time_begin);

    // Compute the tail energy difference due to the deletion
    time_begin = std::chrono::system_clock::now();
    [[maybe_unused]] RunningEnergy tailEnergyDifference =
        Interactions::computeInterMolecularTailEnergyDifference(system.forceField, system.simulationBox,
                                                                system.spanOfMoleculeAtoms(), {}, molecule) +
        Interactions::computeFrameworkMoleculeTailEnergyDifference(system.forceField, system.simulationBox,
                                                                   system.spanOfFrameworkAtoms(), {}, molecule);
    time_end = std::chrono::system_clock::now();
    // Update the CPU time statistics for the tail corrections
    component.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);
    system.mc_moves_cputime[move]["Tail"] += (time_end - time_begin);

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

    // Calculate the correction factor for Ewald summation
    double correctionFactorEwald =
        std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy()));

    // Compute acceptance probability factors
    double fugacity = component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
    double preFactor = correctionFactorEwald * double(system.numberOfIntegerMoleculesPerComponent[selectedComponent]) /
                       (system.beta * component.molFraction * fugacity * system.simulationBox.volume);
    double Pacc = preFactor;  // * idealGasRosenbluthWeight / retraceData.RosenbluthWeight;
    double biasTransitionMatrix = system.tmmc.biasFactor(oldN - 1, oldN);

    // Check if the new macrostate is within the allowed TMMC range
    if (system.tmmc.doTMMC)
    {
      std::size_t newN = oldN - 1;
      if (newN < system.tmmc.minMacrostate)
      {
        return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
      }
    }

    system.mc_moves_cputime[move]["Integration"] += (time_end - time_begin);
    component.mc_moves_statistics.addConstructed(move, 1);

    double drift = std::abs(currentEnergy.conservedEnergy() - referenceEnergy.conservedEnergy());

    if (random.uniform() < biasTransitionMatrix * Pacc * std::exp(-system.beta * drift))
    {
      component.mc_moves_statistics.addAccepted(move, 1);

      system.deleteMolecule(selectedComponent, selectedMolecule, molecule);
      system.moleculeData = moleculeData;
      system.thermostat = thermostat;
      system.timeStep = dt;

      std::copy(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), atomData.begin());
      system.spanOfMoleculeAtoms() = moleculeAtomPositions;
      Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components);
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      return {currentEnergy, double3(Pacc, 1.0 - Pacc, 0.0)};
    }

    return {std::nullopt, double3(Pacc, 1.0 - Pacc, 0.0)};
  }
}
