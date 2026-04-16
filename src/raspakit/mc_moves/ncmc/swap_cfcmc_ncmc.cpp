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

std::pair<std::optional<RunningEnergy>, double3> MC_Moves::swapNCMC(RandomNumber& random, System& system,
                                                                    std::size_t selectedComponent,
                                                                    std::size_t selectedMolecule,
                                                                    bool insertionDisabled, bool deletionDisabled)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::SwapNCMC;
  Component& component = system.components[selectedComponent];

  // Retrieve lambda parameters and select a new lambda bin for the move
  PropertyLambdaProbabilityHistogram& lambda = component.lambdaGC;
  std::size_t oldBin = lambda.currentBin;
  double deltaLambda = lambda.delta;
  double maxChange = component.mc_moves_statistics.getMaxChange(move, 2);
  std::make_signed_t<std::size_t> selectedNewBin = lambda.selectNewBin(random, maxChange);
  std::size_t oldN = system.numberOfIntegerMoleculesPerComponent[selectedComponent];

  std::size_t indexFractionalMolecule = system.indexOfGCFractionalMoleculesPerComponent_CFCMC(selectedComponent);

  // all copied data: moleculeData, moleculeAtomPositions, thermostat, dt
  // all const data: components, forcefield, simulationbox, numberofmoleculespercomponents, fixedFrameworkStoredEik
  // all scratch data: eik_x, eik_y, eik_z, eik_xy
  std::span<Atom> atomData = system.spanOfMoleculeAtoms();
  std::vector<Atom> moleculeAtomPositions(atomData.size());
  std::copy(atomData.begin(), atomData.end(), moleculeAtomPositions.begin());

  std::vector<Molecule> moleculeData(system.moleculeData);
  std::optional<Thermostat> thermostat(system.thermostat);

  std::vector<std::size_t> numberOfMoleculesPerComponent(system.numberOfMoleculesPerComponent);

  RunningEnergy referenceEnergy = system.runningEnergies;
  referenceEnergy.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(moleculeData);
  referenceEnergy.rotationalKineticEnergy =
      Integrators::computeRotationalKineticEnergy(moleculeData, system.components);
  RunningEnergy currentEnergy = referenceEnergy;

  // get Timestep from the max change
  double dt = system.timeStep;
  std::size_t numberOfPerturbations = lambda.numberOfSamplePoints - 1;

  bool insert = (random.uniform() < 0.5);

  if (insert)  // Insertion move
  {
    component.mc_moves_statistics.addTrial(move, 0);

    std::pair<Molecule, std::vector<Atom>> trialMolecule =
        component.equilibratedMoleculeRandomInBox(random, system.simulationBox);

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

    // insert molecule into copied atom positions
    std::size_t index = 0;
    for (std::size_t comp = 0; comp < selectedComponent + 1; ++comp)
    {
      index += system.components[comp].atoms.size() * system.numberOfMoleculesPerComponent[comp];
    }

    std::vector<Atom>::const_iterator iterator =
        moleculeAtomPositions.begin() + static_cast<std::vector<Atom>::difference_type>(index);
    iterator = moleculeAtomPositions.insert(iterator, trialMolecule.second.begin(), trialMolecule.second.end());
    std::span<Atom> inserted(&moleculeAtomPositions[index], trialMolecule.second.size());

    index = std::accumulate(system.numberOfMoleculesPerComponent.begin(),
                            system.numberOfMoleculesPerComponent.begin() + selectedComponent + 1, 0uz);

    std::vector<Molecule>::iterator moleculeIterator =
        moleculeData.begin() + static_cast<std::vector<Atom>::difference_type>(index);
    moleculeData.insert(moleculeIterator, trialMolecule.first);

    numberOfMoleculesPerComponent[selectedComponent] += 1;

    for (Atom& atom : inserted)
    {
      atom.setScaling(1.0 / numberOfPerturbations);
    }

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

    if (intermediateEnergy.potentialEnergy() > system.forceField.energyOverlapCriteria)
      return {std::nullopt, double3(0.0, 1.0, 0.0)};

    intermediateEnergy.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(moleculeData);
    intermediateEnergy.rotationalKineticEnergy =
        Integrators::computeRotationalKineticEnergy(moleculeData, system.components);
    RunningEnergy currentEnergy = intermediateEnergy;

    double work = currentEnergy.potentialEnergy() - system.runningEnergies.potentialEnergy();

    for (std::size_t step = 0; step < system.numberOfHybridMCSteps; ++step)
    {
      currentEnergy = Integrators::velocityVerlet(
          moleculeData, moleculeAtomPositions, system.components, dt, thermostat, system.spanOfFrameworkAtoms(),
          system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
          system.totalEik, system.fixedFrameworkStoredEik, system.interpolationGrids, numberOfMoleculesPerComponent);
    }

    double drift = currentEnergy.conservedEnergy() - intermediateEnergy.conservedEnergy();

    for (std::size_t i = 1; i < numberOfPerturbations; i++)
    {
      double newLambda = (i + 1) / numberOfPerturbations;
      for (Atom& atom : inserted)
      {
        atom.setScaling(newLambda);
      }

      RunningEnergy intermediateEnergy = Integrators::updateGradients(
          moleculeAtomPositions, system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox,
          system.components, system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik,
          system.fixedFrameworkStoredEik, system.interpolationGrids, numberOfMoleculesPerComponent);
      intermediateEnergy.translationalKineticEnergy = currentEnergy.translationalKineticEnergy;
      intermediateEnergy.rotationalKineticEnergy = currentEnergy.rotationalKineticEnergy;

      if (intermediateEnergy.potentialEnergy() > system.forceField.energyOverlapCriteria)
        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      work += intermediateEnergy.potentialEnergy() - currentEnergy.potentialEnergy();

      for (std::size_t step = 0; step < system.numberOfHybridMCSteps; ++step)
      {
        currentEnergy = Integrators::velocityVerlet(
            moleculeData, moleculeAtomPositions, system.components, dt, thermostat, system.spanOfFrameworkAtoms(),
            system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
            system.totalEik, system.fixedFrameworkStoredEik, system.interpolationGrids, numberOfMoleculesPerComponent);
      }
      drift += currentEnergy.conservedEnergy() - intermediateEnergy.conservedEnergy();
    }

    for (Atom& atom : inserted)
    {
      atom.setScalingToInteger();
    }

    // Calculate correction factor for Ewald energy difference
    double correctionFactorEwald = 1.0;
    // std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy()));

    // Compute the acceptance probability pre-factor
    double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
    double preFactor = correctionFactorEwald * system.beta * fugacity * system.simulationBox.volume /
                       double(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);
    double Pacc = preFactor * std::exp(-system.beta * (work + std::abs(drift)));

    component.mc_moves_statistics.addConstructed(move, 0);
    std::cout << "insertion: " << work << " " << drift << " " << Pacc << std::endl;

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
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }

  else  // Deletion move
  {
    component.mc_moves_statistics.addTrial(move, 1);
    if (oldN == 0)
    {
      return {std::nullopt, double3(0.0, 0.0, 0.0)};
    }

    // Get a reference to the molecule being deleted
    std::size_t selectedMolecule = system.randomIntegerMoleculeOfComponent(random, selectedComponent);

    std::size_t index = 0;
    for (std::size_t comp = 0; comp < selectedComponent; ++comp)
    {
      index += system.components[comp].atoms.size() * system.numberOfMoleculesPerComponent[comp];
    }
    index += selectedMolecule * system.components[selectedComponent].atoms.size();

    std::span<Atom> toBeDeleted(&moleculeAtomPositions[index], system.components[selectedComponent].atoms.size());

    for (Atom& atom : toBeDeleted)
    {
      atom.setScaling(1 - (1.0 / numberOfPerturbations));
      atom.groupId = true;
      atom.isFractional = true;
    }

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

    if (intermediateEnergy.potentialEnergy() > system.forceField.energyOverlapCriteria)
      return {std::nullopt, double3(0.0, 1.0, 0.0)};

    intermediateEnergy.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(moleculeData);
    intermediateEnergy.rotationalKineticEnergy =
        Integrators::computeRotationalKineticEnergy(moleculeData, system.components);
    RunningEnergy currentEnergy = intermediateEnergy;

    double work = currentEnergy.potentialEnergy() - system.runningEnergies.potentialEnergy();

    for (std::size_t step = 0; step < system.numberOfHybridMCSteps; ++step)
    {
      currentEnergy = Integrators::velocityVerlet(
          moleculeData, moleculeAtomPositions, system.components, dt, thermostat, system.spanOfFrameworkAtoms(),
          system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
          system.totalEik, system.fixedFrameworkStoredEik, system.interpolationGrids, numberOfMoleculesPerComponent);
    }

    double drift = currentEnergy.conservedEnergy() - intermediateEnergy.conservedEnergy();

    for (std::size_t i = 1; i < numberOfPerturbations; i++)
    {
      double newLambda = 1 - (i + 1) / numberOfPerturbations;
      for (Atom& atom : toBeDeleted)
      {
        atom.setScaling(newLambda);
      }

      RunningEnergy intermediateEnergy = Integrators::updateGradients(
          moleculeAtomPositions, system.spanOfFrameworkAtoms(), system.forceField, system.simulationBox,
          system.components, system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik,
          system.fixedFrameworkStoredEik, system.interpolationGrids, numberOfMoleculesPerComponent);
      intermediateEnergy.translationalKineticEnergy = currentEnergy.translationalKineticEnergy;
      intermediateEnergy.rotationalKineticEnergy = currentEnergy.rotationalKineticEnergy;

      if (intermediateEnergy.potentialEnergy() > system.forceField.energyOverlapCriteria)
        return {std::nullopt, double3(0.0, 1.0, 0.0)};
      work += intermediateEnergy.potentialEnergy() - currentEnergy.potentialEnergy();

      for (std::size_t step = 0; step < system.numberOfHybridMCSteps; ++step)
      {
        currentEnergy = Integrators::velocityVerlet(
            moleculeData, moleculeAtomPositions, system.components, dt, thermostat, system.spanOfFrameworkAtoms(),
            system.forceField, system.simulationBox, system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
            system.totalEik, system.fixedFrameworkStoredEik, system.interpolationGrids, numberOfMoleculesPerComponent);
      }
      drift += currentEnergy.conservedEnergy() - intermediateEnergy.conservedEnergy();
    }

    // Calculate correction factor for Ewald energy difference
    double correctionFactorEwald = 1.0;
    // std::exp(-system.beta * (energyFourierDifference.potentialEnergy() + tailEnergyDifference.potentialEnergy()));

    // Compute the acceptance probability pre-factor
    double fugacity = component.molFraction * component.fugacityCoefficient.value_or(1.0) * system.pressure;
    double idealGasRosenbluthWeight = component.idealGasRosenbluthWeight.value_or(1.0);
    double preFactor = correctionFactorEwald * system.beta * fugacity * system.simulationBox.volume /
                       double(1 + system.numberOfIntegerMoleculesPerComponent[selectedComponent]);
    double Pacc = preFactor * std::exp(-system.beta * (work + std::abs(drift)));

    component.mc_moves_statistics.addConstructed(move, 1);
    std::cout << "deletion: " << work << " " << drift << " " << Pacc << std::endl;

    if (random.uniform() < Pacc)
    {
      component.mc_moves_statistics.addAccepted(move, 1);

      std::copy(moleculeData.begin(), moleculeData.end(), system.moleculeData.begin());
      std::span<Atom> newAtomData = system.spanOfMoleculeAtoms();
      std::copy(moleculeAtomPositions.begin(), moleculeAtomPositions.end(), newAtomData.begin());

      std::span<Atom> molecule = system.spanOfMolecule(selectedComponent, selectedMolecule);
      system.deleteMolecule(selectedComponent, selectedMolecule, molecule);
      system.thermostat = thermostat;
      system.timeStep = dt;
      // system.spanOfMoleculeAtoms() = moleculeAtomPositions;

      system.updateMoleculeAtomInformation();
      Integrators::createCartesianPositions(system.moleculeData, system.spanOfMoleculeAtoms(), system.components);
      Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
      return {currentEnergy, double3(Pacc, 1.0 - Pacc, 0.0)};
    }
    return {std::nullopt, double3(0.0, 1.0, 0.0)};
  }
}