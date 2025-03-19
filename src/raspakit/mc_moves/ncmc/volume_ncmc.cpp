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
#include <numeric>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

module mc_moves_volume_ncmc;

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
import <numeric>;
import <chrono>;
import <cmath>;
import <iostream>;
import <iomanip>;
#endif

import component;
import atom;
import molecule;
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
import running_energy;
import property_lambda_probability_histogram;
import property_widom;
import averages;
import forcefield;
import interactions_framework_molecule;
import interactions_intermolecular;
import interactions_ewald;
import integrators;
import integrators_update;
import integrators_compute;
import thermostat;
import units;
import mc_moves_move_types;

std::optional<RunningEnergy> MC_Moves::volumeMoveNCMC(RandomNumber &random, System &system)
{
  std::chrono::system_clock::time_point time_begin, time_end;
  MoveTypes move = MoveTypes::VolumeNCMC;

  // Update volume move counts
  system.mc_moves_statistics.addTrial(move);

  // Calculate the total number of molecules
  double numberOfMolecules = static_cast<double>(std::reduce(system.numberOfIntegerMoleculesPerComponent.begin(),
                                                             system.numberOfIntegerMoleculesPerComponent.end()));
  double oldVolume = system.simulationBox.volume;
  // double maxVolumeChange = system.mc_moves_statistics.getMaxChange(move);
  double maxVolumeChange = 0.1;

  // Propose a new volume change
  double newVolume = std::exp(std::log(oldVolume) + maxVolumeChange * (2.0 * random.uniform() - 1.0));
  // Compute scaling factor for box dimensions
  double scale = std::pow(newVolume / oldVolume, 1.0 / 3.0);

  SimulationBox newBox = system.simulationBox.scaled(scale);
  std::pair<std::vector<Molecule>, std::vector<Atom>> newPositions = system.scaledCenterOfMassPositions(scale);

  std::vector<Molecule> moleculePositions = newPositions.first;
  std::vector<Atom> moleculeAtomPositions = newPositions.second;
  std::optional<Thermostat> thermostat(system.thermostat);

  double dt = system.timeStep;

  // initialize the velocities according to Boltzmann distribution
  // NOTE: it is important that the reference energy has the initial kinetic energies
  Integrators::initializeVelocities(random, moleculePositions, system.components, system.temperature);

  if (system.numberOfFrameworkAtoms == 0)
  {
    Integrators::removeCenterOfMassVelocityDrift(moleculePositions);
  }

  RunningEnergy referenceEnergy = system.runningEnergies;
  referenceEnergy.translationalKineticEnergy = Integrators::computeTranslationalKineticEnergy(moleculePositions);
  referenceEnergy.rotationalKineticEnergy =
      Integrators::computeRotationalKineticEnergy(moleculePositions, system.components);
  RunningEnergy currentEnergy = referenceEnergy;

  currentEnergy = Integrators::updateGradients(moleculeAtomPositions, system.spanOfFrameworkAtoms(), system.forceField, newBox, system.components,
                               system.eik_x, system.eik_y, system.eik_z, system.eik_xy, system.totalEik,
                               system.fixedFrameworkStoredEik, system.numberOfMoleculesPerComponent);

  // integrate for N steps
  time_begin = std::chrono::system_clock::now();
  for (size_t step = 0; step < system.numberOfHybridMCSteps; ++step)
  {
    currentEnergy = Integrators::velocityVerlet(
        moleculePositions, moleculeAtomPositions, system.components, dt, thermostat, system.spanOfFrameworkAtoms(),
        system.forceField, newBox, system.eik_x, system.eik_y, system.eik_z, system.eik_xy,
        system.totalEik, system.fixedFrameworkStoredEik, system.numberOfMoleculesPerComponent);
  }
  time_end = std::chrono::system_clock::now();

  system.mc_moves_cputime[move]["Integration"] += (time_end - time_begin);
  system.mc_moves_statistics.addConstructed(move);

  double drift = currentEnergy.conservedEnergy() - referenceEnergy.conservedEnergy();

  // Update constructed move counts
  system.mc_moves_statistics.addConstructed(move);

  // Apply acceptance/rejection rule
  if (random.uniform() < std::exp((numberOfMolecules + 1.0) * std::log(newVolume / oldVolume) -
                                  system.beta * (system.pressure * (newVolume - oldVolume) + drift)))
  {
    system.mc_moves_statistics.addAccepted(move);

    system.moleculePositions = moleculePositions;
    system.thermostat = thermostat;
    system.timeStep = dt;
    system.simulationBox = newBox;

    system.spanOfMoleculeAtoms() = moleculeAtomPositions;

    Integrators::createCartesianPositions(system.moleculePositions, system.spanOfMoleculeAtoms(), system.components);
    Interactions::acceptEwaldMove(system.forceField, system.storedEik, system.totalEik);
    return currentEnergy;
  }

  return std::nullopt;
}
