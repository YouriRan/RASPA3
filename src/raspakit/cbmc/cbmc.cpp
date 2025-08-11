module;

#ifdef USE_LEGACY_HEADERS
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <optional>
#include <span>
#include <tuple>
#include <vector>
#endif

module cbmc;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import randomnumbers;
import component;
import atom;
import molecule;
import double3;
import double3x3;
import simulationbox;
import energy_status;
import forcefield;
import energy_factor;
import cbmc_rigid_insertion;
import cbmc_rigid_deletion;
import cbmc_rigid_reinsertion;
import cbmc_flexible_insertion;
import cbmc_flexible_deletion;
//import cbmc_flexible_reinsertion;
import framework;
import component;
import cbmc_chain_data;
import interpolation_energy_grid;

[[nodiscard]] std::optional<ChainData> CBMC::growMoleculeSwapInsertion(
    RandomNumber &random, Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, Component::GrowType growType, double cutOffFrameworkVDW,
    double cutOffMoleculeVDW, double cutOffCoulomb, std::size_t selectedComponent, std::size_t selectedMolecule,
    double scaling, bool groupId, bool isFractional, std::size_t numberOfTrialDirections) noexcept
{
  switch (growType)
  {
    case Component::GrowType::Rigid:
      return CBMC::growRigidMoleculeSwapInsertion(
          random, component, hasExternalField, components, forceField, simulationBox, interpolationGrids, framework,
          frameworkAtoms, moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent,
          selectedMolecule, scaling, groupId, isFractional, numberOfTrialDirections);
    case Component::GrowType::Flexible:
      return CBMC::growFlexibleMoleculeSwapInsertion(
          random, component, hasExternalField, components, forceField, simulationBox, interpolationGrids, framework,
          frameworkAtoms, moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent,
          selectedMolecule, scaling, groupId, isFractional, numberOfTrialDirections);
    default:
      std::unreachable();
  }
}

[[nodiscard]] std::optional<ChainData> CBMC::growMoleculeReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, Component::GrowType growType, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::size_t selectedComponent, std::size_t selectedMolecule, Molecule &molecule,
    std::span<Atom> molecule_atoms, std::size_t numberOfTrialDirections) noexcept
{
  switch (growType)
  {
    case Component::GrowType::Rigid:
      return CBMC::growRigidMoleculeReinsertion(random, component, hasExternalField, components, forceField, simulationBox,
                                                interpolationGrids, framework, frameworkAtoms, moleculeAtoms, beta,
                                                cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent,
                                                selectedMolecule, molecule, molecule_atoms, numberOfTrialDirections);
    case Component::GrowType::Flexible:
      // TODO!
      return CBMC::growRigidMoleculeReinsertion(random, component, hasExternalField, components, forceField, simulationBox,
                                                interpolationGrids, framework, frameworkAtoms, moleculeAtoms, beta,
                                                cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent,
                                                selectedMolecule, molecule, molecule_atoms, numberOfTrialDirections);
    default:
      std::unreachable();
  }
}

[[nodiscard]] ChainData CBMC::retraceMoleculeReinsertion(
    RandomNumber &random, const Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, Component::GrowType growType, 
    double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::size_t selectedComponent, std::size_t selectedMolecule, Molecule &molecule,
    std::span<Atom> molecule_atoms, double storedR, std::size_t numberOfTrialDirections) noexcept
{
  switch (growType)
  {
    case Component::GrowType::Rigid:
      return CBMC::retraceRigidMoleculeReinsertion(
          random, component, hasExternalField, components, forceField, simulationBox, interpolationGrids, framework,
          frameworkAtoms, moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent,
          selectedMolecule, molecule, molecule_atoms, storedR, numberOfTrialDirections);
    case Component::GrowType::Flexible:
      // TODO!
      return CBMC::retraceRigidMoleculeReinsertion(
          random, component, hasExternalField, components, forceField, simulationBox, interpolationGrids, framework,
          frameworkAtoms, moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent,
          selectedMolecule, molecule, molecule_atoms, storedR, numberOfTrialDirections);
    default:
      std::unreachable();
  }
}

[[nodiscard]] ChainData CBMC::retraceMoleculeSwapDeletion(
    RandomNumber &random, const Component &component, bool hasExternalField, const std::vector<Component> &components,
    const ForceField &forceField, const SimulationBox &simulationBox,
    const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
    const std::optional<Framework> &framework, std::span<const Atom> frameworkAtoms,
    std::span<const Atom> moleculeAtoms, double beta, Component::GrowType growType, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
    double cutOffCoulomb, std::size_t selectedComponent, std::size_t selectedMolecule, std::span<Atom> molecule,
    double scaling, std::size_t numberOfTrialDirections) noexcept
{
  switch (growType)
  {
    case Component::GrowType::Rigid:
      return CBMC::retraceRigidMoleculeSwapDeletion(
          random, component, hasExternalField, components, forceField, simulationBox, interpolationGrids, framework,
          frameworkAtoms, moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent,
          selectedMolecule, molecule, scaling, numberOfTrialDirections);
    case Component::GrowType::Flexible:
      return CBMC::retraceFlexibleMoleculeSwapDeletion(
          random, component, hasExternalField, components, forceField, simulationBox, interpolationGrids, framework,
          frameworkAtoms, moleculeAtoms, beta, cutOffFrameworkVDW, cutOffMoleculeVDW, cutOffCoulomb, selectedComponent,
          selectedMolecule, molecule, scaling, numberOfTrialDirections);
    default:
      std::unreachable();
  }

}
