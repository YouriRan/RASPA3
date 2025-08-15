module;

#ifdef USE_LEGACY_HEADERS
#include <optional>
#include <span>
#include <vector>
#endif

export module cbmc_rigid_deletion;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import atom;
import double3x3;
import double3;
import randomnumbers;
import forcefield;
import simulationbox;
import cbmc_chain_data;
import cbmc_interactions;
import framework;
import component;
import interpolation_energy_grid;


export namespace CBMC
{
[[nodiscard]] ChainRetraceData retraceRigidMoleculeChainDeletion(RandomNumber &random, const Component &component, bool hasExternalField,
                                          const ForceField &forcefield, const SimulationBox &simulationBox,
                                          const std::vector<std::optional<InterpolationEnergyGrid>> &interpolationGrids,
                                          const std::optional<Framework> &framework,
                                          std::span<const Atom> frameworkAtomData, std::span<const Atom> moleculeAtomData,
                                          double beta, double cutOffFrameworkVDW, double cutOffMoleculeVDW,
                                          double cutOffCoulomb, std::size_t startingBead, std::span<Atom> molecule_atoms,
                                          std::size_t numberOfTrialDirections) noexcept;
}
