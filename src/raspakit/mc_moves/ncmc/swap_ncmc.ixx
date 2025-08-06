module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <optional>
#include <span>
#include <tuple>
#endif

export module mc_moves_noneq_cbmc;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import double3;
import randomnumbers;
import running_energy;
import atom;
import system;

export namespace MC_Moves
{
/**
 * \brief Performs a volume move in a Monte Carlo simulation.
 *
 * Attempts to change the system volume by a random amount, scaling the simulation box and molecule positions
 * accordingly. Computes the new energies and decides whether to accept the move based on the Metropolis criterion.
 *
 * \param random A random number generator.
 * \param system The simulation system to modify.
 * \return The new total energy if the move is accepted; std::nullopt otherwise.
 */
std::pair<std::optional<RunningEnergy>, double3> NonEqCBMC(RandomNumber& random, System& system,
                                                           std::size_t selectedComponent);
}