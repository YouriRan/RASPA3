module;

#ifdef USE_LEGACY_HEADERS
#include <cstddef>
#include <cmath>
#include <array>
#include <cstring>
#include <fstream>
#include <map>
#include <print>
#include <string>
#include <type_traits>
#include <vector>
#include <tuple>
#endif

export module bend_potential;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import stringutils;
import archive;
import randomnumbers;
import double3;
import units;

/**
 * \brief Maximum number of parameters allowed for bend potentials.
 *
 * Defines the maximum number of parameters that can be associated with a bend potential.
 */
export const size_t maximumNumberOfBendParameters{4};

/**
 * \brief Enumeration of different bend types.
 *
 * Specifies the type of bend potential to be used in simulations.
 */
export enum class BendType : size_t { Rigid = 0, Fixed = 1, Harmonic = 2,
                                      CoreShell = 3, Quartic = 4, CFF_Quartic = 5,
                                      HarmonicCosine = 6, Cosine = 7, Tafipolsky = 8, 
                                      MM3 = 9, MM3_inplane = 10};

/**
 * \brief Represents a bend potential between two particles.
 *
 * The BendPotential struct encapsulates the type of bend and associated parameters between two particles.
 * It includes versioning for serialization, bend type, identifiers of bended particles, and bend parameters.
 */
export struct BendPotential
{
  uint64_t versionNumber{1};  ///< Version number for serialization.

  std::array<size_t, 3> identifiers;                    ///< Identifiers of the two particles forming the bend.
  BendType type;                                             ///< The type of bend potential.
  std::array<double, maximumNumberOfBendParameters> parameters;  ///< Parameters associated with the bend potential.

  /**
   * \brief Default constructor for BendPotential.
   *
   * Initializes a BendPotential object with Undefined bend type and zeroed bend IDs.
   */
  BendPotential() : identifiers({0, 0, 0}), type(BendType::Harmonic) {}

  BendPotential(std::array<size_t, 3> identifiers, BendType type, std::vector<double> vector_parameters);

  /**
   * \brief Constructs a BendPotential with specified type and bend IDs.
   *
   * \param type The type of bend potential.
   * \param identifiers A pair of particle identifiers forming the bend.
   */
  BendPotential(std::array<size_t, 3> identifiers, const BendType type) : identifiers(identifiers), type(type) {}

  bool operator==(BendPotential const &) const = default;

  /**
   * \brief Generates a string representation of the bend potential.
   *
   * Provides a formatted string containing bend type, particle IDs, and parameters.
   *
   * \return A string describing the bend potential.
   */
  std::string print() const;

  /**
   * \brief Number of parameters required for each bend type.
   *
   * A static vector indicating the number of parameters needed for each bend type.
   */
  static inline std::array<size_t, 11> numberOfBendParameters{1, 2, 2, 4, 4, 2, 3, 3, 1, 2, 2};


  /**
   * \brief Mapping of bend type strings to BendType enums.
   *
   * A static map that associates bend type names with their corresponding BendType enumeration values.
   */
  static inline std::map<std::string, BendType, caseInsensitiveComparator> definitionForString{
      {"RIGID", BendType::Rigid}, {"FIXED", BendType::Fixed}, {"HARMONIC", BendType::Harmonic},
      {"CORE_SHELL", BendType::CoreShell}, {"QUARTIC", BendType::Quartic}, {"CFF_QUARTIC", BendType::CFF_Quartic},
      {"HARMONIC_COSINE", BendType::HarmonicCosine}, {"COSINE", BendType::Cosine},
      {"TAFIPOLSKY", BendType::Tafipolsky}, {"MM3", BendType::MM3}, {"MM3_INPLANE", BendType::MM3_inplane}};

  double generateBendAngle(RandomNumber &random, double beta) const;

  double calculateEnergy(const double3 &posA, const double3 &posB, const double3 &posc, std::optional<const double3> &posD) const;

  friend Archive<std::ofstream> &operator<<(Archive<std::ofstream> &archive, const BendPotential &b);
  friend Archive<std::ifstream> &operator>>(Archive<std::ifstream> &archive, BendPotential &b);
};
