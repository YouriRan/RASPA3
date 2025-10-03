module;

#ifdef USE_LEGACY_HEADERS
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <iostream>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#endif

#include "featomic.h"
#include "featomic.hpp"
#include "metatensor.h"
#include "metatensor.hpp"
export module featomic;

#ifndef USE_LEGACY_HEADERS
import std;
#endif

import json;
import atom;
import simulationbox;

export struct FeatomicSystem : featomic::System
{
  FeatomicSystem() {};

  FeatomicSystem(std::span<Atom> atoms, SimulationBox box) : positions_(3 * atoms.size()), atomic_types_(atoms.size())
  {
    for (std::size_t i = 0; i < atoms.size(); i++)
    {
      positions_[3 * i + 0] = atoms[i].position.x;
      positions_[3 * i + 1] = atoms[i].position.y;
      positions_[3 * i + 2] = atoms[i].position.z;
      atomic_types_[i] = static_cast<uint32_t>(atoms[i].type);
    }

    cell_ = box.cell.toArray();
  };

  uintptr_t size() const override { return this->atomic_types_.size(); }

  const int32_t* types() const override { return this->atomic_types_.data(); }

  const double* positions() const override { return this->positions_.data(); }

  CellMatrix cell() const override { return cell_; }

  void compute_neighbors(double /*cutoff*/) override
  {
    throw std::runtime_error("SimpleSystem can only be used with `use_native_system=true`");
  }

  const std::vector<featomic_pair_t>& pairs() const override
  {
    throw std::runtime_error("SimpleSystem can only be used with `use_native_system=true`");
  }

  const std::vector<featomic_pair_t>& pairs_containing(uintptr_t /*atom*/) const override
  {
    throw std::runtime_error("SimpleSystem can only be used with `use_native_system=true`");
  }

  std::vector<double> positions_;
  std::vector<int32_t> atomic_types_;
  CellMatrix cell_;
};

export struct FeatomicCalculator
{
  /// Create a new calculator with the given `name` and `parameters`.
  ///
  /// @throws FeatomicError if `name` is not associated with a known calculator,
  ///         if `parameters` is not valid JSON, or if `parameters` do not
  ///         contains the expected values for the requested calculator.
  ///
  /// @verbatim embed:rst:leading-slashes
  /// The list of available calculators and the corresponding parameters are
  /// in the :ref:`main documentation <userdoc-references>`. The ``parameters``
  /// should be formatted as JSON, according to the requested calculator
  /// schema.
  /// @endverbatim
  FeatomicCalculator(std::string name, std::string parameters)
      : calculator_(featomic_calculator(name.data(), parameters.data()))
  {
    if (this->calculator_ == nullptr)
    {
      throw featomic::FeatomicError(featomic_last_error());
    }
  }

  ~FeatomicCalculator() { featomic_calculator_free(this->calculator_); }

  /// NOTE: this is the same implementation as given in featomic, but implements a copy
  /// constructor. The copy constructor is needed when initializing a simulation with a list of
  /// systems. Beware that the copy constructor is defined as just creating a new calculator
  /// object with exactly the same parameters, losing the internal state of the calculator.
  FeatomicCalculator(const FeatomicCalculator& other)
  {
    auto nm = other.name();
    auto p = other.parameters();

    auto h = featomic_calculator(nm.data(), p.data());
    if (h == nullptr)
    {
      throw featomic::FeatomicError(featomic_last_error());
    }
    this->calculator_ = h;
  }
  FeatomicCalculator& operator=(const FeatomicCalculator& other)
  {
    if (this == &other)
    {
      return *this;
    }
    // First free whatever we already had:
    if (this->calculator_ != nullptr)
    {
      featomic_calculator_free(this->calculator_);
      this->calculator_ = nullptr;
    }
    // Then build a new one from the same name/params:
    auto nm = other.name();
    auto p = other.parameters();
    auto h = featomic_calculator(nm.data(), p.data());
    if (h == nullptr)
    {
      throw featomic::FeatomicError(featomic_last_error());
    }
    this->calculator_ = h;
    return *this;
  }

  /// Calculator is move-constructible
  FeatomicCalculator(FeatomicCalculator&& other) noexcept { *this = std::move(other); }

  /// Calculator can be move-assigned
  FeatomicCalculator& operator=(FeatomicCalculator&& other) noexcept
  {
    this->~FeatomicCalculator();
    this->calculator_ = nullptr;

    std::swap(this->calculator_, other.calculator_);

    return *this;
  }

  /// Get the name used to create this `Calculator`
  std::string name() const
  {
    auto buffer = std::vector<char>(32, '\0');
    while (true)
    {
      auto status = featomic_calculator_name(calculator_, buffer.data(), buffer.size());

      if (status != FEATOMIC_BUFFER_SIZE_ERROR)
      {
        featomic::details::check_status(status);
        return std::string(buffer.data());
      }

      // grow the buffer and retry
      buffer.resize(buffer.size() * 2, '\0');
    }
  }

  /// Get the parameters used to create this `Calculator`
  std::string parameters() const
  {
    auto buffer = std::vector<char>(256, '\0');
    while (true)
    {
      auto status = featomic_calculator_parameters(calculator_, buffer.data(), buffer.size());

      if (status != FEATOMIC_BUFFER_SIZE_ERROR)
      {
        featomic::details::check_status(status);
        return std::string(buffer.data());
      }

      // grow the buffer and retry
      buffer.resize(buffer.size() * 2, '\0');
    }
  }

  /// Get all radial cutoffs used by this `Calculator`'s neighbors lists
  std::vector<double> cutoffs() const
  {
    const double* data = nullptr;
    uintptr_t length = 0;
    featomic::details::check_status(featomic_calculator_cutoffs(calculator_, &data, &length));
    return std::vector<double>(data, data + length);
  }

  /// Runs a calculation with this calculator on the given ``systems``
  metatensor::TensorMap compute(std::vector<featomic_system_t>& systems,
                                featomic::CalculationOptions options = featomic::CalculationOptions()) const
  {
    mts_tensormap_t* descriptor = nullptr;

    featomic::details::check_status(featomic_calculator_compute(
        calculator_, &descriptor, systems.data(), systems.size(), options.as_featomic_calculation_options_t()));

    return metatensor::TensorMap(descriptor);
  }

  /// Runs a calculation for multiple `systems`
  template <typename SystemImpl,
            typename std::enable_if<std::is_base_of<FeatomicSystem, SystemImpl>::value, bool>::type = true>
  metatensor::TensorMap compute(std::vector<SystemImpl>& systems,
                                featomic::CalculationOptions options = featomic::CalculationOptions()) const
  {
    auto featomic_systems = std::vector<featomic_system_t>();
    for (auto& system : systems)
    {
      featomic_systems.push_back(system.as_featomic_system_t());
    }

    return this->compute(featomic_systems, std::move(options));
  }

  /// Runs a calculation for a single `system`
  template <typename SystemImpl,
            typename std::enable_if<std::is_base_of<FeatomicSystem, SystemImpl>::value, bool>::type = true>
  metatensor::TensorMap compute(SystemImpl& system,
                                featomic::CalculationOptions options = featomic::CalculationOptions()) const
  {
    mts_tensormap_t* descriptor = nullptr;

    auto featomic_system = system.as_featomic_system_t();
    featomic::details::check_status(featomic_calculator_compute(calculator_, &descriptor, &featomic_system, 1,
                                                                options.as_featomic_calculation_options_t()));

    return metatensor::TensorMap(descriptor);
  }

  /// Get the underlying pointer to a `featomic_calculator_t`.
  ///
  /// This is an advanced function that most users don't need to call
  /// directly.
  featomic_calculator_t* as_featomic_calculator_t() { return calculator_; }

  /// Get the underlying const pointer to a `featomic_calculator_t`.
  ///
  /// This is an advanced function that most users don't need to call
  /// directly.
  const featomic_calculator_t* as_featomic_calculator_t() const { return calculator_; }

  featomic_calculator_t* calculator_ = nullptr;
};

export FeatomicCalculator getSoapPowerSpectrumCalculator(double cutOff, double smoothingWidth, double gaussianWidth,
                                                         std::size_t numberOfRadialBasisFunctions,
                                                         std::size_t numberOfAngularBasisFunctions)
{
  const std::string parameters =
      std::format(R"json(
{{
  "cutoff": {{
    "radius": {},
    "smoothing": {{
      "type": "ShiftedCosine",
      "width": {}
    }}
  }},
  "density": {{
    "type": "Gaussian",
    "width": {},
    "center_atom_weight": 0
  }},
  "basis": {{
    "type": "TensorProduct",
    "max_angular": {},
    "radial": {{
      "type": "Gto",
      "max_radial": {}
    }}
  }}
}}
)json",
                  cutOff, smoothingWidth, gaussianWidth, numberOfAngularBasisFunctions, numberOfRadialBasisFunctions);
  return FeatomicCalculator("soap_power_spectrum", parameters);
}