#pragma once

#include <cppcrystal/analysis/analyzer.hpp>
#include <cppcrystal/analysis/dataset.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/magnetic_cell.hpp>
#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/data/msg_database.hpp>

#include <utility>
#include <vector>

namespace cppcrystal::analysis {

struct MagneticTraits {
  using CellType = MagneticCell;
  using DatasetType = MagneticDataset;
  using ToleranceType = MagneticTolerance;
};

// The magnetic counterpart of SymmetryAnalyzer: same Facade, same cache, same
// projections; only determine() differs. The rank-1 tensors transform as axial
// or polar vectors according to cell.kind(); tol.moment (unset -> symprec) is
// the moment tolerance.
class MagneticSymmetryAnalyzer
    : public Analyzer<MagneticSymmetryAnalyzer, MagneticTraits> {
public:
  [[nodiscard]] static MagneticSymmetryAnalyzer
  from_cell(MagneticCell cell, MagneticTolerance tol = {}) {
    return MagneticSymmetryAnalyzer{std::move(cell), tol};
  }

  [[nodiscard]] Result<UniNumber const &> uni() const & {
    return project<&MagneticDataset::uni>();
  }
  Result<UniNumber const &> uni() const && = delete;
  [[nodiscard]] Result<HallNumber const &> hall() const & {
    return project<&MagneticDataset::hall>();
  }
  Result<HallNumber const &> hall() const && = delete;
  [[nodiscard]] Result<MagneticOperations const &> operations() const & {
    return project<&MagneticDataset::operations>();
  }
  Result<MagneticOperations const &> operations() const && = delete;
  [[nodiscard]] Result<std::vector<int> const &> equivalent_atoms() const & {
    return project<&MagneticDataset::equivalent_atoms>();
  }
  Result<std::vector<int> const &> equivalent_atoms() const && = delete;
  [[nodiscard]] Result<MagneticCell const &> standardized_cell() const & {
    return project<&MagneticDataset::standardized>();
  }
  Result<MagneticCell const &> standardized_cell() const && = delete;
  [[nodiscard]] Result<data::MagneticSpacegroupType> spacegroup_type() const {
    BOOST_LEAF_AUTO(number, uni());
    return data::magnetic_spacegroup_type(number);
  }

private:
  friend Analyzer;
  MagneticSymmetryAnalyzer(MagneticCell cell, MagneticTolerance tol)
      : Analyzer(std::move(cell), tol) {}

  // Spatial symmetry -> magnetic symmetry -> UNI identification -> idealize
  // -> transform into the standardized setting.
  [[nodiscard]] Result<MagneticDataset> determine() const;
};

} // namespace cppcrystal::analysis
