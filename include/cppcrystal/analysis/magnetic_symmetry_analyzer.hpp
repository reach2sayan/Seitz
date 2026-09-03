#pragma once

#include <cppcrystal/analysis/detail/lazy.hpp>
#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/magnetic_cell.hpp>
#include <cppcrystal/core/magnetic_symmetry_operation.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/data/msg_database.hpp>
#include <cppcrystal/magnetic_dataset.hpp>

#include <optional>
#include <type_traits>
#include <vector>

namespace cppcrystal::analysis {

// The magnetic counterpart of SymmetryAnalyzer: a persistent view over a
// MagneticCell + tolerances that lazily computes and memoizes the magnetic
// space-group analysis via get_magnetic_dataset. Same design as
// SymmetryAnalyzer — immutable after construction, value-owned inputs, the
// success value cached (detail::Lazy). The per-instance caches are not
// thread-safe; call warm() once on a single thread to make the const instance
// shareable read-only across threads.
class MagneticSymmetryAnalyzer {
public:
  // The rank-1 tensors transform as axial or polar vectors according to
  // `cell.kind()`; `tol.moment` (unset -> symprec) is the moment tolerance.
  [[nodiscard]] static MagneticSymmetryAnalyzer
  from_cell(MagneticCell cell, MagneticTolerance tol = {});

  [[nodiscard]] MagneticCell const &cell() const noexcept { return cell_; }
  [[nodiscard]] double symprec() const noexcept { return tol_.symprec; }
  [[nodiscard]] AngleTolerance angle_tolerance() const noexcept {
    return tol_.angle_tolerance;
  }

  [[nodiscard]] Result<MagneticDataset> dataset() const {
    BOOST_LEAF_AUTO(ds, cached_dataset());
    return *ds;
  }

  // Projections of the dataset.
  [[nodiscard]] Result<MagneticOperations> operations() const {
    return project<&MagneticDataset::operations>();
  }
  [[nodiscard]] Result<int> uni_number() const {
    return project<&MagneticDataset::uni_number>();
  }
  [[nodiscard]] Result<int> hall_number() const {
    return project<&MagneticDataset::hall_number>();
  }
  [[nodiscard]] Result<std::vector<int>> equivalent_atoms() const {
    return project<&MagneticDataset::equivalent_atoms>();
  }
  [[nodiscard]] Result<data::MagneticSpacegroupType> spacegroup_type() const {
    BOOST_LEAF_AUTO(uni, uni_number());
    return data::magnetic_spacegroup_type(uni);
  }

  // Standardized magnetic cell assembled from the dataset's std_* fields
  // (idealized lattice, positions, types, and the rotated site tensors).
  [[nodiscard]] Result<MagneticCell> standardized_cell() const;

  // Force the lazy dataset cache. After warm() succeeds this const instance
  // may be shared read-only across threads.
  Result<void> warm() const;

private:
  MagneticSymmetryAnalyzer(MagneticCell cell, MagneticTolerance tol)
      : cell_(std::move(cell)), tol_(tol) {}

  [[nodiscard]] Result<MagneticDataset const *> cached_dataset() const;

  // Copy one field out of the memoized dataset.
  template <auto Member> [[nodiscard]] auto project() const
      -> Result<std::remove_cvref_t<
          decltype(std::declval<MagneticDataset const &>().*Member)>> {
    BOOST_LEAF_AUTO(ds, cached_dataset());
    return ds->*Member;
  }

  MagneticCell cell_;
  MagneticTolerance tol_;

  detail::Lazy<MagneticDataset> dataset_;
};

} // namespace cppcrystal::analysis
