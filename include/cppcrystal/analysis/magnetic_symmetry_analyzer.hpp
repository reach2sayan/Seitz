#pragma once

#include <cppcrystal/analysis/detail/lazy.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/magnetic_cell.hpp>
#include <cppcrystal/core/magnetic_symmetry_operation.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/data/msg_database.hpp>
#include <cppcrystal/magnetic_dataset.hpp>
#include <cppcrystal/spin/spin.hpp>

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
  // `is_axial` selects axial-vector transformation of the rank-1 tensors;
  // `mag_symprec` (std::nullopt -> symprec) is the moment tolerance.
  [[nodiscard]] static MagneticSymmetryAnalyzer
  from_cell(MagneticCell cell, bool is_axial,
            double symprec = kDefaultSymprec,
            AngleTolerance angle_tolerance = std::nullopt,
            std::optional<double> mag_symprec = std::nullopt);

  [[nodiscard]] MagneticCell const &cell() const noexcept { return cell_; }
  [[nodiscard]] bool is_axial() const noexcept { return is_axial_; }
  [[nodiscard]] double symprec() const noexcept { return tol_.symprec; }
  [[nodiscard]] AngleTolerance angle_tolerance() const noexcept {
    return tol_.angle_tolerance;
  }

  [[nodiscard]] Result<MagneticDataset> dataset() const {
    BOOST_LEAF_AUTO(ds, cached_dataset());
    return *ds;
  }

  // Projections of the dataset.
  [[nodiscard]] Result<MagneticSymmetryOperations> operations() const {
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

  // The full magnetic symmetry search of the input cell: the magnetic
  // operations, each atom's magnetic-orbit representative, the atom
  // permutations under every operation, and the primitive lattice implied by
  // the magnetic pure translations. Derives the spatial symmetry
  // (symmetry::find_symmetry of the underlying cell) and runs
  // spin::operations_with_site_tensors with time reversal. Distinct from
  // operations(), which is just the dataset's magnetic operations. Cached
  // independently.
  [[nodiscard]] Result<spin::MagneticSymmetrySearch> symmetry_search() const;

  // Standardized magnetic cell assembled from the dataset's std_* fields
  // (idealized lattice, positions, types, and the rotated site tensors).
  [[nodiscard]] Result<MagneticCell> standardized_cell() const;

  // Force both lazy caches (the magnetic dataset and the independently-cached
  // symmetry search). After warm() succeeds this const instance may be shared
  // read-only across threads.
  Result<void> warm() const;

private:
  MagneticSymmetryAnalyzer(MagneticCell cell, bool is_axial, Tolerance tol,
                           std::optional<double> mag_symprec)
      : cell_(std::move(cell)), is_axial_(is_axial), tol_(tol),
        mag_symprec_(mag_symprec) {}

  [[nodiscard]] Result<MagneticDataset const *> cached_dataset() const;

  // Copy one field out of the memoized dataset.
  template <auto Member> [[nodiscard]] auto project() const
      -> Result<std::remove_cvref_t<
          decltype(std::declval<MagneticDataset const &>().*Member)>> {
    BOOST_LEAF_AUTO(ds, cached_dataset());
    return ds->*Member;
  }

  MagneticCell cell_;
  bool is_axial_ = false;
  Tolerance tol_;
  std::optional<double> mag_symprec_;

  detail::Lazy<MagneticDataset> dataset_;
  detail::Lazy<spin::MagneticSymmetrySearch> symmetry_search_;
};

} // namespace cppcrystal::analysis
