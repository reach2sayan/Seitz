#pragma once

#include <seitz/analysis/analyzer.hpp>
#include <seitz/analysis/dataset.hpp>
#include <seitz/core/cell.hpp>
#include <seitz/core/error.hpp>
#include <seitz/core/keys.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/point_group.hpp>
#include <seitz/core/tolerance.hpp>
#include <seitz/data/spg_database.hpp>
#include <seitz/kpoint/mesh.hpp>

#include <optional>
#include <span>
#include <utility>

#pragma GCC visibility push(default)

namespace seitz::analysis {

// The setting of the standardized cell, and whether its lattice is
// metric-idealized or the input's own (possibly distorted) geometry. Template
// arguments: the four combinations resolve at compile time.
enum class CellSetting { conventional, primitive };
enum class Idealize { yes, no };

struct SpaceGroupTraits {
  using CellType = Cell;
  using DatasetType = Dataset;
  using ToleranceType = Tolerance;
};

class SymmetryAnalyzer : public Analyzer<SymmetryAnalyzer, SpaceGroupTraits> {
public:
  // An unset `setting` searches every Hall setting of the cell's family; a set
  // one fixes it.
  [[nodiscard]] static SymmetryAnalyzer
  from_cell(Cell cell, Tolerance tol = {},
            std::optional<HallNumber> setting = std::nullopt) {
    return SymmetryAnalyzer{std::move(cell), tol, setting};
  }

  // Projections of the dataset.
  [[nodiscard]] Result<HallNumber const &> hall() const & {
    return project<&Dataset::hall>();
  }
  Result<HallNumber const &> hall() const && = delete;
  [[nodiscard]] Result<Operations const &> operations() const & {
    return project<&Dataset::operations>();
  }
  Result<Operations const &> operations() const && = delete;
  [[nodiscard]] Result<std::vector<Site> const &> sites() const & {
    return project<&Dataset::sites>();
  }
  Result<std::vector<Site> const &> sites() const && = delete;
  [[nodiscard]] Result<data::SpacegroupType const &> spacegroup_type() const {
    BOOST_LEAF_AUTO(setting, hall());
    return data::spacegroup_type(setting);
  }

  // The standardized conventional idealized cell, as the dataset holds it.
  [[nodiscard]] Result<Cell const &> standardized_cell() const & {
    return project<&Dataset::standardized>();
  }
  Result<Cell const &> standardized_cell() const && = delete;

  // The standardized cell in another setting: primitive vs conventional,
  // idealized vs the input's geometry. Not memoized -- keyed by its template
  // arguments, and the conventional/idealized case is the accessor above.
  template <CellSetting S, Idealize I>
  [[nodiscard]] Result<Cell> standardized_cell() const;

  // References into their own memos, &-qualified like the dataset projections:
  // by value, lattice_symmetry() copies a static_vector<Matrix3i, 48> per
  // call.
  [[nodiscard]] Result<Operations const &> cell_operations() const &;
  Result<Operations const &> cell_operations() const && = delete;

  // The lattice point group: the rotations (in the cell basis) that map the
  // Delaunay-reduced lattice metric onto itself.
  [[nodiscard]] Result<PointSymmetry const &> lattice_symmetry() const &;
  Result<PointSymmetry const &> lattice_symmetry() const && = delete;

  // The primitive cell, cached independently of the full dataset
  [[nodiscard]] Result<Cell const &> primitive_cell() const &;
  Result<Cell const &> primitive_cell() const && = delete;

  // The irreducible reciprocal mesh: the determination's rotations made
  // reciprocal (with the inversion partner under time reversal), reducing
  // `mesh`.
  [[nodiscard]] Result<kpoint::ReciprocalMesh>
  reciprocal_mesh(kpoint::Mesh mesh, TimeReversal time_reversal) const;

private:
  friend Analyzer;
  SymmetryAnalyzer(Cell cell, Tolerance tol, std::optional<HallNumber> setting)
      : Analyzer(std::move(cell), tol), setting_(setting) {}

  // The pipeline itself: find primitive -> match Hall setting -> refine
  [[nodiscard]] Result<Dataset> determine() const;

  std::optional<HallNumber> setting_;

  detail::Lazy<Cell> primitive_cell_;
  detail::Lazy<Operations> cell_operations_;
  detail::Lazy<PointSymmetry> lattice_symmetry_;
};

extern template Result<Cell>
SymmetryAnalyzer::standardized_cell<CellSetting::conventional, Idealize::yes>()
    const;
extern template Result<Cell>
SymmetryAnalyzer::standardized_cell<CellSetting::conventional, Idealize::no>()
    const;
extern template Result<Cell>
SymmetryAnalyzer::standardized_cell<CellSetting::primitive, Idealize::yes>()
    const;
extern template Result<Cell>
SymmetryAnalyzer::standardized_cell<CellSetting::primitive, Idealize::no>()
    const;

} // namespace seitz::analysis

#pragma GCC visibility pop
