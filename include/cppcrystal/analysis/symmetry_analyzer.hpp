#pragma once

#include <cppcrystal/analysis/analyzer.hpp>
#include <cppcrystal/analysis/dataset.hpp>
#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/point_group.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/data/spg_database.hpp>
#include <cppcrystal/kpoint/mesh.hpp>

#include <optional>
#include <span>
#include <utility>

namespace cppcrystal::analysis {

// Which setting the standardized cell is expressed in, and whether it carries
// the metric-idealized (symmetrized) lattice or keeps the input's real —
// possibly distorted — geometry. Template arguments, so the four combinations
// are resolved at compile time.
enum class CellSetting { conventional, primitive };
enum class Idealize { yes, no };

struct SpaceGroupTraits {
  using CellType = Cell;
  using DatasetType = Dataset;
  using ToleranceType = Tolerance;
};

class SymmetryAnalyzer : public Analyzer<SymmetryAnalyzer, SpaceGroupTraits> {
public:
  // Named factory (no overloaded constructors). An unset `setting` searches
  // every Hall setting of the cell's family; a set one fixes it.
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

  // The standardized conventional, idealized cell — the one the dataset
  // already holds.
  [[nodiscard]] Result<Cell const &> standardized_cell() const & {
    return project<&Dataset::standardized>();
  }
  Result<Cell const &> standardized_cell() const && = delete;

  // The standardized cell in another setting: primitive vs conventional,
  // idealized vs the input's own geometry. Not memoized — keyed by its
  // template arguments, and the conventional/idealized case is the accessor
  // above.
  template <CellSetting S, Idealize I>
  [[nodiscard]] Result<Cell> standardized_cell() const;

  // All space-group operations of the input cell exactly as given, including
  // the centering translations of a non-primitive cell. Distinct from
  // operations(), which are the dataset's operations in the input basis.
  //
  // A reference into the memo, and &-qualified for it, like the dataset
  // projections above: these three used to copy their cache on every call, and
  // for lattice_symmetry() that copy is a static_vector<Matrix3i, 48>.
  [[nodiscard]] Result<Operations const &> cell_operations() const &;
  Result<Operations const &> cell_operations() const && = delete;

  // The lattice point group: the rotations (in the cell basis) that map the
  // Delaunay-reduced lattice metric onto itself.
  [[nodiscard]] Result<PointSymmetry const &> lattice_symmetry() const &;
  Result<PointSymmetry const &> lattice_symmetry() const && = delete;

  // The primitive cell, cached independently of the full dataset so a caller
  // that only wants it does not pay for standardization.
  [[nodiscard]] Result<Cell const &> primitive_cell() const &;
  Result<Cell const &> primitive_cell() const && = delete;

  // The irreducible reciprocal mesh of this crystal: the determination's
  // rotations, made reciprocal (adding the inversion partner with time
  // reversal), reducing `mesh`. The Builder that replaces
  // ReciprocalMeshBuilder — the analyzer already owns the cell and tolerance.
  [[nodiscard]] Result<kpoint::ReciprocalMesh>
  reciprocal_mesh(kpoint::Mesh mesh, TimeReversal time_reversal) const;

private:
  friend Analyzer;
  SymmetryAnalyzer(Cell cell, Tolerance tol, std::optional<HallNumber> setting)
      : Analyzer(std::move(cell), tol), setting_(setting) {}

  // The pipeline itself: find primitive -> match Hall setting -> refine, at
  // progressively tighter tolerances until one attempt yields a consistent
  // cell. One runtime branch on the family at the top; everything below it is
  // compile-time specialised.
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

} // namespace cppcrystal::analysis
