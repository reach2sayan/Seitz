#include <seitz/analysis/symmetry_analyzer.hpp>

#include "core/family_dispatch.hpp"
#include "core/validation.hpp"
#include "refine/refinement.hpp"
#include "spacegroup/spacegroup.hpp"
#include "symmetry/primitive.hpp"
#include "symmetry/search.hpp"

#include <algorithm>
#include <iterator>
#include <ranges>
#include <utility>
#include <vector>

namespace seitz::analysis {

namespace {

constexpr int kNumAttemptOuter = 10;
constexpr double kReduceRateOuter = 0.9; // tolerance shrink factor per attempt

// One determination attempt at one tolerance, for one family. std::nullopt
// means "this tolerance did not yield a consistent cell"; the caller tightens
// and retries.
template <GroupFamily F>
[[nodiscard]] std::optional<Dataset>
attempt(Cell const &cell, Tolerance const &tol,
        std::optional<HallNumber> setting) {
  symmetry::PrimitiveFinder<F> const finder(cell, tol);
  auto const primitive = finder.find();
  if (!primitive) {
    return std::nullopt;
  }
  spacegroup::SpacegroupMatcher<F> const matcher(*primitive, setting);
  auto const matched = matcher.search();
  if (!matched) {
    return std::nullopt;
  }

  // Orient the bravais lattice, recover the exact operations in the input
  // cell, and build the idealized cell + Wyckoff data.
  refine::Refinement<F> const refinement =
      refine::Refinement<F>{*matched, primitive->cell, cell,
                            primitive->tolerance}
          .similar_bravais();
  // Not const: moved into the Dataset below, after its last read.
  auto operations = refinement.operations();
  if (!operations) {
    return std::nullopt;
  }
  auto standardized =
      refinement.standardize(*operations, primitive->mapping_table);
  if (!standardized) {
    return std::nullopt;
  }

  SpacegroupMatch const &sg = refinement.matched();
  Lattice const bravais{sg.bravais_lattice};
  Lattice const std_lattice = standardized->bravais.lattice();

  // Zipped, not indexed: these five per-atom tables are parallel, and a shared
  // subscript is the one way assembling a Site can silently go wrong.
  std::vector<Site> sites;
  sites.reserve(standardized->wyckoffs.size());
  std::ranges::transform(
      std::views::zip(
          standardized->wyckoffs, standardized->site_symmetry_symbols,
          standardized->equivalent_atoms, standardized->crystallographic_orbits,
          primitive->mapping_table),
      std::back_inserter(sites), [](auto const &fields) {
        auto const &[wyckoff, symbol, equivalent, orbit, primitive_atom] =
            fields;
        return Site{wyckoff, symbol, equivalent, orbit, primitive_atom};
      });

  return Dataset{
      .hall = sg.hall,
      .bravais = bravais,
      .setting = {.transformation =
                      sg.bravais_lattice.inverse() * cell.lattice().matrix(),
                  .origin_shift = sg.origin_shift,
                  .rigid_rotation = bravais.rigid_rotation_to(std_lattice)},
      .operations = std::move(*operations),
      .sites = std::move(sites),
      .standardized = std::move(standardized->bravais),
      .std_mapping_to_primitive =
          std::move(standardized->std_mapping_to_primitive),
      .primitive = primitive->cell.lattice(),
  };
}

} // namespace

Result<Dataset> SymmetryAnalyzer::determine() const {
  if (auto valid = validate_cell(cell_); !valid) {
    return valid.error();
  }
  // The one runtime branch on the family; everything below is templated. The
  // outer loop progressively tightens the tolerance in case a given one yields
  // an inconsistent cell.
  return dispatch_family(cell_.periodicity(),
                         [&]<GroupFamily F>() -> Result<Dataset> {
                           Tolerance tol = tol_;
                           for (int i = 0; i < kNumAttemptOuter;
                                ++i, tol.symprec *= kReduceRateOuter) {
                             if (auto ds = attempt<F>(cell_, tol, setting_)) {
                               return std::move(*ds);
                             }
                           }
                           return leaf::new_error(e_spacegroup_search_failed{});
                         });
}

template <CellSetting S, Idealize I>
Result<Cell> SymmetryAnalyzer::standardized_cell() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return dispatch_family(
      cell_.periodicity(), [&]<GroupFamily F>() -> Result<Cell> {
        SpacegroupMatch const match{ds->hall, ds->bravais.matrix(),
                                    ds->setting.origin_shift};
        refine::Refinement<F> const refinement{match, cell_, cell_, tol_};

        if constexpr (I == Idealize::yes) {
          // The dataset already holds the standardized conventional cell.
          if constexpr (S == CellSetting::conventional) {
            return ds->standardized;
          } else {
            return refinement.to_primitive(ds->standardized,
                                           Matrix3d::Identity());
          }
        } else {
          // Transform the input cell instead, preserving its real geometry.
          BOOST_LEAF_AUTO(primitive, refinement.to_primitive(
                                         cell_, ds->setting.transformation));
          if constexpr (S == CellSetting::primitive) {
            return primitive;
          } else {
            return (match.type().centering == data::Centering::primitive)
                       ? primitive
                       : refinement.from_primitive(primitive);
          }
        }
      });
}

// clang-format off
template Result<Cell> SymmetryAnalyzer::standardized_cell<CellSetting::conventional, Idealize::yes>() const;
template Result<Cell> SymmetryAnalyzer::standardized_cell<CellSetting::conventional, Idealize::no>() const;
template Result<Cell> SymmetryAnalyzer::standardized_cell<CellSetting::primitive, Idealize::yes>() const;
template Result<Cell> SymmetryAnalyzer::standardized_cell<CellSetting::primitive, Idealize::no>() const;
// clang-format on


Result<Operations const &> SymmetryAnalyzer::cell_operations() const & {
  BOOST_LEAF_AUTO(ops, cell_operations_.get([&] {
    return dispatch_family(cell_.periodicity(), [&]<GroupFamily F>() {
      return symmetry::SymmetrySearch<F>{cell_, tol_}.operations();
    });
  }));
  return *ops;
}

Result<PointSymmetry const &> SymmetryAnalyzer::lattice_symmetry() const & {
  BOOST_LEAF_AUTO(ps, lattice_symmetry_.get([&] {
    return dispatch_family(cell_.periodicity(), [&]<GroupFamily F>() {
      return symmetry::SymmetrySearch<F>{cell_, tol_}.lattice_symmetry();
    });
  }));
  return *ps;
}

Result<Cell const &> SymmetryAnalyzer::primitive_cell() const & {
  BOOST_LEAF_AUTO(cell, primitive_cell_.get([&]() -> Result<Cell> {
    return dispatch_family(
        cell_.periodicity(), [&]<GroupFamily F>() -> Result<Cell> {
          symmetry::PrimitiveFinder<F> const finder(cell_, tol_);
          BOOST_LEAF_AUTO(prim, finder.find());
          return prim.cell;
        });
  }));
  return *cell;
}

Result<kpoint::ReciprocalMesh>
SymmetryAnalyzer::reciprocal_mesh(kpoint::Mesh mesh,
                                  TimeReversal time_reversal) const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return kpoint::ReciprocalMesh::from_rotations(
      mesh, ds->operations.rotations(), time_reversal);
}

} // namespace seitz::analysis
