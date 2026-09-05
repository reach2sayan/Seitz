#include <cppcrystal/analysis/dataset.hpp>
#include <cppcrystal/analysis/symmetry_analyzer.hpp>
#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/point_group.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/spacegroup_match.hpp>

#include "casters.hpp" // to_str
#include "errors.hpp"  // unwrap, memo

#include <pybind11/eigen.h>
#include <pybind11/native_enum.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace cppcrystal::python {

namespace {

using analysis::CellSetting;
using analysis::Dataset;
using analysis::Idealize;
using analysis::Setting;
using analysis::Site;
using analysis::SymmetryAnalyzer;

// standardized_cell is a template over (CellSetting, Idealize) with four
// explicit instantiations. Python cannot pass template arguments, so the
// runtime pair is dispatched here -- once, in the one place that knows the
// mapping.
[[nodiscard]] Cell standardized_cell(SymmetryAnalyzer const &self,
                                     CellSetting setting, bool idealize) {
  return unwrap([&]() -> Result<Cell> {
    py::gil_scoped_release const unlocked;
    if (setting == CellSetting::primitive) {
      return idealize ? self.standardized_cell<CellSetting::primitive,
                                               Idealize::yes>()
                      : self.standardized_cell<CellSetting::primitive,
                                               Idealize::no>();
    }
    return idealize ? self.standardized_cell<CellSetting::conventional,
                                             Idealize::yes>()
                    : self.standardized_cell<CellSetting::conventional,
                                             Idealize::no>();
  });
}

// The per-atom answers as parallel arrays, with no per-site Python object built
// at all. The vectorized route: a 5000-atom cell is 5000 model constructions
// through .sites, and this is what the docstring points large cells at.
[[nodiscard]] py::dict site_arrays(SymmetryAnalyzer const &self) {
  std::vector<Site> sites = memo(&SymmetryAnalyzer::sites)(self);
  auto const count = static_cast<Index>(sites.size());

  py::array_t<int> wyckoff(count);
  py::array_t<int> equivalent_atom(count);
  py::array_t<int> orbit(count);
  py::array_t<int> primitive_atom(count);
  py::list site_symmetry;

  auto w = wyckoff.mutable_unchecked<1>();
  auto e = equivalent_atom.mutable_unchecked<1>();
  auto o = orbit.mutable_unchecked<1>();
  auto p = primitive_atom.mutable_unchecked<1>();
  for (Index i = 0; i < count; ++i) {
    Site const &site = sites[static_cast<std::size_t>(i)];
    w(i) = site.wyckoff;
    e(i) = site.equivalent_atom;
    o(i) = site.orbit;
    p(i) = site.primitive_atom;
    site_symmetry.append(to_str(site.site_symmetry));
  }

  py::dict out;
  out["wyckoff"] = wyckoff;
  out["site_symmetry"] = site_symmetry;
  out["equivalent_atom"] = equivalent_atom;
  out["orbit"] = orbit;
  out["primitive_atom"] = primitive_atom;
  return out;
}

} // namespace

void bind_analysis(py::module_ &m) {
  py::native_enum<CellSetting>(m, "CellSetting", "enum.IntEnum",
                               "Which setting a standardized cell is expressed "
                               "in.")
      .value("conventional", CellSetting::conventional)
      .value("primitive", CellSetting::primitive)
      .finalize();

  py::native_enum<LatticeSetting>(
      m, "LatticeSetting", "enum.IntEnum",
      "Whether a lattice handed to a spacegroup search is the conventional "
      "cell or already a primitive one.")
      .value("conventional", LatticeSetting::conventional)
      .value("primitive", LatticeSetting::primitive)
      .finalize();

  // Idealize is deliberately NOT an enum here: Python has no idiomatic
  // two-valued enum, and `idealize=False` reads better than `Idealize.no`.

  py::class_<Setting>(m, "Setting",
                      "How the input cell maps onto the standardized setting.")
      .def_readonly("transformation", &Setting::transformation)
      .def_readonly("origin_shift", &Setting::origin_shift)
      .def_readonly("rigid_rotation", &Setting::rigid_rotation)
      .def("__repr__", [](Setting const &) { return "Setting(...)"; });

  py::class_<Site>(m, "Site", "The per-atom result of a determination.")
      .def_readonly("wyckoff", &Site::wyckoff, "Wyckoff letter index, 0 = 'a'.")
      // A copy to str. site_symmetry is a string_view into a constexpr table
      // that outlives the interpreter, but a Python str has to own its bytes.
      .def_property_readonly(
          "site_symmetry",
          [](Site const &self) { return to_str(self.site_symmetry); },
          "Tabulated site-symmetry symbol.")
      .def_readonly("equivalent_atom", &Site::equivalent_atom)
      .def_readonly("orbit", &Site::orbit)
      .def_readonly("primitive_atom", &Site::primitive_atom)
      .def("__repr__", [](Site const &self) {
        return "Site(wyckoff=" + std::to_string(self.wyckoff) + ", '" +
               std::string(self.site_symmetry) + "')";
      });

  py::class_<Dataset>(
      m, "Dataset",
      "The result of a space-group determination. The group's "
      "identity is the Hall setting alone -- number and symbols "
      "are one lookup away through spacegroup_type().")
      .def_readonly("hall", &Dataset::hall)
      .def_readonly("bravais", &Dataset::bravais)
      .def_readonly("setting", &Dataset::setting)
      .def_readonly("operations", &Dataset::operations)
      .def_readonly("sites", &Dataset::sites)
      .def_readonly("standardized", &Dataset::standardized)
      .def_readonly("std_mapping_to_primitive",
                    &Dataset::std_mapping_to_primitive)
      .def_readonly("primitive", &Dataset::primitive)
      .def("__repr__", [](Dataset const &self) {
        return "Dataset(hall=" + std::to_string(self.hall.index()) + ", " +
               std::to_string(self.operations.size()) + " operations)";
      });

  // ---- the analyzer ------------------------------------------------------
  //
  // No py::init: the analyzer is non-copyable (detail::Lazy holds a mutex) and
  // reached through a named factory, exactly as the C++ intends. Every
  // projection goes through memo(), which resolves the deleted-rvalue overload
  // pair and copies out of the memo rather than handing Python a reference that
  // would silently pin an analyzer holding two Cells behind every HallNumber.
  py::class_<SymmetryAnalyzer>(
      m, "SymmetryAnalyzer",
      "A persistent view over a cell plus tolerances that lazily computes and "
      "memoizes its determination. The object you keep, not a call you repeat. "
      "Every const query is thread-safe.")
      .def_static(
          "from_cell",
          [](Cell cell, Tolerance tol, std::optional<HallNumber> setting) {
            return SymmetryAnalyzer::from_cell(std::move(cell), tol, setting);
          },
          py::arg("cell"), py::arg("tolerance") = Tolerance{},
          py::arg("setting") = std::nullopt,
          py::doc("An unset `setting` searches every Hall setting of the "
                  "cell's family; a set one fixes it."))
      // Lambdas rather than member pointers: both are declared noexcept on the
      // Analyzer base, and pybind11's method_adaptor has no noexcept overload,
      // so a member pointer would bind with the unregistered CRTP base as its
      // self type and fail at call time. (See the same note in group.cpp.)
      .def_property_readonly(
          "cell", [](SymmetryAnalyzer const &self) { return self.cell(); })
      .def_property_readonly(
          "tolerance",
          [](SymmetryAnalyzer const &self) { return self.tolerance(); })
      .def_property_readonly(
          "dataset", memo_as<SymmetryAnalyzer>(&SymmetryAnalyzer::dataset),
          py::doc("The full determination."))
      .def_property_readonly("hall", memo(&SymmetryAnalyzer::hall))
      .def_property_readonly("operations", memo(&SymmetryAnalyzer::operations),
                             py::doc("Space-group operations of the input cell "
                                     "as given."))
      .def_property_readonly("sites", memo(&SymmetryAnalyzer::sites))
      .def_property_readonly("spacegroup_type",
                             memo(&SymmetryAnalyzer::spacegroup_type))
      // Explicit template arguments, unlike its siblings: standardized_cell
      // also names the (CellSetting, Idealize) member template below, and an
      // overload set holding a function template is a non-deduced context
      // ([temp.deduct.call]/6), so Self and T have to be given rather than
      // deduced. Naming them makes memo's parameter type concrete, which is
      // what picks the `const &` accessor back out of the set.
      .def_property_readonly("standardized_cell",
                             memo<SymmetryAnalyzer, Cell>(
                                 &SymmetryAnalyzer::standardized_cell),
                             py::doc("The standardized conventional, idealized "
                                     "cell."))
      .def_property_readonly("cell_operations",
                             memo(&SymmetryAnalyzer::cell_operations),
                             py::doc("All operations of the input cell, "
                                     "including the centering translations of "
                                     "a non-primitive cell."))
      // PointSymmetry is boost::container::static_vector<Matrix3i, 48>, which
      // pybind11/stl.h knows nothing about -- so unlike its sibling
      // projections this one cannot go straight through memo(). Copied into a
      // list of (3, 3) int32 arrays, which is what a caller wants anyway.
      .def_property_readonly(
          "lattice_symmetry",
          [](SymmetryAnalyzer const &self) {
            PointSymmetry const rotations =
                memo(&SymmetryAnalyzer::lattice_symmetry)(self);
            py::list out;
            for (Matrix3i const &rotation : rotations) {
              out.append(py::cast(rotation));
            }
            return out;
          },
          py::doc("The lattice point group: the rotations, in the cell basis, "
                  "that map the reduced lattice metric onto itself."))
      .def_property_readonly("primitive_cell",
                             memo(&SymmetryAnalyzer::primitive_cell))
      .def("standardized_cell_in", &standardized_cell, py::arg("setting"),
           py::arg("idealize") = true,
           py::doc("The standardized cell in another setting: conventional or "
                   "primitive, idealized or keeping the input's own geometry."))
      .def("site_arrays", &site_arrays,
           py::doc("The per-atom answers as parallel arrays, with no per-site "
                   "object built. Prefer this to .sites for large cells."))
      // Non-copyable by construction. Saying so with a sentence beats the
      // RuntimeError pybind11 would otherwise raise from a missing copy ctor.
      .def("__copy__",
           [](py::object const &) -> py::object {
             throw py::type_error("SymmetryAnalyzer is not copyable: it owns a "
                                  "memo. Share the object instead -- every "
                                  "const query on it is thread-safe.");
           })
      .def(
          "__deepcopy__",
          [](py::object const &, py::dict const &) -> py::object {
            throw py::type_error("SymmetryAnalyzer is not copyable: it owns a "
                                 "memo. Share the object instead -- every "
                                 "const query on it is thread-safe.");
          },
          py::arg("memo"))
      .def("__repr__", [](SymmetryAnalyzer const &self) {
        return "SymmetryAnalyzer(" + std::to_string(self.cell().size()) +
               " atoms, symprec=" + std::to_string(self.tolerance().symprec) +
               ")";
      });

  py::class_<SpacegroupMatch>(m, "SpacegroupMatch",
                              "The matched space group of a cell or an "
                              "operation set.")
      .def_readonly("hall", &SpacegroupMatch::hall)
      .def_readonly("bravais_lattice", &SpacegroupMatch::bravais_lattice)
      .def_readonly("origin_shift", &SpacegroupMatch::origin_shift)
      .def_property_readonly("type", &SpacegroupMatch::type,
                             py::return_value_policy::reference);
}

} // namespace cppcrystal::python
