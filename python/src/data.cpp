#include <seitz/core/keys.hpp>
#include <seitz/data/element_data.hpp>
#include <seitz/data/spg_database.hpp>

#include "casters.hpp" // to_str

#include <pybind11/native_enum.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace seitz::python {

namespace {

using data::Centering;
using data::SpacegroupType;

} // namespace

void bind_data(py::module_ &m) {
  py::native_enum<Centering>(m, "Centering", "enum.IntEnum",
                             "Bravais centering.")
      .value("error", Centering::error)
      .value("primitive", Centering::primitive)
      .value("body", Centering::body)
      .value("face", Centering::face)
      .value("a_face", Centering::a_face)
      .value("b_face", Centering::b_face)
      .value("c_face", Centering::c_face)
      .value("base", Centering::base)
      .value("r_center", Centering::r_center)
      .finalize();

  // Every string field is a string_view into a constexpr table. They outlive
  // the interpreter, but a Python str has to own its bytes, so each property
  // copies -- that is correct, not an optimization left undone.
  py::class_<SpacegroupType>(m, "SpacegroupType",
                             "One Hall setting's metadata. The Hall number "
                             "itself is not a field: it is the key that "
                             "addresses the row.")
      .def_readonly("number", &SpacegroupType::number,
                    "International space-group or layer-group number.")
      .def_property_readonly(
          "schoenflies",
          [](SpacegroupType const &self) { return to_str(self.schoenflies); })
      .def_property_readonly(
          "hall_symbol",
          [](SpacegroupType const &self) { return to_str(self.hall_symbol); })
      .def_property_readonly(
          "international",
          [](SpacegroupType const &self) { return to_str(self.international); })
      .def_property_readonly("international_full",
                             [](SpacegroupType const &self) {
                               return to_str(self.international_full);
                             })
      .def_property_readonly("international_short",
                             [](SpacegroupType const &self) {
                               return to_str(self.international_short);
                             })
      .def_property_readonly(
          "choice",
          [](SpacegroupType const &self) { return to_str(self.choice); })
      .def_readonly("centering", &SpacegroupType::centering)
      .def_readonly("pointgroup_number", &SpacegroupType::pointgroup_number,
                    "1..32.")
      .def("__repr__", [](SpacegroupType const &self) {
        return "SpacegroupType(" + std::to_string(self.number) + ", '" +
               std::string(self.international_short) + "')";
      });

  // Total -- a HallNumber has no invalid state, so this cannot fail and needs
  // no Result. Bound through a lambda because the C++ name is an overload set
  // (a template and a non-template), which &data::spacegroup_type cannot name.
  m.def(
      "spacegroup_type",
      [](HallNumber hall) -> SpacegroupType const & {
        return data::spacegroup_type(hall);
      },
      py::arg("hall"), py::return_value_policy::reference,
      py::doc("The metadata of a Hall setting."));

  // Lambdas, not &data::f: each name is an overload set (the family template
  // and its runtime-family form). The spans become list[int].
  m.def(
      "halls_with_number",
      [](GroupFamily family, int number) {
        return std::vector<int>{std::from_range,
                                data::halls_with_number(family, number)};
      },
      py::arg("family"), py::arg("number"),
        py::doc("Every Hall setting index of an international number, "
                "ascending; empty if out of range."));
  m.def(
      "default_hall",
      [](GroupFamily family, int number) {
        return data::default_hall(family, number);
      },
      py::arg("family"), py::arg("number"),
        py::doc("The default (first) Hall setting of an international number, "
                "or None."));
  m.def(
      "default_halls_with_pointgroup",
      [](GroupFamily family, int pointgroup_number) {
        return std::vector<int>{std::from_range,
                                data::default_halls_with_pointgroup(
                                    family, pointgroup_number)};
      },
      py::arg("family"), py::arg("pointgroup_number"),
        py::doc("The default Hall setting index of every group with that point "
                "group, ascending by group number."));
  m.def("operations_from_database", &data::operations_from_database,
        py::arg("hall"), py::return_value_policy::copy,
        py::doc("The symmetry operations of a Hall setting."));

  m.attr("K_NUM_SPACEGROUPS") = data::kNumSpacegroups;
  m.attr("K_NUM_LAYER_GROUPS") = data::kNumLayerGroups;
  m.attr("K_NUM_POINTGROUPS") = data::kNumPointgroups;

  // ---- elements ----------------------------------------------------------
  //
  // All of these answer with None rather than a sentinel, which is the
  // library's own rule and maps straight onto Python.
  py::module_ elements = m.def_submodule(
      "elements", "Tabulated per-element data, keyed by atomic number.");
  elements.def("is_known_element", &data::is_known_element, py::arg("z"));
  elements.def("covalent_radius", &data::covalent_radius, py::arg("z"),
               py::doc("Single-bond covalent radius in angstrom, or None."));
  elements.def("atomic_volume", &data::atomic_volume, py::arg("z"),
               py::doc("Sphere volume from the covalent radius, or None."));
  elements.def("element_symbol", &data::element_symbol, py::arg("z"),
               py::doc("Chemical symbol, or None."));
  elements.def("atomic_number", &data::atomic_number, py::arg("symbol"),
      py::doc("Atomic number of a chemical symbol (case-sensitive), or None."));
  elements.attr("K_NUM_ELEMENTS") = data::kNumElements;
}

} // namespace seitz::python
