#include <seitz/core/operation_set.hpp>
#include <seitz/core/point_group.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/core/types.hpp>

#include "casters.hpp" // to_str
#include "errors.hpp"  // detail::raise, unwrap

#include <pybind11/eigen.h>
#include <pybind11/native_enum.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstddef>
#include <string>
#include <string_view>

namespace seitz::python {

namespace {

[[nodiscard]] std::string operation_repr(SymmetryOperation const &self) {
  return "SymmetryOperation(det=" +
         std::to_string(self.rotation.determinant()) + ", identity_rotation=" +
         (self.is_identity_rotation() ? "True" : "False") + ")";
}

} // namespace

void bind_core_symmetry(py::module_ &m) {
  py::native_enum<TimeReversal>(
      m, "TimeReversal", "enum.IntEnum",
      "Whether an operation set includes its time-reversal partners.")
      .value("off", TimeReversal::off)
      .value("on", TimeReversal::on)
      .finalize();

  py::native_enum<Holohedry>(m, "Holohedry", "enum.IntEnum",
                             "Crystal system / holohedry.")
      .value("none", Holohedry::none)
      .value("triclinic", Holohedry::triclinic)
      .value("monoclinic", Holohedry::monoclinic)
      .value("orthorhombic", Holohedry::orthorhombic)
      .value("tetragonal", Holohedry::tetragonal)
      .value("trigonal", Holohedry::trigonal)
      .value("hexagonal", Holohedry::hexagonal)
      .value("cubic", Holohedry::cubic)
      .finalize();

  py::native_enum<Laue>(m, "Laue", "enum.IntEnum", "Laue class.")
      .value("none", Laue::none)
      .value("laue_1", Laue::laue_1)
      .value("laue_2m", Laue::laue_2m)
      .value("laue_mmm", Laue::laue_mmm)
      .value("laue_4m", Laue::laue_4m)
      .value("laue_4mmm", Laue::laue_4mmm)
      .value("laue_3", Laue::laue_3)
      .value("laue_3m", Laue::laue_3m)
      .value("laue_6m", Laue::laue_6m)
      .value("laue_6mmm", Laue::laue_6mmm)
      .value("laue_m3", Laue::laue_m3)
      .value("laue_m3m", Laue::laue_m3m)
      .finalize();

  // The 32 crystal classes. Enumerator values are the international point-group
  // numbering, so CrystalClass(n) round-trips with pointgroup_by_number(n).
  py::native_enum<CrystalClass>(m, "CrystalClass", "enum.IntEnum",
                                "The 32 crystallographic point groups, named "
                                "by Schoenflies symbol.")
      .value("none", CrystalClass::none)
      .value("c1", CrystalClass::c1)
      .value("ci", CrystalClass::ci)
      .value("c2", CrystalClass::c2)
      .value("cs", CrystalClass::cs)
      .value("c2h", CrystalClass::c2h)
      .value("d2", CrystalClass::d2)
      .value("c2v", CrystalClass::c2v)
      .value("d2h", CrystalClass::d2h)
      .value("c4", CrystalClass::c4)
      .value("s4", CrystalClass::s4)
      .value("c4h", CrystalClass::c4h)
      .value("d4", CrystalClass::d4)
      .value("c4v", CrystalClass::c4v)
      .value("d2d", CrystalClass::d2d)
      .value("d4h", CrystalClass::d4h)
      .value("c3", CrystalClass::c3)
      .value("c3i", CrystalClass::c3i)
      .value("d3", CrystalClass::d3)
      .value("c3v", CrystalClass::c3v)
      .value("d3d", CrystalClass::d3d)
      .value("c6", CrystalClass::c6)
      .value("c3h", CrystalClass::c3h)
      .value("c6h", CrystalClass::c6h)
      .value("d6", CrystalClass::d6)
      .value("c6v", CrystalClass::c6v)
      .value("d3h", CrystalClass::d3h)
      .value("d6h", CrystalClass::d6h)
      .value("t", CrystalClass::t)
      .value("th", CrystalClass::th)
      .value("o", CrystalClass::o)
      .value("td", CrystalClass::td)
      .value("oh", CrystalClass::oh)
      .finalize();

  // ---- symmetry operations -----------------------------------------------
  //
  // No __eq__. tolerance.hpp is emphatic that this library has exactly one
  // definition of "equal within tolerance", and an exact == on floating-point
  // translations would quietly compete with it. same_operation() below is the
  // comparison, and it takes the symprec that makes the answer meaningful.
  py::class_<SymmetryOperation>(
      m, "SymmetryOperation",
      "x -> rotation . x + translation, on fractional coordinates.")
      .def(py::init([](Matrix3i const &rotation, Vector3d const &translation) {
             return SymmetryOperation{rotation, translation};
           }),
           py::arg("rotation") = Matrix3i(Matrix3i::Identity()),
           py::arg("translation") = Vector3d(Vector3d::Zero()))
      .def_readonly("rotation", &SymmetryOperation::rotation)
      .def_readonly("translation", &SymmetryOperation::translation)
      .def("apply", &SymmetryOperation::apply, py::arg("x"))
      .def(py::self * py::self)
      .def_property_readonly(
          "inverse",
          [](SymmetryOperation const &self) { return self.inverse(); },
          py::doc("The inverse, or None when the rotation is not unimodular."))
      .def_property_readonly("is_identity_rotation",
                             &SymmetryOperation::is_identity_rotation)
      .def(py::pickle(
          [](SymmetryOperation const &self) {
            return py::make_tuple(self.rotation, self.translation);
          },
          [](py::tuple const &t) {
            return SymmetryOperation{t[0].cast<Matrix3i>(),
                                     t[1].cast<Vector3d>()};
          }))
      .def("__copy__", [](SymmetryOperation const &self) { return self; })
      .def(
          "__deepcopy__",
          [](SymmetryOperation const &self, py::dict const &) { return self; },
          py::arg("memo"))
      .def("to_xyz", &to_xyz,
           py::doc("The Jones-faithful coordinate triplet ('x,y,z', "
                   "'-y,x-y,z+1/2') CIF's symop loop carries."))
      .def_static(
          "from_xyz",
          [](std::string_view text) {
            return unwrap([&] { return from_xyz(text); });
          },
          py::arg("text"),
          py::doc("The inverse reading, permissive about spacing, case and "
                  "whether translations are fractions or decimals. Raises "
                  "InvalidXyzError on text that is not three coordinates."))
      .def("__repr__", &operation_repr);

  m.def("same_operation", &same_operation, py::arg("a"), py::arg("b"),
        py::arg("symprec"),
        py::doc("Whether two operations agree to within `symprec`. This is the "
                "comparison to use -- SymmetryOperation has no __eq__, because "
                "an exact float comparison would compete with the library's "
                "one definition of equal-within-tolerance."));

  m.def(
      "conjugated_by",
      [](SymmetryOperation const &op, Matrix3d const &t,
         Matrix3d const &t_inv) { return conjugated_by(op, t, t_inv); },
      py::arg("op"), py::arg("t"), py::arg("t_inv"),
      py::doc("Change of basis of one operation: (T,0)(R,t)(T,0)^-1."));

  // ---- operation sets ----------------------------------------------------
  //
  // The one range in this module where py::make_iterator is genuinely safe:
  // ops_ is an owned std::vector, so keep_alive<0, 1> on the set is enough.
  // Contrast Cell.atoms and Mesh.addresses, which are lazy views over a
  // temporary and have to be materialised.
  py::class_<Operations>(m, "Operations",
                         "An immutable set of space-group operations.")
      .def(py::init<std::vector<SymmetryOperation>>(), py::arg("operations"))
      .def("__len__", &Operations::size)
      .def(
          "__iter__",
          [](Operations const &self) {
            return py::make_iterator(self.begin(), self.end());
          },
          py::keep_alive<0, 1>())
      .def("__getitem__",
           [](Operations const &self, Index i) {
             auto const count = static_cast<Index>(self.size());
             if (i < 0) {
               i += count;
             }
             if (i < 0 || i >= count) {
               throw py::index_error("operation index out of range");
             }
             return self[static_cast<std::size_t>(i)];
           })
      .def_property_readonly("empty", &Operations::empty)
      .def_property_readonly("rotations", &Operations::rotations,
                             py::doc("The rotation parts, in order."))
      .def_property_readonly("pure_translations",
                             &Operations::pure_translations,
                             py::doc("The identity-rotation translations, "
                                     "including zero."))
      .def("conjugated_by", &Operations::conjugated_by, py::arg("t"),
           py::arg("t_inv"))
      .def("__repr__", [](Operations const &self) {
        return "Operations(" + std::to_string(self.size()) + ")";
      });

  // ---- point-group metadata ----------------------------------------------
  //
  // seitz::PointGroup, the plain metadata row -- distinct from
  // seitz::group::PointGroup, the 0D group class that arrives in phase 2.
  // Named PointGroupType here so the two never collide in one namespace.
  py::class_<PointGroup>(m, "PointGroupType",
                         "Metadata for one of the 32 crystallographic point "
                         "groups. Not the 0D group object.")
      .def_readonly("number", &PointGroup::number)
      .def_property_readonly(
          "symbol", [](PointGroup const &self) { return to_str(self.symbol); })
      .def_property_readonly(
          "schoenflies",
          [](PointGroup const &self) { return to_str(self.schoenflies); })
      .def_readonly("holohedry", &PointGroup::holohedry)
      .def_readonly("laue", &PointGroup::laue)
      .def_readonly("crystal_class", &PointGroup::crystal_class)
      .def("__repr__", [](PointGroup const &self) {
        return "PointGroupType(" + std::to_string(self.number) + ", '" +
               std::string(self.symbol) + "')";
      });

  m.def("pointgroup_by_number", &pointgroup_by_number, py::arg("number"),
        py::doc("Metadata for point group 1..32; number 0 or out of range "
                "gives an empty row."));
}

} // namespace seitz::python
