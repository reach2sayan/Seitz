#include <seitz/core/lattice.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/core/types.hpp>
#include <seitz/kpoint/mesh.hpp>

#include "casters.hpp" // array_1d, array_n3

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace seitz::python {

namespace {

using kpoint::Address;
using kpoint::BrillouinZone;
using kpoint::Mesh;
using kpoint::ReciprocalMesh;

// Mesh's checked door, as Lattice has: `Mesh(...)` raises, `Mesh.of(...)`
// answers None.
[[nodiscard]] Mesh mesh_from(Address divisions, std::array<bool, 3> shift) {
  if (auto const mesh = Mesh::of(divisions, shift)) {
    return *mesh;
  }
  throw py::value_error("mesh divisions must be strictly positive");
}

[[nodiscard]] std::string address_text(Address const &a) {
  return "(" + std::to_string(a[0]) + ", " + std::to_string(a[1]) + ", " +
         std::to_string(a[2]) + ")";
}

} // namespace

void bind_kpoint(py::module_ &m) {
  py::class_<Mesh>(m, "Mesh",
                   "Reciprocal-space grid geometry: the map between integer "
                   "addresses and linear grid-point indices. Double-mesh "
                   "convention: q = (2 * address + shift) / (2 * divisions).")
      .def(py::init(&mesh_from), py::arg("divisions"),
           py::arg("shift") = std::array<bool, 3>{},
           py::doc("Raises ValueError unless every division is positive."))
      .def_static(
          "of",
          [](Address divisions, std::array<bool, 3> shift) {
            return Mesh::of(divisions, shift);
          },
          py::arg("divisions"), py::arg("shift") = std::array<bool, 3>{},
          py::doc("The mesh, or None when a division is not positive."))
      .def_property_readonly("divisions", &Mesh::divisions)
      .def_property_readonly("shift", &Mesh::shift)
      .def("__len__", &Mesh::size)
      .def("index_of", &Mesh::index_of, py::arg("address"),
           py::doc("Linear index of an address, folded into the grid first."))
      .def("index_of_doubled", &Mesh::index_of_doubled, py::arg("doubled"))
      .def("address_of", &Mesh::address_of, py::arg("index"),
           py::doc("The address at a linear index, folded onto the "
                   "parallelepiped."))
      // addresses() is a lazy view over `this`, so it is materialised: one
      // (N, 3) int32 array, in grid-point-index order.
      .def_property_readonly(
          "addresses",
          [](Mesh const &self) {
            std::vector<Address> const rows(std::from_range, self.addresses());
            return array_n3<int>(rows);
          },
          py::doc("Every address as an (N, 3) int32 array, in index order."))
      .def("doubled_address", &Mesh::doubled_address, py::arg("address"))
      .def("doubled", &Mesh::doubled,
           py::doc("The mesh at twice the resolution."))
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def("__repr__", [](Mesh const &self) {
        return "Mesh(" + address_text(self.divisions()) + ")";
      });

  py::class_<ReciprocalMesh>(
      m, "ReciprocalMesh",
      "The symmetry reduction of a mesh: the reciprocal point group and, per "
      "grid point, the index of its irreducible representative.")
      .def_static(
          "from_rotations",
          [](Mesh mesh, std::vector<Matrix3i> const &rotations,
             TimeReversal time_reversal) {
            py::gil_scoped_release const unlocked;
            return ReciprocalMesh::from_rotations(mesh, rotations,
                                                  time_reversal);
          },
          py::arg("mesh"), py::arg("rotations"),
          py::arg("time_reversal") = TimeReversal::on,
          py::doc("Reduce `mesh` by real-space rotations (transposed to "
                  "reciprocal space, plus the inversion partner under time "
                  "reversal)."))
      .def(
          "stabilized",
          [](ReciprocalMesh const &self, std::vector<Vector3d> const &qpoints) {
            py::gil_scoped_release const unlocked;
            return self.stabilized(qpoints);
          },
          py::arg("qpoints"),
          py::doc("Reduced only by the rotations that map the q-point set onto "
                  "itself."))
      .def_property_readonly("mesh", &ReciprocalMesh::mesh)
      .def_property_readonly(
          "rotations",
          [](ReciprocalMesh const &self) {
            return std::vector<Matrix3i>(std::from_range, self.rotations());
          },
          py::doc("The reciprocal point group."))
      .def_property_readonly(
          "mapping",
          [](ReciprocalMesh const &self) {
            return array_1d<std::size_t>(self.mapping());
          },
          py::doc("mapping[i] == i exactly when i is an irreducible "
                  "representative; an (N,) uint64 array."))
      .def_property_readonly("num_irreducible", &ReciprocalMesh::num_irreducible)
      .def("images_of", &ReciprocalMesh::images_of, py::arg("address"),
           py::doc("The grid-point index one address maps to under each "
                   "rotation."))
      .def(
          "brillouin_zone",
          [](ReciprocalMesh const &self, Lattice const &reciprocal) {
            py::gil_scoped_release const unlocked;
            return self.brillouin_zone(reciprocal);
          },
          py::arg("reciprocal"),
          py::doc("Relocate the grid into the first Brillouin zone of "
                  "`reciprocal` (columns = reciprocal basis vectors)."))
      .def("__repr__", [](ReciprocalMesh const &self) {
        return "ReciprocalMesh(" + address_text(self.mesh().divisions()) +
               ", " + std::to_string(self.num_irreducible()) + " irreducible)";
      });

  py::class_<BrillouinZone>(
      m, "BrillouinZone",
      "A grid in the first Brillouin zone: each point at its image nearest "
      "the origin, boundary points duplicated. The map is on the doubled mesh.")
      .def_property_readonly(
          "addresses",
          [](BrillouinZone const &self) {
            return array_n3<int>(self.addresses());
          },
          py::doc("The in-BZ points as an (M, 3) int32 array: the first "
                  "len(mesh) at their original index, then the duplicates."))
      .def("map", &BrillouinZone::map, py::arg("doubled_index"),
           py::doc("Doubled-mesh index -> index into addresses, or None."))
      .def("images_of", &BrillouinZone::images_of, py::arg("address"))
      .def("__repr__", [](BrillouinZone const &self) {
        return "BrillouinZone(" + std::to_string(self.addresses().size()) +
               " points)";
      });
}

} // namespace seitz::python
