#include <seitz/core/cell.hpp>
#include <seitz/core/keys.hpp>
#include <seitz/core/magnetic_cell.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/types.hpp>
#include <seitz/data/msg_database.hpp>

#include "casters.hpp" // array_1d, to_str
#include "errors.hpp"  // detail::raise

#include <pybind11/eigen.h>
#include <pybind11/native_enum.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <optional>
#include <string>
#include <utility>
#include <variant>

namespace seitz::python {

namespace {

// The rank is the array's shape: (N,) is collinear, (N, 3) is non-collinear.
// The variant alternative IS the rank on the C++ side, so nothing here can
// fall out of step with it.
[[nodiscard]] SiteTensors tensors_from(py::array_t<double> const &moments,
                                       Index atoms) {
  if (moments.ndim() == 1 && moments.shape(0) == atoms) {
    auto const view = moments.unchecked<1>();
    CollinearTensors out(static_cast<std::size_t>(atoms));
    for (Index i = 0; i < atoms; ++i) {
      out[static_cast<std::size_t>(i)] = view(i);
    }
    return out;
  }
  if (moments.ndim() == 2 && moments.shape(0) == atoms && moments.shape(1) == 3) {
    auto const view = moments.unchecked<2>();
    NoncollinearTensors out(atoms, 3);
    for (Index i = 0; i < atoms; ++i) {
      for (Index k = 0; k < 3; ++k) {
        out(i, k) = view(i, k);
      }
    }
    return out;
  }
  throw py::value_error("moments must be shaped (N,) or (N, 3) for N atoms");
}

// The tensors back as the array they came from.
[[nodiscard]] py::object moments_of(MagneticCell const &self) {
  return std::visit(
      [](auto const &tensors) -> py::object {
        if constexpr (std::same_as<std::decay_t<decltype(tensors)>,
                                   CollinearTensors>) {
          return array_1d<double>(tensors);
        } else {
          return py::cast(tensors);
        }
      },
      self.tensors());
}

} // namespace

void bind_magnetic(py::module_ &m) {
  py::native_enum<TensorKind>(m, "TensorKind", "enum.IntEnum",
                              "How a rank-1 site tensor transforms under an "
                              "improper operation.")
      .value("polar", TensorKind::polar)
      .value("axial", TensorKind::axial)
      .finalize();

  py::native_enum<SiteTensor>(m, "SiteTensor", "enum.IntEnum",
                              "Rank of the per-site tensors.")
      .value("none", SiteTensor::none)
      .value("collinear", SiteTensor::collinear)
      .value("noncollinear", SiteTensor::noncollinear)
      .finalize();

  py::class_<MagneticCell>(m, "MagneticCell",
                           "A Cell plus per-site magnetic moments: (N,) scalars "
                           "for a collinear structure, (N, 3) vectors for a "
                           "non-collinear one.")
      .def(py::init([](Cell cell, py::array_t<double> const &moments,
                       TensorKind kind) {
             SiteTensors tensors = tensors_from(moments, cell.size());
             return MagneticCell(std::move(cell), std::move(tensors), kind);
           }),
           py::arg("cell"), py::arg("moments"),
           py::arg("kind") = TensorKind::polar)
      .def_property_readonly("cell", &MagneticCell::cell,
                             py::return_value_policy::copy)
      .def_property_readonly("moments", &moments_of,
                             py::doc("The moments as given: (N,) or (N, 3)."))
      .def_property_readonly("kind", &MagneticCell::kind)
      .def_property_readonly("rank", &MagneticCell::rank)
      .def("__len__", &MagneticCell::size)
      .def("__repr__", [](MagneticCell const &self) {
        return "MagneticCell(" + std::to_string(self.size()) + " atoms, " +
               (self.rank() == SiteTensor::collinear ? "collinear"
                                                     : "noncollinear") +
               ")";
      });

  // ---- the magnetic database -----------------------------------------------
  py::class_<data::MagneticSpacegroupType>(
      m, "MagneticSpacegroupType",
      "One UNI number's metadata: the BNS and OG symbols, the family "
      "space-group number and the construction type.")
      .def_readonly("uni_number", &data::MagneticSpacegroupType::uni_number)
      .def_readonly("litvin_number",
                    &data::MagneticSpacegroupType::litvin_number)
      .def_property_readonly("bns_number",
                             [](data::MagneticSpacegroupType const &self) {
                               return to_str(self.bns_number);
                             })
      .def_property_readonly("og_number",
                             [](data::MagneticSpacegroupType const &self) {
                               return to_str(self.og_number);
                             })
      .def_readonly("number", &data::MagneticSpacegroupType::number,
                    "Family space-group international number, 1..230.")
      .def_readonly("type", &data::MagneticSpacegroupType::type,
                    "Construction type 1..4.")
      .def("__repr__", [](data::MagneticSpacegroupType const &self) {
        return "MagneticSpacegroupType(" + std::to_string(self.uni_number) +
               ", '" + std::string(self.bns_number) + "')";
      });

  m.def(
      "magnetic_spacegroup_type",
      [](UniNumber uni) -> data::MagneticSpacegroupType const & {
        return data::magnetic_spacegroup_type(uni);
      },
      py::arg("uni"), py::return_value_policy::reference,
      py::doc("The metadata of a UNI number."));
  m.def(
      "uni_candidates",
      [](HallNumber hall) { return data::uni_candidates(hall); },
      py::arg("hall"),
      py::doc("The (first, last) UNI numbers a 3D Hall setting can carry."));
  m.def("magnetic_operations_from_database",
        &data::magnetic_operations_from_database, py::arg("uni"),
        py::arg("hall") = std::nullopt, py::return_value_policy::copy,
        py::doc("The operations of a UNI number in a Hall setting (None = its "
                "first); empty when the pairing is invalid."));
  m.def("magnetic_std_transformations", &data::magnetic_std_transformations,
        py::arg("uni"), py::arg("hall") = std::nullopt,
        py::return_value_policy::copy,
        py::doc("The alternative standardized-setting transformations of a "
                "UNI number, identity first."));
  m.attr("K_NUM_UNI_NUMBERS") = data::kNumUniNumbers;
}

} // namespace seitz::python
