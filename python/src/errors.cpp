#include "errors.hpp"

#include <pybind11/pybind11.h>

#include <string>

namespace cppcrystal::python {

namespace {

// The classes themselves, built with the C API rather than py::exception<T>.
// PyErr_NewExceptionWithDoc makes a plain heap type -- so its instances accept
// arbitrary attributes, which is how the payloads below reach Python -- and it
// takes the docstring at creation instead of needing one patched on afterwards.
//
// The returned reference is deliberately never released: these are
// interpreter-lifetime type objects, and letting a static py::object decref one
// at static-destruction time would run after the GIL is already gone.
[[nodiscard]] py::handle make(py::module_ &m, char const *name, char const *doc,
                              py::handle base) {
  std::string const qualified = std::string("cppcrystal._core.") + name;
  py::handle const type =
      PyErr_NewExceptionWithDoc(qualified.c_str(), doc, base.ptr(), nullptr);
  if (!type) {
    throw py::error_already_set();
  }
  m.add_object(name, type);
  return type;
}

ErrorTypes g_types;

} // namespace

namespace detail {

void raise(py::handle type, char const *text) {
  py::set_error(type, text);
  throw py::error_already_set();
}

void raise(py::handle type, char const *text, char const *name, double value) {
  // Built by hand rather than through py::set_error(type, text), so the payload
  // can be attached to the instance the caller will actually catch.
  py::object const error = py::reinterpret_steal<py::object>(
      PyObject_CallFunction(type.ptr(), "s", text));
  if (!error) {
    throw py::error_already_set();
  }
  py::setattr(error, name, py::float_(value));
  PyErr_SetObject(type.ptr(), error.ptr());
  throw py::error_already_set();
}

} // namespace detail

ErrorTypes const &error_types() noexcept { return g_types; }

void register_errors(py::module_ &m) {
  g_types.base =
      make(m, "CppCrystalError",
           "Base of every error this library reports. The C++ side returns "
           "Result<T> and never throws; this is what that becomes.",
           PyExc_Exception);

  g_types.spacegroup_search_failed =
      make(m, "SpacegroupSearchFailedError", "No Hall setting fits this cell.",
           g_types.base);
  g_types.cell_standardization_failed =
      make(m, "CellStandardizationFailedError",
           "The cell could not be brought into its standardized setting.",
           g_types.base);
  g_types.symmetry_operation_search_failed = make(
      m, "SymmetryOperationSearchFailedError",
      "The search for the cell's symmetry operations failed.", g_types.base);
  g_types.magnetic_symmetry_search_failed =
      make(m, "MagneticSymmetrySearchFailedError",
           "The magnetic symmetry search failed.", g_types.base);
  g_types.pointgroup_not_found = make(
      m, "PointgroupNotFoundError",
      "The operations match no crystallographic point group.", g_types.base);
  g_types.niggli_failed =
      make(m, "NiggliFailedError", "Niggli reduction did not converge.",
           g_types.base);
  g_types.delaunay_failed =
      make(m, "DelaunayFailedError",
           "Delaunay reduction failed: a degenerate cell, or a change of basis "
           "that is not unimodular.",
           g_types.base);
  g_types.empty_cell =
      make(m, "EmptyCellError", "The cell has no atoms.", g_types.base);

  // The two payload-carrying classes. The attribute is the whole reason to
  // catch one of these rather than the base, so the docstring names it.
  g_types.invalid_lattice =
      make(m, "InvalidLatticeError",
           "A (near-)singular lattice basis, rejected before its inverse could "
           "propagate NaN through the pipeline. Carries .determinant.",
           g_types.base);
  g_types.atoms_too_close =
      make(m, "AtomsTooCloseError",
           "Two atoms closer than the distance tolerance allows. Carries "
           ".distance.",
           g_types.base);
}

} // namespace cppcrystal::python
