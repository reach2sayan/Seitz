#include <cppcrystal/core/version.hpp>

#include "errors.hpp" // register_errors

#include <pybind11/pybind11.h>

// The extension module. Nothing is bound here: each subsystem's translation
// unit binds its own headers, mirroring include/cppcrystal/, and this file only
// fixes the order they run in.
//
// _core is private. The public surface is the pure-Python cppcrystal package,
// which re-exports what belongs in it and adds the pydantic models -- the
// things a C++ binding cannot express well (keyword-only arguments, overloads,
// a stable __all__, JSON round-trips).
namespace cppcrystal::python {

void bind_core(py::module_ &m);
void bind_core_symmetry(py::module_ &m);
void bind_analysis(py::module_ &m);
void bind_group(py::module_ &m);
void bind_data(py::module_ &m);

} // namespace cppcrystal::python

PYBIND11_MODULE(_core, m) {
  namespace ccp = cppcrystal::python;

  m.doc() = "Raw bindings for CppCrystal. Import cppcrystal instead.";
  m.attr("__version__") = cppcrystal::version_string();

  // Errors first: every other binding's unwrap() reaches for these types, so
  // they have to exist before any of them can raise.
  ccp::register_errors(m);

  // Then vocabulary before the things phrased in it -- Cell and Tolerance have
  // to be registered before SymmetryAnalyzer names them in a signature, or the
  // generated docstrings fall back to raw C++ type names.
  ccp::bind_core(m);
  ccp::bind_core_symmetry(m);
  ccp::bind_data(m);
  ccp::bind_group(m);
  ccp::bind_analysis(m);
}
