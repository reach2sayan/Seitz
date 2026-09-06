#include <seitz/core/version.hpp>

#include "errors.hpp" // register_errors

#include <pybind11/pybind11.h>

// The extension module. Nothing is bound here: each subsystem's TU binds its
// own headers, mirroring include/seitz/, and this file only fixes their order.
//
// _core is private; the public surface is the pure-Python seitz package, which
// re-exports from it and adds what a C++ binding expresses poorly (keyword-only
// arguments, overloads, a stable __all__, JSON round-trips).
namespace seitz::python {

void bind_core(py::module_ &m);
void bind_core_symmetry(py::module_ &m);
void bind_analysis(py::module_ &m);
void bind_group(py::module_ &m);
void bind_data(py::module_ &m);
void bind_io(py::module_ &m);

} // namespace seitz::python

PYBIND11_MODULE(_core, m) {
  namespace sp = seitz::python;

  m.doc() = "Raw bindings for Seitz. Import seitz instead.";
  m.attr("__version__") = seitz::version_string();

  // Errors first: every other binding's unwrap() reaches for these types, so
  // they have to exist before any of them can raise.
  sp::register_errors(m);

  // Then vocabulary before the things phrased in it -- Cell and Tolerance have
  // to be registered before SymmetryAnalyzer names them in a signature, or the
  // generated docstrings fall back to raw C++ type names.
  sp::bind_core(m);
  sp::bind_core_symmetry(m);
  sp::bind_data(m);
  sp::bind_group(m);
  sp::bind_analysis(m);

  // Last: the CIF layer names Cell, Tolerance and SymmetryAnalyzer in its
  // signatures, so all three have to be registered before it.
  sp::bind_io(m);
}
