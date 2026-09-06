#pragma once

#include <seitz/core/error.hpp>

#include <pybind11/pybind11.h>

#include <concepts>
#include <type_traits>
#include <utility>

// The one place a Result<T> becomes a Python exception, and the one helper that
// binds an Analyzer projection. Everything else in this module is ordinary
// pybind11.
namespace seitz::python {

namespace py = pybind11;

// The exception classes, built once at import. Held by handle rather than
// py::object: these are interpreter-lifetime type objects, and a static
// py::object would decref them during static destruction, after the interpreter
// -- and the GIL -- are already gone.
struct ErrorTypes {
  py::handle base; // SeitzError, which every other one derives from
  py::handle spacegroup_search_failed;
  py::handle cell_standardization_failed;
  py::handle symmetry_operation_search_failed;
  py::handle magnetic_symmetry_search_failed;
  py::handle pointgroup_not_found;
  py::handle niggli_failed;
  py::handle delaunay_failed;
  py::handle empty_cell;
  py::handle invalid_lattice;
  py::handle atoms_too_close;
  py::handle cif_syntax;
  py::handle cif_missing_tag;
  py::handle invalid_xyz;
  py::handle unknown_element;
  py::handle unknown_spacegroup_symbol;
};

void register_errors(py::module_ &m);
[[nodiscard]] ErrorTypes const &error_types() noexcept;

namespace detail {

// Set `type` as the pending Python exception and throw. The payload overload
// also puts the number on the instance, so `except AtomsTooCloseError as e:
// e.distance` reads the double the search computed instead of parsing a
// sentence back apart.
[[noreturn]] void raise(py::handle type, char const *text);
[[noreturn]] void raise(py::handle type, char const *text, char const *name,
                        double value);
// The same, for a payload that is not a number: a tag, a symbol, the text that
// would not parse. Takes an already-built object so one call site can attach a
// string and another an int.
[[noreturn]] void raise(py::handle type, char const *text, char const *name,
                        py::object value);

// The message LEAF carried, or the tag's own default when it carried none.
[[nodiscard]] inline char const *message_or(e_message const *m,
                                            char const *fallback) noexcept {
  return m != nullptr ? m->text.c_str() : fallback;
}

} // namespace detail

// The success value of `make()`, or a raised Python exception carrying the
// error.
//
// Takes the CALL, not its Result: LEAF keeps an error's payloads in
// thread-local slots live only while a matching context is active, so a Result
// produced before try_handle_all is entered arrives with its payloads gone and
// every typed error degrades to the bare fallback.
//
// Nothing in the library throws -- every fallible entry point returns Result<T>
// with a typed tag -- so this is the one translation the binding layer needs,
// raising directly rather than through C++ exception classes that would exist
// only to be caught two lines later.
//
// The value comes back BY VALUE, always: Result<T const &> holds a
// reference_wrapper into an analyzer's memo, and returning that reference from
// try_handle_all re-binds it through its own frame. Result<Cell const &> ->
// Result<Cell> is LEAF's converting move constructor.
//
// Handlers take `e_message const *`, not `const &`: LEAF reads a pointer
// parameter as "and this one if also present", so the specific class and the
// specific sentence both survive instead of one shadowing the other.
template <ResultProducer F> [[nodiscard]] auto unwrap(F &&make) {
  using Value =
      std::remove_cvref_t<typename std::invoke_result_t<F &>::value_type>;
  using detail::message_or;
  using detail::raise;
  ErrorTypes const &types = error_types();

  return leaf::try_handle_all(
      [&]() -> Result<Value> { return make(); },
      [&](e_invalid_lattice const &e, e_message const *m) -> Value {
        raise(types.invalid_lattice,
              message_or(m, "the lattice basis is singular"), "determinant",
              e.determinant);
      },
      [&](e_atoms_too_close const &e, e_message const *m) -> Value {
        raise(types.atoms_too_close,
              message_or(m, "atoms overlap within the tolerance"), "distance",
              e.distance);
      },
      [&](e_cif_syntax const &e, e_message const *m) -> Value {
        py::object const error = py::reinterpret_steal<py::object>(
            PyObject_CallFunction(types.cif_syntax.ptr(), "s",
                                  message_or(m, "the CIF could not be read")));
        if (!error) {
          throw py::error_already_set();
        }
        py::setattr(error, "line", py::int_(e.line));
        py::setattr(error, "column", py::int_(e.column));
        PyErr_SetObject(types.cif_syntax.ptr(), error.ptr());
        throw py::error_already_set();
      },
      [&](e_cif_missing const &e, e_message const *m) -> Value {
        raise(types.cif_missing_tag,
              message_or(m, "the block does not carry a tag the reader needs"),
              "tag", py::str(e.tag));
      },
      [&](e_invalid_xyz const &e, e_message const *m) -> Value {
        raise(types.invalid_xyz,
              message_or(m, "not a coordinate triplet"), "text",
              py::str(e.text));
      },
      [&](e_unknown_element const &e, e_message const *m) -> Value {
        raise(types.unknown_element,
              message_or(m, "no tabulated element has this symbol"), "symbol",
              py::str(e.symbol));
      },
      [&](e_unknown_spacegroup_symbol const &e, e_message const *m) -> Value {
        raise(types.unknown_spacegroup_symbol,
              message_or(m, "no tabulated setting has this symbol"), "symbol",
              py::str(e.symbol));
      },
      [&](e_empty_cell const &, e_message const *m) -> Value {
        raise(types.empty_cell, message_or(m, "the cell has no atoms"));
      },
      [&](e_spacegroup_search_failed const &, e_message const *m) -> Value {
        raise(types.spacegroup_search_failed,
              message_or(m, "no Hall setting fits this cell"));
      },
      [&](e_cell_standardization_failed const &, e_message const *m) -> Value {
        raise(types.cell_standardization_failed,
              message_or(m, "the cell could not be standardized"));
      },
      [&](e_symmetry_operation_search_failed const &,
          e_message const *m) -> Value {
        raise(types.symmetry_operation_search_failed,
              message_or(m, "the symmetry operation search failed"));
      },
      [&](e_magnetic_symmetry_search_failed const &,
          e_message const *m) -> Value {
        raise(types.magnetic_symmetry_search_failed,
              message_or(m, "the magnetic symmetry search failed"));
      },
      [&](e_pointgroup_not_found const &, e_message const *m) -> Value {
        raise(types.pointgroup_not_found,
              message_or(m, "the operations match no crystallographic point "
                            "group"));
      },
      [&](e_niggli_failed const &, e_message const *m) -> Value {
        raise(types.niggli_failed,
              message_or(m, "Niggli reduction did not converge"));
      },
      [&](e_delaunay_failed const &, e_message const *m) -> Value {
        raise(types.delaunay_failed,
              message_or(m, "Delaunay reduction failed"));
      },
      // A message with no tag beside it: the group catalogs raise these.
      [&](e_message const &m) -> Value { raise(types.base, m.text.c_str()); },
      [&](leaf::error_info const &) -> Value {
        raise(types.base, "seitz: unclassified failure");
      });
}

// A const&-qualified Result<T const &> accessor, as something pybind11 can
// bind. The parameter type is the point: every Analyzer projection is a pair --
// `f() const &` and `f() const && = delete` -- so a bare &SymmetryAnalyzer::hall
// is an ambiguous overload set. Naming the const& signature resolves it once
// instead of a static_cast per binding site.
//
// The GIL is dropped around the call and the copy out: determination is the
// expensive operation, detail::Lazy makes a const analyzer shareable across
// threads, and nothing below touches a Python object. unwrap() then runs with
// the GIL held, since raising needs it. py::call_guard<gil_scoped_release>
// would be wrong -- it wraps the return cast too, and unwrap() raises.
template <class Self, class T>
[[nodiscard]] auto memo(Result<T const &> (Self::*accessor)() const &) {
  return [accessor](Self const &self) -> T {
    return unwrap([&]() -> Result<T> {
      // Released around the call itself, inside the LEAF context. The handlers
      // run after this scope closes, so they hold the GIL again by the time one
      // of them needs to raise.
      py::gil_scoped_release const unlocked;
      return Result<T>{(self.*accessor)()};
    });
  };
}

// The same, for the one projection that was never ref-qualified because it
// hands back a reference into a constexpr catalog rather than into the memo.
template <class Self, class T>
[[nodiscard]] auto memo(Result<T const &> (Self::*accessor)() const) {
  return [accessor](Self const &self) -> T {
    return unwrap([&]() -> Result<T> {
      // Released around the call itself, inside the LEAF context. The handlers
      // run after this scope closes, so they hold the GIL again by the time one
      // of them needs to raise.
      py::gil_scoped_release const unlocked;
      return Result<T>{(self.*accessor)()};
    });
  };
}

// The same again, naming the class to bind against, for a projection the
// analyzer INHERITS -- dataset() is declared on the CRTP base, so plain memo()
// would deduce Self as Analyzer<Derived, Traits> and hand pybind11 a lambda
// taking a base-class reference it has no caster for. (A member pointer would
// have been fine: class_::def runs method_adaptor over one of those. A lambda
// gets no such help, so the derived type is named here instead.)
template <class Derived, class Base, class T>
  requires std::derived_from<Derived, Base>
[[nodiscard]] auto memo_as(Result<T const &> (Base::*accessor)() const &) {
  return [accessor](Derived const &self) -> T {
    return unwrap([&]() -> Result<T> {
      // Released around the call itself, inside the LEAF context. The handlers
      // run after this scope closes, so they hold the GIL again by the time one
      // of them needs to raise.
      py::gil_scoped_release const unlocked;
      return Result<T>{(self.*accessor)()};
    });
  };
}

} // namespace seitz::python
