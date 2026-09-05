#include <seitz/core/cell.hpp>
#include <seitz/core/keys.hpp>
#include <seitz/core/lattice.hpp>
#include <seitz/core/periodicity.hpp>
#include <seitz/core/tolerance.hpp>
#include <seitz/core/types.hpp>
#include <seitz/core/version.hpp>
#include <seitz/warmup.hpp>

#include "casters.hpp" // borrowed_list, types_array, to_str
#include "errors.hpp"  // unwrap

#include <pybind11/eigen.h>
#include <pybind11/native_enum.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace seitz::python {

namespace {

// Lattice's checked door. The unchecked `explicit Lattice(Matrix3d)` is not
// bound: from_basis is what rejects a singular basis, and a Python constructor
// that silently accepted one would hand NaN to every later query.
[[nodiscard]] Lattice lattice_from(Matrix3d const &basis) {
  if (auto lattice = Lattice::from_basis(basis)) {
    return *std::move(lattice);
  }
  detail::raise(error_types().invalid_lattice, "the lattice basis is singular",
                "determinant", basis.determinant());
}

[[nodiscard]] std::string lattice_repr(Lattice const &self) {
  return "Lattice(volume=" + std::to_string(self.volume()) + ")";
}

[[nodiscard]] std::string cell_repr(Cell const &self) {
  return "Cell(" + std::to_string(self.size()) + " atoms)";
}

} // namespace

void bind_core(py::module_ &m) {
  // ---- enums -------------------------------------------------------------
  //
  // Real enum.IntEnum classes, not pybind11's enum type: Python prints, hashes,
  // aliases, pickles and pattern-matches those already. Member names keep the
  // C++ lower_snake_case against Python's UPPER_CASE convention -- a rename
  // would be a second vocabulary for one concept.
  py::native_enum<GroupFamily>(m, "GroupFamily", "enum.IntEnum",
                               "The two families the determination handles: "
                               "3D space groups, and 2D-periodic layer groups.")
      .value("space", GroupFamily::space)
      .value("layer", GroupFamily::layer)
      .finalize();

  py::native_enum<AxisKind>(m, "AxisKind", "enum.IntEnum",
                            "Whether a cell axis repeats.")
      .value("periodic", AxisKind::periodic)
      .value("aperiodic", AxisKind::aperiodic)
      .finalize();

  // IntFlag, so `Warm.space_groups | Warm.layer_groups` works in Python
  // without binding the C++ operator|.
  py::native_enum<Warm>(m, "Warm", "enum.IntFlag",
                        "Which per-setting caches warmup() should prime.")
      .value("space_groups", Warm::space_groups)
      .value("layer_groups", Warm::layer_groups)
      .value("all", Warm::all)
      .finalize();

  // ---- periodicity -------------------------------------------------------
  //
  // CellPeriodicity is std::array<AxisKind, 3>, which pybind11/stl.h already
  // converts to and from a 3-tuple, so only the named constructors need
  // binding. They are what a caller should reach for -- the tuple spelling is
  // the escape hatch, not the API.
  m.def("all_periodic", &all_periodic,
        py::doc("The fully periodic descriptor: a 3D space-group cell."));
  m.def("aperiodic_along", &aperiodic_along, py::arg("axis"),
        py::doc("Periodic in the plane, aperiodic along `axis` -- a layer."));
  m.def("periodic_along", &periodic_along, py::arg("axis"),
        py::doc("Periodic along `axis` only -- a rod."));
  m.def("none_periodic", &none_periodic,
        py::doc("No periodic axis at all -- a 0D cluster."));
  m.def("family_of", &family_of, py::arg("periodicity"),
        py::doc("The family a periodicity puts a cell in."));
  m.def("aperiodic_axis", &aperiodic_axis, py::arg("periodicity"),
        py::doc("The single aperiodic axis, or None when there is not exactly "
                "one."));
  m.def("minimal_image", &minimal_image, py::arg("diff"),
        py::arg("periodicity"),
        py::doc("`diff` folded to its nearest-image residue along every "
                "periodic axis."));
  m.def("wrap", &wrap, py::arg("vector"), py::arg("periodicity"),
        py::doc("A fractional coordinate folded into [0, 1) along every "
                "periodic axis."));

  // ---- validated keys ----------------------------------------------------
  //
  // HallNumber has a private constructor and a static of() returning optional.
  // Both Python idioms are bound, uniformly for every such type here:
  // `X.of(...)` answers None, `X(...)` raises -- as re.match and int("x") do.
  py::class_<HallNumber>(m, "HallNumber",
                         "A validated Hall setting: a family plus a 1-based "
                         "index into that family's settings.")
      .def(py::init([](GroupFamily family, int index) {
             if (auto const hall = HallNumber::of(family, index)) {
               return *hall;
             }
             detail::raise(error_types().base,
                           "hall number out of range for this family");
           }),
           py::arg("family"), py::arg("index"))
      .def_static(
          "of",
          [](GroupFamily family, int index) {
            return HallNumber::of(family, index);
          },
          py::arg("family"), py::arg("index"),
          py::doc("The setting, or None when `index` is out of range."))
      .def_property_readonly("family", &HallNumber::family)
      .def_property_readonly("index", &HallNumber::index)
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self < py::self)
      .def(py::self <= py::self)
      .def(py::self > py::self)
      .def(py::self >= py::self)
      .def("__hash__",
           [](HallNumber const &self) {
             return py::hash(py::make_tuple(self.family(), self.index()));
           })
      // setstate returns by value, so it can go through of() and never needs
      // access to the private constructor.
      .def(py::pickle(
          [](HallNumber const &self) {
            return py::make_tuple(self.family(), self.index());
          },
          [](py::tuple const &t) {
            auto const hall =
                HallNumber::of(t[0].cast<GroupFamily>(), t[1].cast<int>());
            if (!hall) {
              detail::raise(error_types().base, "unpickling an invalid hall "
                                                "number");
            }
            return *hall;
          }))
      .def("__repr__", [](HallNumber const &self) {
        return "HallNumber(GroupFamily." +
               std::string(self.family() == GroupFamily::layer ? "layer"
                                                               : "space") +
               ", " + std::to_string(self.index()) + ")";
      });

  py::class_<UniNumber>(m, "UniNumber",
                        "A validated magnetic space-group (UNI) number, "
                        "1..1651.")
      .def(py::init([](int number) {
             if (auto const uni = UniNumber::of(number)) {
               return *uni;
             }
             detail::raise(error_types().base, "uni number out of range");
           }),
           py::arg("number"))
      .def_static(
          "of", [](int number) { return UniNumber::of(number); },
          py::arg("number"), py::doc("The number, or None when out of range."))
      .def_property_readonly("value", &UniNumber::value)
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self < py::self)
      .def(py::self <= py::self)
      .def(py::self > py::self)
      .def(py::self >= py::self)
      .def("__hash__",
           [](UniNumber const &self) {
             return py::hash(py::int_(self.value()));
           })
      .def(py::pickle(
          [](UniNumber const &self) { return py::make_tuple(self.value()); },
          [](py::tuple const &t) {
            auto const uni = UniNumber::of(t[0].cast<int>());
            if (!uni) {
              detail::raise(error_types().base, "unpickling an invalid uni "
                                                "number");
            }
            return *uni;
          }))
      .def("__repr__", [](UniNumber const &self) {
        return "UniNumber(" + std::to_string(self.value()) + ")";
      });

  m.attr("K_SPACE_HALL_SETTINGS") = kSpaceHallSettings;
  m.attr("K_LAYER_HALL_SETTINGS") = kLayerHallSettings;
  m.attr("K_UNI_NUMBERS") = kUniNumbers;

  // ---- tolerances --------------------------------------------------------
  //
  // Bound as a plain struct even though the pure-Python layer models it with
  // pydantic: the model validates and then hands one of these across, so this
  // is the boundary type, not the user-facing one.
  py::class_<Tolerance>(m, "Tolerance",
                        "The tolerances threaded through the symmetry search.")
      .def(py::init([](double symprec, AngleTolerance angle_tolerance) {
             return Tolerance{symprec, angle_tolerance};
           }),
           py::arg("symprec") = kDefaultSymprec,
           py::arg("angle_tolerance") = std::nullopt)
      .def_readwrite("symprec", &Tolerance::symprec)
      .def_readwrite("angle_tolerance", &Tolerance::angle_tolerance)
      .def("__repr__", [](Tolerance const &self) {
        return "Tolerance(symprec=" + std::to_string(self.symprec) + ")";
      });

  py::class_<MagneticTolerance, Tolerance>(
      m, "MagneticTolerance",
      "Tolerance plus the magnetic search's moment comparison.")
      .def(py::init([](double symprec, AngleTolerance angle_tolerance,
                       std::optional<double> moment) {
             MagneticTolerance tol;
             tol.symprec = symprec;
             tol.angle_tolerance = angle_tolerance;
             tol.moment = moment;
             return tol;
           }),
           py::arg("symprec") = kDefaultSymprec,
           py::arg("angle_tolerance") = std::nullopt,
           py::arg("moment") = std::nullopt)
      .def_readwrite("moment", &MagneticTolerance::moment)
      .def_property_readonly("moment_or_symprec",
                             &MagneticTolerance::moment_or_symprec);

  m.attr("K_DEFAULT_SYMPREC") = kDefaultSymprec;
  m.attr("K_ZERO_PREC") = kZeroPrec;

  // ---- lattice -----------------------------------------------------------
  py::class_<Lattice>(m, "Lattice",
                      "A crystal lattice: the three basis vectors as the "
                      "COLUMNS of a 3x3 matrix.")
      .def(py::init(&lattice_from), py::arg("basis"),
           py::doc("Build from a 3x3 basis whose columns are the basis "
                   "vectors. Raises InvalidLatticeError if it is singular."))
      .def_static(
          "from_basis",
          [](Matrix3d const &basis) { return Lattice::from_basis(basis); },
          py::arg("basis"),
          py::doc("The lattice, or None for a (near-)singular basis."))
      // A copy, not a view into basis_: Lattice is immutable on the C++ side,
      // and reference_internal would hand Python a writable alias of it.
      // Matrix3d is column-major, so this arrives Fortran-ordered and
      // `lattice.matrix[:, 0]` really is basis vector a, with no transpose.
      .def_property_readonly(
          "matrix", [](Lattice const &self) { return self.matrix(); },
          py::doc("The 3x3 basis, columns = basis vectors (F-ordered)."))
      .def_property_readonly(
          "a",
          [](Lattice const &self) { return Vector3d(self.matrix().col(0)); },
          py::doc("First basis vector."))
      .def_property_readonly(
          "b",
          [](Lattice const &self) { return Vector3d(self.matrix().col(1)); },
          py::doc("Second basis vector."))
      .def_property_readonly(
          "c",
          [](Lattice const &self) { return Vector3d(self.matrix().col(2)); },
          py::doc("Third basis vector."))
      .def_property_readonly("metric", &Lattice::metric,
                             py::doc("The metric (Gram) tensor L^T L."))
      .def_property_readonly("volume", &Lattice::volume)
      .def("to_cartesian", &Lattice::to_cartesian, py::arg("fractional"))
      .def("to_fractional", &Lattice::to_fractional, py::arg("cartesian"))
      .def("transformed", &Lattice::transformed, py::arg("t"),
           py::doc("The lattice whose columns are L . t."))
      .def("rigid_rotation_to", &Lattice::rigid_rotation_to, py::arg("ideal"))
      .def(
          "niggli",
          [](Lattice const &self, double eps) {
            return unwrap([&] { return self.niggli(eps); });
          },
          py::arg("eps") = kDefaultSymprec,
          py::doc("Krivy-Gruber Niggli reduction."))
      .def(
          "delaunay",
          [](Lattice const &self, double symprec) {
            return unwrap([&] { return self.delaunay(symprec); });
          },
          py::arg("symprec") = kDefaultSymprec,
          py::doc("Delaunay reduction; the result is right-handed."))
      .def(
          "delaunay_in_plane",
          [](Lattice const &self, int unique_axis, double symprec) {
            return unwrap(
                [&] { return self.delaunay_in_plane(unique_axis, symprec); });
          },
          py::arg("unique_axis"), py::arg("symprec") = kDefaultSymprec,
          py::doc("2D Delaunay reduction in the plane normal to "
                  "`unique_axis`."))
      .def(py::pickle(
          [](Lattice const &self) { return py::make_tuple(self.matrix()); },
          [](py::tuple const &t) {
            return lattice_from(t[0].cast<Matrix3d>());
          }))
      .def("__copy__", [](Lattice const &self) { return self; })
      .def(
          "__deepcopy__",
          [](Lattice const &self, py::dict const &) { return self; },
          py::arg("memo"))
      .def("__repr__", &lattice_repr);

  // ---- cell --------------------------------------------------------------
  py::class_<Cell>(m, "Cell",
                   "A crystal cell: a lattice, fractional positions (row i is "
                   "atom i) and integer types, plus the per-axis periodicity.")
      .def(py::init<Lattice, Positions, Types, CellPeriodicity>(),
           py::arg("lattice"), py::arg("positions"), py::arg("types"),
           py::arg_v("periodicity", all_periodic(), "all_periodic()"))
      .def_property_readonly("lattice", &Cell::lattice,
                             py::return_value_policy::copy)
      // Positions is Dynamic x 3 RowMajor, so this is a C-contiguous (N, 3)
      // float64 array. One copy, because the accessor hands back a reference
      // into an immutable Cell.
      .def_property_readonly(
          "positions", [](Cell const &self) { return self.positions(); },
          py::doc("Fractional positions as a C-contiguous (N, 3) float64 "
                  "array."))
      .def_property_readonly(
          "types", [](Cell const &self) { return types_array(self.types()); },
          py::doc("Atom types as an (N,) int32 array."))
      .def_property_readonly("periodicity", &Cell::periodicity)
      .def("position", &Cell::position, py::arg("i"))
      .def("type", &Cell::type, py::arg("i"))
      // Cell::atoms() is an iota|transform view that captures `this`, and its
      // iterators hold a pointer to the view object -- so py::make_iterator
      // over the temporary returned here would hand Python an iterator into
      // storage that is already gone. Materialised instead; a caller who wants
      // the columns rather than the rows has .positions and .types, which are
      // the arrays themselves.
      .def_property_readonly(
          "atoms",
          [](Cell const &self) {
            py::list atoms;
            for (auto const &[position, type] : self.atoms()) {
              atoms.append(py::make_tuple(position, type));
            }
            return atoms;
          },
          py::doc("The atoms as a list of (fractional position, type) pairs."))
      .def("with_lattice", &Cell::with_lattice, py::arg("lattice"))
      .def("with_periodicity", &Cell::with_periodicity, py::arg("periodicity"))
      .def("__len__", &Cell::size)
      .def("__getitem__",
           [](Cell const &self, Index i) {
             if (i < 0) {
               i += self.size();
             }
             if (i < 0 || i >= self.size()) {
               throw py::index_error("cell index out of range");
             }
             return py::make_tuple(self.position(i), self.type(i));
           })
      .def("__iter__",
           [](py::object const &self) { return py::iter(self.attr("atoms")); })
      .def(py::pickle(
          [](Cell const &self) {
            return py::make_tuple(self.lattice(), self.positions(),
                                  self.types(), self.periodicity());
          },
          [](py::tuple const &t) {
            return Cell{t[0].cast<Lattice>(), t[1].cast<Positions>(),
                        t[2].cast<Types>(), t[3].cast<CellPeriodicity>()};
          }))
      .def("__copy__", [](Cell const &self) { return self; })
      .def(
          "__deepcopy__",
          [](Cell const &self, py::dict const &) { return self; },
          py::arg("memo"))
      .def("__repr__", &cell_repr);

  m.def(
      "to_positions",
      [](std::vector<Vector3d> const &rows) {
        return to_positions(std::span<Vector3d const>{rows});
      },
      py::arg("rows"),
      py::doc("Pack a sequence of 3-vectors into an (N, 3) block."));

  // ---- version and warmup ------------------------------------------------
  py::class_<Version>(m, "Version", "A major.minor.patch triple.")
      .def_readonly("major", &Version::major)
      .def_readonly("minor", &Version::minor)
      .def_readonly("patch", &Version::patch)
      .def("__repr__", [](Version const &v) {
        return "Version(" + std::to_string(v.major) + ", " +
               std::to_string(v.minor) + ", " + std::to_string(v.patch) + ")";
      });

  m.attr("K_VERSION") = kVersion;
  m.attr("K_REFERENCE_SPGLIB_VERSION") = kReferenceSpglibVersion;
  m.def("version_string", &version_string,
        py::doc("\"major.minor.patch\" of this port."));

  // warmup_async is deliberately not bound: it returns std::future<void>, which
  // has no caster, and waiting on one would hold the GIL for the duration.
  // concurrent.futures.ThreadPoolExecutor().submit(warmup) is the Python answer
  // and composes with the rest of the language.
  m.def(
      "warmup", [](Warm what) { warmup(what); }, py::arg("what") = Warm::all,
      py::call_guard<py::gil_scoped_release>(),
      py::doc("Build the requested per-setting caches ahead of first use. An "
              "optimization, never a precondition."));
}

} // namespace seitz::python
