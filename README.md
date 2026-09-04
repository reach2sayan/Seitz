# CppCrystal

[![C++23](https://img.shields.io/badge/C%2B%2B-23-blue.svg)](https://en.cppreference.com/w/cpp/23)
[![CMake](https://img.shields.io/badge/CMake-3.28%2B-064F8C.svg?logo=cmake&logoColor=white)](https://cmake.org)
[![Compilers](https://img.shields.io/badge/compilers-GCC%2015%2B%20%7C%20Clang-brightgreen.svg)](#building)
[![License](https://img.shields.io/badge/license-BSD--3--Clause-green.svg)](LICENSE)
[![Validated against spglib](https://img.shields.io/badge/validated%20against-spglib%20v2.7.0-orange.svg)](https://github.com/spglib/spglib)

A modern C++23 crystallography library for symmetry analysis and structure
generation. It determines the space group of a crystal, standardizes and refines
cells, classifies magnetic structures, reduces reciprocal-space meshes to their
irreducible wedge, and generates random crystals, layers, rods, and clusters
from a chosen symmetry group.

> **License & attribution.** CppCrystal is released under the BSD-3-Clause
> license (see [`LICENSE`](LICENSE)). Its algorithms are derived from
> [**spglib**](https://github.com/spglib/spglib) (BSD-3-Clause) and its
> object-oriented, generation-focused API is modeled on
> [**PyXtal**](https://github.com/qzhu2017/PyXtal) (MIT). Full credits are in
> [Attribution](#attribution).

---

## The API

One umbrella header reaches everything, and it is also everything that ships:

```cpp
#include <cppcrystal/cppcrystal.hpp>
```

| Capability | Entry point |
|---|---|
| Space-group determination | `analysis::SymmetryAnalyzer::from_cell(cell)` |
| Layer groups | the same analyzer — a `Cell` with one aperiodic axis takes the layer path |
| Cell standardization | `analyzer.standardized_cell<CellSetting, Idealize>()` |
| Magnetic structures | `analysis::MagneticSymmetryAnalyzer::from_cell(mcell)` |
| Lattice reduction | `Lattice::niggli`, `::delaunay`, `::delaunay_in_plane` |
| Group catalogs | `group::SpaceGroup::of` / `::from_number`, `group::PointGroup`, `group::RodGroup` |
| Wyckoff positions | `group.wyckoffs()` → `std::span<group::Wyckoff const>` |
| Subgroup relations | `group::SubgroupGraph::maximal_subgroups`, `::is_subgroup`, `::path` |
| Structure generation | `generate::Generator<G>` |
| Reciprocal-space meshes | `kpoint::Mesh`, `kpoint::ReciprocalMesh`, `kpoint::BrillouinZone` |
| Symmetry tables | `data::spacegroup_type`, `data::element` |

---

## Quick start

### Determine a space group

`SymmetryAnalyzer` owns the cell and the tolerances and memoizes each stage, so
determination runs once and every later query is served from that cache:

```cpp
#include <cppcrystal/cppcrystal.hpp>
#include <print>
using namespace cppcrystal;

Cell const rutile{Lattice{basis}, positions, types};  // lattice columns = basis vectors
auto const sa = analysis::SymmetryAnalyzer::from_cell(rutile);

leaf::try_handle_all(
    [&]() -> Result<void> {
      BOOST_LEAF_AUTO(hall, sa.hall());        // runs determination
      BOOST_LEAF_AUTO(ops, sa.operations());   // reuses it
      BOOST_LEAF_AUTO(prim, sa.primitive_cell());

      auto const &type = data::spacegroup_type(hall);
      std::println("space group {} ({}), {} operations, {} primitive atoms",
                   type.number, type.international_short, ops.size(),
                   prim.size());
      return {};
    },
    [](e_spacegroup_search_failed) {
      std::println("no space group found at this tolerance");
    });
```

The analyzer is immutable after construction — to analyze at a different
tolerance, build a new one. Its queries hand back references into the memo, so
it must outlive them: the rvalue overloads are deleted, which makes
`from_cell(cell).hall()` a compile error rather than a dangling read.

`sa.determine()` returns a `Dataset` carrying the Hall key, the Bravais lattice,
the conventional setting, the operations of the input cell, one `Site` per atom
(Wyckoff letter, site-symmetry symbol, equivalence class) and the fully
standardized `Cell`.

### Query a group without a structure

```cpp
BOOST_LEAF_AUTO(fm3m, group::SpaceGroup::from_number(GroupFamily::space, 225));

for (group::Wyckoff const &w : fm3m->wyckoffs()) {
  std::println("{}{}  site symmetry {}, {} free parameter(s)",
               w.multiplicity(), w.letter(), w.site_symmetry(),
               w.degrees_of_freedom());
}
```

### Generate a random crystal

```cpp
generate::Composition const nacl{{11, 4}, {17, 4}};   // Na4 Cl4
BOOST_LEAF_AUTO(xtal, generate::Generator(*fm3m, {.seed = 42})(nacl));
// xtal.cell        -> the generated Cell
// xtal.assignment  -> which Wyckoff position each atom was seated on
```

`Generator` enumerates the valid ways to seat the composition on the group's
Wyckoff positions, then resamples free coordinates and the lattice until the
minimum-distance criterion is met. The same generator covers every family:
substitute a layer `SpaceGroup`, a `RodGroup` or a `PointGroup` and it supplies
that family's periodicity, seed window and random metric. The result is always a
`Cell` that describes itself — a cluster is simply a cell with no periodic axis.

### Irreducible k-points

```cpp
auto const mesh = *kpoint::Mesh::of({8, 8, 8});   // optional; rejects a non-positive mesh
BOOST_LEAF_AUTO(rmesh, sa.reciprocal_mesh(mesh, TimeReversal::on));

rmesh.num_irreducible();                  // count of irreducible points
rmesh.mapping();                          // each grid point's representative
rmesh.brillouin_zone(reciprocal_lattice); // Brillouin-zone folding
```

---

## What you can rely on

**Errors are values.** Nothing throws and no function returns a magic sentinel.
A fallible call returns `Result<T>` (`= boost::leaf::result<T>`); failures carry
typed context — `e_spacegroup_search_failed`, `e_invalid_lattice{det}`,
`e_magnetic_symmetry_search_failed`, and so on. "Absent" is `std::optional`,
never `-1`.

**Invalid states are unrepresentable, so most errors never arise.** A
`HallNumber`, a `UniNumber` and a `kpoint::Mesh` are each built through an
`of()` that is the only way in and the only place the check lives. A function
taking one of these keys is total — it cannot fail on a bad key, because you
cannot hold one.

**Cells describe periodicity per axis.** A `Cell` holds a column-major lattice,
fractional positions, integer types, and a `CellPeriodicity`. The same type
represents a 3D crystal, a 2D layer (one aperiodic axis), a 1D rod (two), and a
0D cluster (three) — the periodicity selects which symmetry path runs, so there
is no separate entry point per dimensionality.

**Groups are first-class objects.** `SpaceGroup`, `PointGroup` and `RodGroup`
are built from a group number and own their operations and Wyckoff positions as
queryable objects, not parallel arrays of letters and multiplicities.

**Value semantics throughout.** Every result type is a plain owning value; there
is no manual memory management and no ownership to reason about.

**Everything is safe to share.** A `const` analyzer can be used from any number
of threads from the moment it is built, the group catalogs are immutable and
race-free, and no query mutates observable state. `cppcrystal::warmup()` builds
every group setting up front (about 30 ms) to move first-use cost off the query
path; it is an optimization, never a precondition.

For the pipelines behind the API see [`docs/WORKFLOW.md`](docs/WORKFLOW.md), and
for the group theory [`docs/MATHEMATICS.md`](docs/MATHEMATICS.md).

---

## Building

Requires **GCC 15+ or Clang** with libstdc++ (the presets use `g++-15`) and
CMake ≥ 3.28; configuring with any other compiler is a hard error rather than an
untested branch. Everything else is fetched, unconditionally — there is no
system lookup and no `find_package` fallback to go stale: **Eigen 5.0.0**,
**Boost 1.88** (`container`, `flyweight`, `geometry`, `leaf`) and **Catch2 3**
for the tests. Nothing needs to be installed on the system.

`cppcrystal` builds as a **shared library**, versioned `0.1.0` with
`SONAME 0.1` — while the API is pre-1.0, every minor version may break ABI.

```bash
cmake --preset default              # configure (Debug)
cmake --build cmake-build-debug     # build the library + demo
ctest --test-dir cmake-build-debug  # run the unit tests
```

### Options

| Option | Default | Effect |
|---|---|---|
| `CPPCRYSTAL_BUILD_TESTS` | ON | unit test suite |
| `CPPCRYSTAL_BUILD_DEMO` | ON | `cppcrystal_demo` executable |
| `CPPCRYSTAL_BUILD_ORACLE_TESTS` | OFF | validate against reference spglib v2.7.0 |
| `CPPCRYSTAL_ENABLE_SANITIZERS` | OFF | Address + UndefinedBehavior sanitizers |
| `CPPCRYSTAL_BUILD_TOOLS` | OFF | offline data generators |

The `oracle` preset builds reference spglib via `FetchContent` and links it only
into the oracle tests — never into the library — to cross-check every result
(space groups, magnetic groups, reductions, k-point meshes) against it:

```bash
cmake --preset oracle && ctest --test-dir cmake-build-debug
```

### Use it in your project

As a subdirectory:

```cmake
add_subdirectory(CppCrystal)
target_link_libraries(your_target PRIVATE cppcrystal::cppcrystal)
```

Or installed:

```bash
cmake --install cmake-build-release --prefix /your/prefix
```

```cmake
find_package(CppCrystal REQUIRED)
target_link_libraries(your_target PRIVATE cppcrystal::cppcrystal)
```

Only `include/` ships, so an installed consumer reaches exactly what the
umbrella header reaches and no more. Eigen, Boost.Container and Boost.LEAF are
public dependencies, and because this project fetches them rather than taking
them from the system, `cmake --install` writes their CMake packages into the
same prefix — the result is self-contained and needs nothing installed
system-wide.

---

## Attribution

CppCrystal stands on two existing projects and is grateful to both.

- **[spglib](https://github.com/spglib/spglib)** (Atsushi Togo and contributors),
  BSD-3-Clause. The symmetry algorithms — space-group determination, lattice
  reduction, Wyckoff refinement, magnetic space-group machinery, and k-point
  reduction — and the symmetry databases are derived from spglib's C core. The
  reference implementation is also used (built separately, never linked into the
  library) as a validation oracle in the test suite. This library targets
  spglib **v2.7.0**.

- **[PyXtal](https://github.com/qzhu2017/PyXtal)** (Qiang Zhu, Scott Fredericks,
  and contributors), MIT. The structure-generation layer and the object-oriented
  framing — groups as standalone objects that own their Wyckoff positions, named
  factory functions, and random crystal/cluster generation by seating a
  composition on Wyckoff sites — follow PyXtal's design.

Both upstream licenses (BSD-3-Clause and MIT) are permissive and impose no
copyleft, so this combined and derived work is distributed under the
BSD-3-Clause license. The original spglib and PyXtal copyright notices are
retained for the portions derived from them — the symmetry core and the
`spacegroup_*` / `msg_*` / `sitesym_*` / `hall_generators*` data tables from
spglib, and the rod-group tables from PyXtal. The full verbatim notices are in
the [THIRD-PARTY NOTICES](LICENSE) section of the license file. The element
covalent-radius data is from Cordero et al. (2008), cited there.
