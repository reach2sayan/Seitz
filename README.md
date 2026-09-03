# CppCrystal

A modern C++23 crystallography library for symmetry analysis and structure
generation. It determines the space group of a crystal, standardizes and refines
cells, classifies magnetic structures, reduces reciprocal-space meshes to their
irreducible wedge, and generates random crystals, layers, rods, and clusters
from a chosen symmetry group.

The library is value-semantic and exception-free: every fallible operation
returns a `boost::leaf::result<T>` carrying either the answer or a structured
error, and every result type is a plain owning value with no manual memory
management. Symmetry tables (space, layer, rod, point, and magnetic groups) are
compiled into the binary as `constexpr` catalogs and decoded on demand.

> **License & attribution.** CppCrystal is released under the BSD-3-Clause
> license (see [`LICENSE`](LICENSE)). Its algorithms are derived from
> [**spglib**](https://github.com/spglib/spglib) (BSD-3-Clause) and its
> object-oriented, generation-focused API is modeled on
> [**PyXtal**](https://github.com/qzhu2017/PyXtal) (MIT). Full credits are in
> [Attribution](#attribution).

---

## What it does

| Capability | Entry point | Notes |
|---|---|---|
| Space-group determination | `analysis::SymmetryAnalyzer::from_cell(cell)` | owns the cell, memoizes every derived query |
| Layer groups | the same analyzer | a `Cell` with one aperiodic axis takes the layer path |
| Cell standardization | `.standardized_cell<Setting, Idealize>()` | primitive/conventional × idealized/as-input |
| Magnetic structures | `analysis::MagneticSymmetryAnalyzer::from_cell(mcell)` | magnetic space-group identification, standardized magnetic cell |
| Lattice reduction | `Lattice::niggli`, `::delaunay`, `::delaunay_in_plane` | on the lattice itself; the in-plane form for layers and monoclinic cells |
| Group catalogs | `group::SpaceGroup::of`, `PointGroup`, `RodGroup` | structure-free, queryable Wyckoff positions |
| Subgroup relations | `group::SubgroupGraph` | translationengleiche maximal-subgroup lattice of all 230 groups |
| Structure generation | `generate::Generator<G>` | one generator over space, layer, rod and point groups |
| Reciprocal-space meshes | `kpoint::Mesh`, `ReciprocalMesh`, `BrillouinZone` | irreducible k-points, Brillouin-zone folding |

A single umbrella header pulls in the whole API:

```cpp
#include <cppcrystal/cppcrystal.hpp>
```

---

## Quick start

### Determine a space group

There is one entry point, and it is an object. `SymmetryAnalyzer` owns the cell
and the tolerances and memoizes each stage of the pipeline, so determination runs
once and every projection is served from that cache:

```cpp
#include <cppcrystal/cppcrystal.hpp>
using namespace cppcrystal;

Cell const rutile{Lattice{basis}, positions, types}; // lattice columns = basis vectors
auto const sa = analysis::SymmetryAnalyzer::from_cell(rutile);

leaf::try_handle_all(
    [&]() -> Result<void> {
      BOOST_LEAF_AUTO(hall, sa.hall());          // runs determination
      BOOST_LEAF_AUTO(ops, sa.operations());     // reuses it
      BOOST_LEAF_AUTO(prim, sa.primitive_cell()); // reuses it
      std::printf("space group %d (%s), %zu operations\n",
                  data::spacegroup_type(hall).number,
                  data::spacegroup_type(hall).international_short.data(),
                  ops.size());
      return {};
    },
    [](e_spacegroup_search_failed) {
      std::puts("no space group found at this tolerance");
    });
```

The analyzer is immutable after construction — to analyze at a different
tolerance, build a new one. Its projections hand back references into the memo,
so it must outlive them: the rvalue overloads are deleted, which makes
`from_cell(cell).hall()` a compile error rather than a dangling read. A layer
cell (one aperiodic axis) goes through the same object; the periodicity picks the
path.

`Dataset` carries the Hall key, the Bravais lattice, the conventional `Setting`,
the operations of the input cell, one `Site` per atom (Wyckoff letter,
site-symmetry symbol, equivalence class) and the fully standardized `Cell`.

### Generate a random crystal

```cpp
BOOST_LEAF_AUTO(fm3m, group::SpaceGroup::from_number(GroupFamily::space, 225));
generate::Composition const nacl{{11, 4}, {17, 4}};   // Na4 Cl4

BOOST_LEAF_AUTO(xtal, generate::Generator(*fm3m, {.seed = 42})(nacl));
// xtal.cell        -> the generated Cell
// xtal.assignment  -> which Wyckoff position each atom was seated on
```

`Generator` enumerates the valid ways to seat the composition on the group's
Wyckoff positions, then resamples free coordinates and the lattice until the
minimum-distance criterion is met. The same generator covers every family:
substitute a layer `SpaceGroup`, a `RodGroup` or a `PointGroup` and the traits
supply that family's periodicity, seed window and random metric. The result is
always a `Cell` that describes itself — a cluster is simply a cell with no
periodic axis.

### Irreducible k-points

```cpp
auto const mesh = *kpoint::Mesh::of({8, 8, 8});   // rejects a non-positive mesh
auto rmesh = sa.reciprocal_mesh(mesh, TimeReversal::yes);
// rmesh->num_irreducible()  -> count of irreducible points
// rmesh->mapping()          -> each grid point's irreducible representative
// rmesh->brillouin_zone(lattice) -> Brillouin-zone folding
```

---

## Design

**Errors are values.** Nothing throws and no function returns a magic sentinel.
A fallible call returns `Result<T>` (`= boost::leaf::result<T>`); failures carry
typed context — `e_spacegroup_search_failed`, `e_invalid_lattice{det}`,
`e_magnetic_symmetry_search_failed`, and so on. "Absent" is `std::optional`,
never `-1`.

**Invalid states are unrepresentable, so most errors never arise.** A
`HallNumber` carries its family and a validated index (there is no negative-Hall
convention and no sentinel row in any catalog); a `UniNumber` and a
`kpoint::Mesh` are likewise built through an `of()` that is the only way in and
the only place the check lives. A function taking one of these keys is total.

**Cells describe periodicity per axis.** A `Cell` holds a column-major lattice,
fractional positions, integer types, and a `CellPeriodicity`. The same type
represents a 3D crystal, a 2D layer (one aperiodic axis), a 1D rod (two), and a
0D cluster (three) — the periodicity drives which symmetry path runs.

**Groups are first-class objects.** `SpaceGroup`, `PointGroup`, and `RodGroup`
are built from a group number and own their operations and Wyckoff positions as
queryable objects, not parallel arrays of letters and multiplicities. This is the
catalog face of the symmetry database and the foundation for generation and
subgroup work, independent of any atomic structure.

**Symmetry tables are compiled in.** Each database (space, magnetic, site
symmetry, rod, elements) is transcribed offline into a generated header by the
scripts in [`tools/`](tools), then decoded once at compile time into a
`constexpr` catalog indexed directly by its key (Hall # 1..530, UNI # 1..1651).
Operation matrices, which cannot be `constexpr` (they are `Eigen` types), are
decoded lazily at runtime and cached. Generated tables are split into a public
metadata header and a private operation header so the encoded packing never
reaches a public interface.

**The family is a compile-time parameter.** Space and layer groups run the same
pipeline, instantiated as `if constexpr` on a `GroupFamily` template argument
rather than branching on the sign of a Hall number at each stage; a single
`dispatch_family()` is the one runtime branch, at the top.

**Thread-safety.** Everything is safe to share. `group::SpaceGroup::of` is a
Boost.Flyweight — one immutable object per Hall setting, built on first use and
shared thereafter; the generated tables are read-only after their one-time
decode; and every analyzer memo is race-free, so a const analyzer can be shared
across threads from the moment it is built (the first caller of each query pays
for it). `cppcrystal::warmup()` builds every setting up front (both families, in
about 30 ms) to move that cost off the query path.

### Source layout

The public interface and the implementation are separated by directory, not by
convention: `#include <cppcrystal/…>` is the API, `#include "module/x.hpp"` is
internal, and only `include/` is installed.

```
include/cppcrystal/       the public API — everything the umbrella header reaches
  core/                 Cell, Lattice, keys, operations, errors, tolerance, periodicity
  analysis/             SymmetryAnalyzer, MagneticSymmetryAnalyzer (memoizing facades)
  group/                SpaceGroup, PointGroup, RodGroup, Wyckoff, SubgroupGraph
  generate/             Generator<G>, Wyckoff assignments, distance checking
  kpoint/               Mesh, ReciprocalMesh, BrillouinZone
  data/                 constexpr symmetry / magnetic / element catalogs
src/                      the implementation, headers beside their .cpp files
  symmetry/             symmetry search, point-group matching, primitive cell
  spacegroup/           Hall-symbol matching, space-group determination
  reduce/               Delaunay and Niggli lattice reduction
  refine/               standardization and Wyckoff refinement
  magnetic/  spin/      magnetic space-group identification and search
  data/                 the generated operation tables and their decoders
  math/                 integer matrices, fractional coordinates
tools/                    offline data-table generators (Python + one C++ tool)
tests/                    unit tests and oracle tests against reference spglib
docs/                     workflow flowchart, mathematical background
```

See [`docs/WORKFLOW.md`](docs/WORKFLOW.md) for a flowchart of the determination,
generation, magnetic, and k-point pipelines, and
[`docs/MATHEMATICS.md`](docs/MATHEMATICS.md) for the group theory behind them
(orbits, stabilisers, Wyckoff positions, space-group identification).

---

## Building

Requires a C++23 compiler (the presets use **g++-15**) and CMake ≥ 3.28.
Everything else is fetched: **Eigen 5.0.0** (the first release whose fixed-size
matrices are literal types, so they can be built in `constexpr` context),
**Boost 1.88** (`container`, `leaf`, `parser`, `preprocessor`), and **Catch2 3**
for the tests. Nothing needs to be installed on the system.

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
| `CPPCRYSTAL_USE_MKL` | OFF | Intel MKL as Eigen's BLAS/LAPACK backend |
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

Or installed, with only `include/` shipping — the pipeline headers under `src/`
are on the library's `PRIVATE` include path and are not part of the exported
interface, so an installed consumer reaches exactly what the umbrella header
reaches and no more:

```bash
cmake --install cmake-build-release --prefix /your/prefix
```

```cmake
find_package(CppCrystal REQUIRED)
target_link_libraries(your_target PRIVATE cppcrystal::cppcrystal)
```

Eigen and Boost are public dependencies of the exported target, and this project
fetches both rather than taking them from the system, so `cmake --install` writes
their CMake packages into the same prefix. The result is self-contained: a
consumer needs nothing installed system-wide.

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
spglib, and the rod-group tables (`src/data/rod_group_tables.hpp`) from PyXtal. The full
verbatim notices are in the [THIRD-PARTY NOTICES](LICENSE) section of the
license file. The element covalent-radius data is from Cordero et al. (2008),
cited there.
