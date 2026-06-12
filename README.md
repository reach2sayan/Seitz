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
| Space-group determination | `get_dataset(cell)` | identity, operations, standardized cell, per-atom Wyckoff data |
| Layer / rod / point groups | `get_layer_dataset`, `RodGroup`, `PointGroup` | 2D-, 1D-, and 0D-periodic crystals |
| Stateful analysis | `analysis::SymmetryAnalyzer` | owns a cell, memoizes every derived query |
| Cell standardization | `standardize_cell`, `find_primitive`, `refine_cell` | primitive/conventional × idealized/as-input |
| Lattice reduction | `reduce::delaunay_reduce`, `reduce::niggli_reduce` | also 2D Delaunay for layers |
| Group catalogs | `group::SpaceGroup`, `PointGroup`, `RodGroup` | structure-free, queryable Wyckoff positions |
| Subgroup relations | `group::SubgroupGraph` | translationengleiche maximal-subgroup lattice of all 230 groups |
| Crystal generation | `generate::random_crystal`, `random_layer_crystal`, `random_cluster` | place a composition on Wyckoff sites at random |
| Magnetic structures | `get_magnetic_dataset(mcell)` | magnetic space-group identification, standardized magnetic cell |
| Reciprocal-space meshes | `kpoint::ir_reciprocal_mesh`, `stabilized_reciprocal_mesh` | irreducible k-points, Brillouin-zone folding |

A single umbrella header pulls in the whole API:

```cpp
#include <spglib/spglib.hpp>
```

---

## Quick start

### Determine a space group

```cpp
#include <spglib/spglib.hpp>

spglib::Cell rutile{lattice, positions, types}; // lattice columns = basis vectors
auto result = spglib::get_dataset(rutile);

spglib::leaf::try_handle_all(
    [&]() -> spglib::Result<void> {
      BOOST_LEAF_AUTO(ds, spglib::get_dataset(rutile));
      std::printf("space group %d (%s)\n",
                  ds.spacegroup_number,
                  std::string(ds.international_symbol).c_str());
      return {};
    },
    [](spglib::e_spacegroup_search_failed) {
      std::puts("no space group found at this tolerance");
    });
```

`Dataset` carries the international and Hall numbers, the operations of the input
cell, the transformation onto the conventional setting, the fully standardized
cell, and per-atom Wyckoff letters, site-symmetry symbols, and equivalence
classes.

### Reuse work with the analyzer

When you need several derived quantities from one cell, build a
`SymmetryAnalyzer`. It owns the cell and the tolerances and caches each stage of
the pipeline, so the dataset is computed once and every projection
(`operations()`, `wyckoffs()`, `primitive_cell()`, `standardized_cell()`, …) is
served from that cache:

```cpp
auto sa = spglib::analysis::SymmetryAnalyzer::from_cell(rutile);
auto number = sa.spacegroup_number();   // runs determination
auto prim   = sa.primitive_cell();      // reuses it
auto std    = sa.standardized_cell();   // reuses it
```

The analyzer is immutable after construction — to analyze at a different
tolerance, build a new one.

### Generate a random crystal

```cpp
auto sg   = spglib::group::SpaceGroup::from_number(225).value();   // Fm-3m
auto comp = spglib::Composition{{ {/*Z=*/11, 4}, {/*Z=*/17, 4} }}; // NaCl
auto xtal = spglib::generate::random_crystal(sg, comp);
```

`random_crystal` enumerates the valid ways to seat the composition on the group's
Wyckoff positions, then resamples free coordinates and the lattice until the
minimum-distance criterion is met. The same shape generates 2D layers
(`random_layer_crystal`), 1D rods (rod-group overload), and 0D clusters
(`random_cluster` over a `PointGroup`).

### Irreducible k-points

```cpp
auto mesh = spglib::kpoint::ir_reciprocal_mesh(
    cell, /*mesh=*/{8, 8, 8}, /*is_shift=*/{0, 0, 0},
    /*time_reversal=*/true);
// mesh.num_ir  -> count of irreducible points
// mesh.mapping -> each grid point's irreducible representative
```

---

## Design

**Errors are values.** Nothing throws and no function returns a magic sentinel.
A fallible call returns `Result<T>` (`= boost::leaf::result<T>`); failures carry
typed context — `e_spacegroup_search_failed`, `e_invalid_lattice{det}`,
`e_invalid_mesh`, `e_magnetic_symmetry_search_failed`, and so on. "Absent" is
`std::optional`, never `-1`.

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

**Thread-safety.** The shared global tables are primed by `spglib::warmup()` and
are race-free for concurrent reads afterward. Per-instance analyzer caches are
not; call `.warm()` on one thread before sharing an analyzer read-only.

### Source layout

```
include/spglib/         public headers
  core/                 Cell, symmetry operations, errors, tolerance, periodicity
  analysis/             SymmetryAnalyzer, MagneticSymmetryAnalyzer (memoizing facades)
  group/                SpaceGroup, PointGroup, RodGroup, Wyckoff, SubgroupGraph
  generate/             random_crystal / random_cluster + Wyckoff combinatorics
  symmetry/             symmetry search, point-group matching, primitive cell
  spacegroup/           Hall-symbol matching, space-group determination
  reduce/               Delaunay and Niggli lattice reduction
  refine/               standardization and Wyckoff refinement
  magnetic/  spin/      magnetic space-group identification and search
  kpoint/               grids, irreducible meshes, Brillouin-zone mapping
  data/                 constexpr symmetry / magnetic / element catalogs
  math/                 integer matrices, fractional coordinates
src/                    implementations, mirroring the header tree
tools/                  offline data-table generators (Python + one C++ tool)
tests/                  unit tests and oracle tests against reference spglib
docs/                   workflow flowchart and porting status
```

See [`docs/WORKFLOW.md`](docs/WORKFLOW.md) for a flowchart of the determination,
generation, magnetic, and k-point pipelines.

---

## Building

Requires a C++23 compiler (GCC 13+/Clang 17+), CMake ≥ 3.28, **Eigen 3.4**, and
**Boost 1.83** (headers only, for `boost::leaf`). **Catch2 3** is fetched for the
tests.

```bash
cmake --preset default              # configure (Debug)
cmake --build cmake-build-debug     # build the library + demo
ctest --test-dir cmake-build-debug  # run the unit tests
```

### Options

| Option | Default | Effect |
|---|---|---|
| `SPGLIB_BUILD_TESTS` | ON | unit test suite |
| `SPGLIB_BUILD_DEMO` | ON | `cppcrystal_demo` executable |
| `SPGLIB_BUILD_ORACLE_TESTS` | OFF | validate against reference spglib v2.7.0 |
| `SPGLIB_USE_MKL` | OFF | Intel MKL as Eigen's BLAS/LAPACK backend |
| `SPGLIB_ENABLE_SANITIZERS` | OFF | Address + UndefinedBehavior sanitizers |
| `SPGLIB_BUILD_TOOLS` | OFF | offline data generators |

The `oracle` preset builds reference spglib via `FetchContent` and links it only
into the oracle tests — never into the library — to cross-check every result
(space groups, magnetic groups, reductions, k-point meshes) against it:

```bash
cmake --preset oracle && ctest --test-dir cmake-build-debug
```

### Use it in your project

```cmake
add_subdirectory(CppCrystal)
target_link_libraries(your_target PRIVATE cppcrystal::cppcrystal)
```

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
spglib, and the rod-group tables (`rod_group_tables.hpp`) from PyXtal. The full
verbatim notices are in the [THIRD-PARTY NOTICES](LICENSE) section of the
license file. The element covalent-radius data is from Cordero et al. (2008),
cited there.
