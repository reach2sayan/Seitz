# Seitz

**Crystallographic symmetry: determination, classification, and construction, in C++23.**

[![C++23](https://img.shields.io/badge/C%2B%2B-23-blue.svg)](https://en.cppreference.com/w/cpp/23)
[![CMake](https://img.shields.io/badge/CMake-3.28%2B-064F8C.svg?logo=cmake&logoColor=white)](https://cmake.org)
[![Compilers](https://img.shields.io/badge/compilers-GCC%2015%2B%20%7C%20Clang%20%7C%20MSVC-brightgreen.svg)](#building)
[![License](https://img.shields.io/badge/license-BSD--3--Clause-green.svg)](LICENSE)
[![Validated against spglib](https://img.shields.io/badge/validated%20against-spglib%20v2.7.0-orange.svg)](https://github.com/spglib/spglib)

Seitz computes the symmetry group of a periodic structure and everything that
follows from it: the standardized cell, the Wyckoff decomposition of the atomic
orbits, the magnetic (Shubnikov) group of a spin arrangement, the irreducible
wedge of a reciprocal-space mesh, the maximal-subgroup graph of the 230 space
groups, random structures carrying a prescribed group, and the
symmetry-distinct cluster orbits of a substitutional alloy.

Everything ships behind one header:

```cpp
#include <seitz/seitz.hpp>
```

---

## Contents

| | |
|---|---|
| [Objects](#objects) | lattices, cells, Seitz operators, operation sets |
| [Determination](#determination) | structure $\to$ space group |
| [Standardization and reduction](#standardization-and-reduction) | Niggli, Delaunay, idealization |
| [Groups and Wyckoff positions](#groups-and-wyckoff-positions) | the catalogs as objects |
| [The subgroup graph](#the-subgroup-graph) | maximal $t$- and $k$-subgroups |
| [Magnetic symmetry](#magnetic-symmetry) | Shubnikov groups, BNS/UNI |
| [Reciprocal space](#reciprocal-space) | irreducible meshes, Brillouin zone |
| [Structure generation](#structure-generation) | random crystals on a group |
| [Alloy configurations](#alloy-configurations) | cluster orbits, CVM |
| [Error model and invariants](#error-model-and-invariants) | `Result<T>`, total functions |
| [Building](#building) · [Python](#python) · [License and attribution](#license-and-attribution) | |

---

## Objects

A **lattice** is a rank-3 $\mathbb{Z}$-module $\Lambda = A\mathbb{Z}^3$, stored
as the matrix $A \in GL(3,\mathbb{R})$ whose *columns* are the basis vectors
$\mathbf{a}, \mathbf{b}, \mathbf{c}$. All geometry follows from the Gram (metric)
tensor $G = A^{\mathsf{T}}A$, and the cell volume is $|\det A|$.

```cpp
Lattice const lat{basis};                      // columns are a, b, c
lat.metric();                                  // G = AᵀA
lat.volume();                                  // |det A|
lat.to_cartesian(x);  lat.to_fractional(r);    // A x  and  A⁻¹ r
lat.transformed(P);                            // A P — change of basis
Lattice::from_basis(A);                        // nullopt if |det A| < kZeroPrec
```

A **cell** is $(A, X, \tau, \pi)$: a lattice, fractional coordinates
$X \in [0,1)^{N \times 3}$, integer types $\tau$, and a per-axis periodicity
$\pi \in \{\text{periodic}, \text{aperiodic}\}^3$. The periodicity is what
selects the symmetry path, so one type carries all four dimensionalities — a 3D
crystal, a 2D layer (one aperiodic axis, layer groups), a 1D rod (two), and a 0D
cluster (three) — and there is no separate entry point per case.

```cpp
Cell const c{lat, positions, types};                       // all_periodic() by default
Cell const slab = c.with_periodicity(aperiodic_along(2));  // layer-group path
c.transformed(P);  c.supercell({2, 2, 1});  c.translated(s);
```

A **symmetry operation** is a Seitz operator $\{R \mid t\}$ acting on
fractional coordinates,

$$\{R \mid t\}\,x = Rx + t, \qquad R \in GL(3,\mathbb{Z}),\ t \in \mathbb{R}^3/\mathbb{Z}^3 ,$$

with composition $\{R_1|t_1\}\{R_2|t_2\} = \{R_1R_2 \mid R_1t_2 + t_1\}$ and
inverse $\{R^{-1} \mid -R^{-1}t\}$, defined exactly when $R$ is unimodular.
Under a change of basis $P$ an operator conjugates as
$\{R|t\} \mapsto \{PRP^{-1} \mid Pt\}$.

```cpp
SymmetryOperation const op{.rotation = R, .translation = t};
op.apply(x);                       // Rx + t
auto const g = a * b;              // (a * b).apply(x) == a.apply(b.apply(x))
op.inverse();                      // optional — nullopt if R is not unimodular
conjugated_by(op, P, P_inv);       // P R P⁻¹, P t
same_operation(a, b, symprec);     // equality modulo Z³, within tolerance
```

An **`OperationSet<Op>`** is an immutable, range-modelling set of these — the
group itself. It is the object the determination hands back, and it can be
interrogated on its own:

```cpp
ops.rotations();                   // the point group, as integer matrices
ops.pure_translations();           // the centring vectors of the coset {I | t}
ops.conjugated_by(P, P_inv);       // the group in another basis
BOOST_LEAF_AUTO(m, ops.spacegroup(lattice, LatticeSetting::conventional));
// m.hall, m.bravais, m.origin_shift, m.type()
```

---

## Determination

Given a cell and a tolerance $\varepsilon$, determination recovers the space
group $\mathcal{G}$ of the structure. The pipeline is:

1. **Translation subgroup.** Find the candidate pure translations
   $t$ with $\tau(x_i + t) = \tau(x_i)$ for all $i$ up to $\varepsilon$; they
   generate $\Lambda_{\text{prim}} \supseteq \Lambda$, and the primitive cell is
   its reduced basis.
2. **Lattice point group.** Delaunay-reduce the primitive metric $G$ and collect
   the integer matrices $R$ with $R^{\mathsf{T}} G R = G$ — the holohedry, a
   finite group of order $\le 48$.
3. **Space-group operations.** For each such $R$, solve for a translation part
   $t$ that maps the decorated structure onto itself; the pairs that survive are
   the operations of the input cell.
4. **Matching.** Identify the resulting group against the Hall-symbol database
   by conjugating it into each candidate setting, which yields the Hall number,
   the Bravais lattice, and the origin shift.
5. **Refinement.** Transport the atoms into the standardized setting, assign
   each to a Wyckoff position by its stabilizer, and record the orbits.

`SymmetryAnalyzer` is the interface. It owns the cell and the tolerances,
memoizes each stage, and is immutable after construction — so determination runs
once and every subsequent query is served from the memo. A `const` analyzer is
safe to share across threads from the moment it exists.

```cpp
auto const sa = analysis::SymmetryAnalyzer::from_cell(cell, Tolerance{},
                                                     /*setting=*/std::nullopt);

leaf::try_handle_all(
    [&]() -> Result<void> {
      BOOST_LEAF_AUTO(hall, sa.hall());          // runs the pipeline
      BOOST_LEAF_AUTO(ops,  sa.operations());    // reuses it
      BOOST_LEAF_AUTO(prim, sa.primitive_cell());
      auto const &type = data::spacegroup_type(hall);
      std::println("{} ({}), |G| = {}, Z' cell of {} atoms",
                   type.international_short, type.number, ops.size(), prim.size());
      return {};
    },
    [](e_spacegroup_search_failed) { std::println("no group at this tolerance"); });
```

| Query | Yields |
|---|---|
| `sa.hall()` | `HallNumber` — the key everything else is indexed on |
| `sa.spacegroup_type()` | number, international/Schoenflies symbols, Hall symbol, point group |
| `sa.operations()` | $\{R\mid t\}$ of the **input** cell |
| `sa.cell_operations()` | operations of the input cell as given, without standardization |
| `sa.lattice_symmetry()` | the lattice point group $\{R : R^{\mathsf{T}}GR = G\}$ |
| `sa.sites()` | one `Site` per atom: Wyckoff letter, site-symmetry symbol, orbit representative |
| `sa.primitive_cell()` | the primitive cell found in step 1 |
| `sa.standardized_cell()` | conventional, idealized |
| `sa.standardized_cell<CellSetting::primitive, Idealize::no>()` | the other three settings, resolved at compile time |
| `sa.determine()` | the whole `Dataset` at once |

The `Dataset` also carries the `Setting` that relates the input to the
standard description — the change of basis $P$, the origin shift $s$, and the
rigid rotation onto the idealized lattice — so any tabulated quantity can be
pulled back onto the user's own coordinates.

Because the analyzer's queries return references into its memo, it must outlive
them: the rvalue overloads are deleted, which makes `from_cell(cell).hall()` a
compile error rather than a dangling read.

---

## Standardization and reduction

Two classical reductions of a lattice basis, both returning a new `Lattice`:

- **Niggli** (Křivý–Gruber): the unique reduced basis of a lattice, obtained by
  normalizing the six metric parameters $(a^2, b^2, c^2, bc, ac, ab)$ under the
  reduction conditions.
- **Delaunay**: reduction of the extended basis
  $\{\mathbf{a}, \mathbf{b}, \mathbf{c}, -(\mathbf{a}+\mathbf{b}+\mathbf{c})\}$
  until all pairwise scalar products are non-positive. This is the reduction the
  holohedry search runs on, since it exposes the lattice's full point symmetry.

```cpp
BOOST_LEAF_AUTO(n, lat.niggli(eps));
BOOST_LEAF_AUTO(d, lat.delaunay(symprec));
BOOST_LEAF_AUTO(p, lat.delaunay_in_plane(unique_axis, symprec));  // layer case
lat.rigid_rotation_to(ideal);   // R with ideal = R · lat
```

Idealization is the second half: a standardized cell may either keep the input's
own (possibly strained) metric or be replaced by the exact metric its Bravais
class demands. The two axes — `CellSetting` $\times$ `Idealize` — are template
parameters of `standardized_cell`, so the four combinations are distinct
instantiations rather than runtime flags.

---

## Groups and Wyckoff positions

The catalogs are objects, not parallel arrays. `SpaceGroup`, `PointGroup` and
`RodGroup` share a face — number, symbol, order, operations, Wyckoff positions —
without a runtime hierarchy.

```cpp
BOOST_LEAF_AUTO(g, group::SpaceGroup::from_number(GroupFamily::space, 225));
group::SpaceGroup const &h = group::SpaceGroup::of(hall);   // by Hall setting
g->order();  g->symbol();  g->centering();  g->operations();
```

A **Wyckoff position** is an orbit type: the set of points whose stabilizer
(site-symmetry group) is conjugate to a given subgroup $S \le \mathcal{G}$. Its
locus is an affine subspace $x(\lambda) = x_0 + \sum_i \lambda_i b_i$ of
dimension $d$ (the degrees of freedom), and the orbit–stabilizer theorem fixes
the multiplicity:

$$m \cdot |S| = |\mathcal{G}| .$$

```cpp
for (group::Wyckoff const &w : g->wyckoffs()) {
  w.multiplicity();  w.letter();  w.site_symmetry();  w.degrees_of_freedom();
  w.operations();                 // the stabilizer S; |S| · m = |G|
  w.sample(params);               // x₀ + Σ λᵢ bᵢ
  w.canonical(xyz);               // idempotent projection onto the locus
  w.orbit(xyz, periodicity);      // the full orbit, one row per image
}
BOOST_LEAF_AUTO(wa, g->wyckoff('a'));
```

`canonical` is the projector onto the locus, so an approximate coordinate still
generates the exact orbit — the numerically robust way to seat an atom on a
special position.

---

## The subgroup graph

Maximal subgroups of the 230 space-group types form a directed graph. An edge
$\mathcal{G} \to \mathcal{H}$ records that $\mathcal{H}$ is a **maximal**
subgroup of $\mathcal{G}$, of one of two kinds:

- **translationengleiche** ($t$): same translation lattice, smaller point group,
  $[\mathcal{G}:\mathcal{H}] = |P_{\mathcal{G}}| / |P_{\mathcal{H}}|$;
- **klassengleiche** ($k$): same point group, sublattice of index
  $[\Lambda_{\mathcal{G}} : \Lambda_{\mathcal{H}}]$.

Each edge carries the transformation taking the supergroup's conventional
description to the subgroup's — a basis matrix $P$ and an origin shift $s$:

$$(\mathbf{a}_H\ \mathbf{b}_H\ \mathbf{c}_H) = (\mathbf{a}_G\ \mathbf{b}_G\ \mathbf{c}_G)\,P, \qquad x_G = P x_H + s ,$$

with $s$ defined modulo $P\mathbb{Z}^3$, so an operation of $\mathcal{G}$ can be
pushed into $\mathcal{H}$'s frame as
$\{P^{-1}RP \mid P^{-1}(Rs + t - s)\}$ — defined exactly when $P^{-1}RP$ is
integral, i.e. when $R$ preserves $\Lambda_H$.

```cpp
using group::SubgroupGraph;
for (int id : SubgroupGraph::maximal_subgroups(225, SubgroupKind::translationengleiche)) {
  group::SubgroupEdge const e = SubgroupGraph::edge(id);
  e.super; e.sub; e.index; e.kind; e.hall; e.basis; e.origin;
  e.in_subgroup_frame(op);            // optional<SymmetryOperation>
}
SubgroupGraph::minimal_supergroups(number);      // the reverse relation
SubgroupGraph::is_subgroup(sub, super);          // reflexive reachability
SubgroupGraph::path(super, sub);                 // shortest chain of edges
```

The relation table is `consteval`, and so is its transitive closure: reachability
is a bit test in a compile-time Warshall closure, not a search. There is one
storage — the graph is a Boost.Graph model (`bidirectional_graph`,
`vertex_list_graph`, `edge_list_graph`) over that same table, so
`breadth_first_search`, `filtered_graph`, `reverse_graph` and the BGL visitors
apply to it unchanged, with no second representation to keep in sync.

---

## Magnetic symmetry

A magnetic (Shubnikov) group extends the Seitz operator with the time-reversal
generator $\theta$: elements are $\{R \mid t\}$ and $\{R \mid t\}\theta$, the
latter reversing axial site tensors. Classification follows the
Bärnighausen/BNS construction types:

| Type | Structure |
|---|---|
| I | $\mathcal{M} = \mathcal{G}$, no anti-operation |
| II | $\mathcal{M} = \mathcal{G} + \mathcal{G}\theta$ (grey) |
| III | $\mathcal{M} = \mathcal{H} + (\mathcal{G} \setminus \mathcal{H})\theta$, $\mathcal{H}$ a $t$-subgroup of index 2 |
| IV | as III with $\mathcal{H}$ a $k$-subgroup: $\theta$ paired with an anti-translation |

```cpp
auto const ma = analysis::MagneticSymmetryAnalyzer::from_cell(mcell, MagneticTolerance{});
BOOST_LEAF_AUTO(uni,  ma.uni());               // UNI number (1651 types)
BOOST_LEAF_AUTO(ops,  ma.operations());        // MagneticOperations
BOOST_LEAF_AUTO(std,  ma.standardized_cell()); // tensors rotated into the standard basis
ma.magnetic_spacegroup_type();                 // BNS/OG symbols, type I–IV
```

Collinear and non-collinear moments are both handled; the moment tolerance is
separate from the positional one, so spin and geometry are compared on their own
scales.

---

## Reciprocal space

A sampling mesh is a quotient of the reciprocal lattice. Points are addressed on
the doubled grid so that shifted (Monkhorst–Pack) meshes stay integral:

$$q = \frac{2a + \sigma}{2d}, \qquad a \in \mathbb{Z}^3,\ \sigma \in \{0,1\}^3,\ d \in \mathbb{Z}_{>0}^3 .$$

`Mesh` is pure grid arithmetic — the bijection between addresses and linear
indices, `constexpr` throughout, with its round-trip identities asserted at
compile time.

The reduction acts with the **reciprocal** point group: real-space rotations
transposed, together with the inversion partner when time reversal is imposed.
Each grid point is mapped to the smallest index in its orbit, so
`mapping()[i] == i` characterizes the irreducible representatives.

```cpp
auto const mesh = *kpoint::Mesh::of({8, 8, 8});     // nullopt on a non-positive mesh
BOOST_LEAF_AUTO(rm, sa.reciprocal_mesh(mesh, TimeReversal::on));
rm.num_irreducible();                 // |IBZ|
rm.mapping();                         // grid point → representative
rm.images_of(address);                // the orbit of one address
auto const bz = rm.brillouin_zone(reciprocal_lattice);
```

`BrillouinZone` relocates every point to the reciprocal-lattice image nearest
the origin, duplicating boundary points that are equidistant from several images
within tolerance — the convention that makes zone-boundary weights come out
right.

---

## Structure generation

The inverse problem: given a group $\mathcal{G}$ and a composition
$\{(Z_k, n_k)\}$, produce a structure whose symmetry is $\mathcal{G}$. An
assignment is a multiset of Wyckoff positions per species satisfying

$$\sum_{i \in \text{assignment}(k)} m_i = n_k ,$$

with positions of zero degrees of freedom usable at most once (two atoms cannot
share a fixed point). The generator enumerates valid assignments — pruned by a
reachability table over the remaining counts — then, for each, samples the free
parameters $\lambda$ and a random lattice of the correct Bravais class until the
minimum-distance criterion is met.

```cpp
generate::Composition const comp{{11, 4}, {17, 4}};          // Na₄Cl₄
generate::Generator const gen{*g, {.seed = 42, .attempts_per_combination = 50}};

gen.compatible(comp);                 // is any assignment possible at all?
gen.assignments(comp);                // enumerate them
BOOST_LEAF_AUTO(x, gen(comp));        // x.cell, x.assignment (with generating coordinates)
```

`GenerateOptions` fixes the search: `seed` (fully deterministic), `scale` on the
element-aware size estimate, `distance` acceptance, `placement`
(`general_only` forbids special positions, so no accidental symmetry is added),
a caller-supplied `lattice` whose metric must be invariant under $\mathcal{G}$,
and `sites` — atoms pinned to a Wyckoff letter, and optionally to a coordinate
on it, before the search begins.

The generator is templated on the group family. `SpaceGroup` (3D and, when its
Hall key names the layer family, 2D), `RodGroup` (1D) and `PointGroup` (0D) each
supply their periodicity, seed window and random metric through
`GroupTraits<G>`; assignment enumeration, the distance test and the search are
written once. The result is always a `Cell` that describes its own periodicity —
a cluster is a cell with no periodic axis.

---

## Alloy configurations

For substitutional disorder on a fixed parent lattice, the observables are
symmetry orbits of decorated clusters. A cluster expansion writes a
configurational property as

$$E(\sigma) = \sum_\alpha m_\alpha J_\alpha \, \xi_\alpha(\sigma) ,$$

with $\xi_\alpha$ the orbit-averaged correlation functions of the point-function
basis and $m_\alpha$ the orbit multiplicities. Seitz enumerates the orbits and
builds the Cluster Variation Method ingredients for a chosen set of maximal
clusters — the configurations $j$ on each subcluster $c$ with multiplicities
$m_{jc}$, the matrix $V_c$ mapping correlations to cluster probabilities
$\rho_c = V_c\xi$, and the Kikuchi–Barker coefficients $k_c$ of the entropy

$$S = -k_B \sum_c k_c \sum_j m_{jc}\, \rho_{jc} \ln \rho_{jc} .$$

```cpp
BOOST_LEAF_AUTO(pool, alloy::ClustersPool::generate(parent, {.radii = {{2, 6.0}, {3, 4.5}}}));
for (alloy::Orbit const &o : pool) { o.multiplicity(); }

BOOST_LEAF_AUTO(cvm, alloy::Cvm::create(parent, maximal_clusters, basis));
cvm.clusters();    // per subcluster: configurations, V-matrix, Kikuchi–Barker coefficient
cvm.functions();   // the basis functions indexing ξ
```

This layer generates the ingredients; fitting the $J_\alpha$ and minimizing the
free energy are left to the caller.

---

## Error model and invariants

**Errors are values.** Nothing throws and nothing returns a magic sentinel. A
fallible call returns `Result<T>` (`= boost::leaf::result<T>`), and failures
carry typed context: `e_spacegroup_search_failed`, `e_invalid_lattice{det}`,
`e_atoms_too_close{distance}`, `e_incompatible_lattice`,
`e_invalid_transformation`, `e_niggli_failed`, and so on. "Absent" is
`std::optional`, never `-1`.

**Invalid states are unrepresentable, so most errors never arise.** A
`HallNumber`, a `UniNumber` and a `kpoint::Mesh` are each built through an
`of()` that is the only way in and the only place the check lives. A function
taking such a key is total: it cannot fail on a bad key, because no bad key can
be held.

**Value semantics throughout**, no manual memory management, no ownership to
reason about. **Everything is safe to share**: a `const` analyzer is usable from
any number of threads the moment it is built, the catalogs are immutable, and no
query mutates observable state. `seitz::warmup()` (about 30 ms) builds every
group setting up front to move first-use cost off the query path; it is an
optimization, never a precondition.

Every result is cross-checked against reference spglib v2.7.0 by an oracle test
suite that builds the reference implementation separately and never links it
into the library.

---

## Building

One line, and the presets carry the flags — CI runs exactly these:

```bash
cmake --preset release && cmake --build --preset release && ctest --preset release
```

| Preset | |
|---|---|
| `debug` | Debug, unit tests + spglib oracle |
| `release` | Release, unit tests + spglib oracle |
| `asan` | Debug, oracle, Address + UndefinedBehavior sanitizers |
| `python` | Release, the Python extension module |

Needs GCC 15+, Clang with libstdc++, or MSVC 19.43+; CMake $\ge$ 3.28; and a
Python 3.9+ interpreter, used at configure time to transcribe the symmetry
tables. Everything else is fetched rather than found — Eigen 5.0.0, Boost 1.88,
Catch2 3 — so nothing needs to be installed first. The presets default to
`gcc-15`/`g++-15`; pass `-DCMAKE_CXX_COMPILER=` for anything else. Off the
presets — `pip install .`, an IDE, a bare `cmake` — set `CXX`, since a default
`c++` is GCC 13 on Ubuntu 24.04 and too old for this tree:

```bash
CXX=g++-15 pip install .
```

| Option | Default | Effect |
|---|---|---|
| `SEITZ_BUILD_TESTS` | ON | unit test suite |
| `SEITZ_BUILD_DEMO` | ON | `seitz_demo` executable |
| `SEITZ_BUILD_ORACLE_TESTS` | OFF | cross-check every result against reference spglib v2.7.0 (all presets turn this ON) |
| `SEITZ_ENABLE_SANITIZERS` | OFF | Address + UndefinedBehavior sanitizers |
| `SEITZ_BUILD_PYTHON` | OFF | the Python extension module |

### Consuming it

```cmake
add_subdirectory(Seitz)                 # or: find_package(Seitz REQUIRED)
target_link_libraries(your_target PRIVATE seitz::seitz)
```

`cmake --install build/release --prefix /your/prefix` writes the vendored Eigen
and Boost packages into the same prefix, so an installed tree is self-contained
and `find_package(Seitz REQUIRED)` needs only `CMAKE_PREFIX_PATH`. Prebuilt
prefixes for Linux and Windows are attached to each GitHub release.

Only `include/` ships, so a consumer reaches exactly what the umbrella header
reaches. `seitz` is a shared library, `0.1.0` with `SONAME 0.1`; while the API is
pre-1.0 every minor version may break ABI. On Windows, and inside the Python
extension, it is a static archive instead.

---

## Python

The same API, NumPy-first: the C++ library does the work, and the Python layer
adds validated inputs, serializable records via **pydantic**, and an exception
hierarchy in place of `Result<T>`.

Wheels for CPython 3.11-3.14 (Linux x86-64, Windows x64) are attached to each
GitHub release, and need no compiler or GCC 15 on the target machine:

```bash
uv pip install https://github.com/reach2sayan/Seitz/releases/latest/download/seitz-<version>-cp313-cp313-manylinux_2_28_x86_64.whl
```

```python
import numpy as np
import seitz as sz

cell = sz.Cell(
    sz.Lattice(3.0 * np.eye(3)),          # columns are the basis vectors
    [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    [26, 26],
)

analyzer = sz.analyze(cell)                            # nothing computed yet
print(analyzer.spacegroup_type.international_short)    # Im-3m
print(len(analyzer.operations))                        # 96

group = sz.SpaceGroup.of(analyzer.hall)
print(group.wyckoffs[-1].multiplicity)                 # the general position

for e in sz.subgroups.maximal_subgroups(225, sz.SubgroupKind.translationengleiche):
    print(e.sub, e.index)

print(sz.DatasetRecord.from_analyzer(analyzer).model_dump_json(indent=2))
```

The analyzer memoizes, so it is the object you keep rather than a call you
repeat, and every query on it is thread-safe. "Absent" is `None`; failures raise
a `seitz.errors.SeitzError` subclass carrying its context —
`InvalidLatticeError.determinant`, `AtomsTooCloseError.distance`. Layer groups
are not a separate entry point: a cell built with
`periodicity=sz.aperiodic_along(2)` goes through the same analyzer.

From a checkout, the `python` preset builds the extension and runs the pytest
suite against the build tree:

```bash
cmake --preset python && cmake --build --preset python && ctest --preset python
```

---

## License and attribution

Seitz is released under the **BSD-3-Clause** license; see [`LICENSE`](LICENSE).

It stands on two existing projects and is grateful to both.

- **[spglib](https://github.com/spglib/spglib)** (Atsushi Togo and contributors),
  BSD-3-Clause. The symmetry algorithms — space-group determination, lattice
  reduction, Wyckoff refinement, magnetic space-group machinery, and k-point
  reduction — and the symmetry databases are derived from spglib's C core. The
  reference implementation is also used, built separately and never linked into
  the library, as a validation oracle in the test suite. This library targets
  spglib **v2.7.0**.

- **[PyXtal](https://github.com/qzhu2017/PyXtal)** (Qiang Zhu, Scott Fredericks,
  and contributors), MIT. The structure-generation layer and the object-oriented
  framing — groups as standalone objects owning their Wyckoff positions, named
  factory functions, and random crystal/cluster generation by seating a
  composition on Wyckoff sites — follow PyXtal's design.

The alloy layer follows the cluster-expansion and CVM formulations of ATAT
(A. van de Walle) and clusterX.

Both upstream licenses (BSD-3-Clause and MIT) are permissive and impose no
copyleft, so this combined and derived work is distributed under the
BSD-3-Clause license. The original spglib and PyXtal copyright notices are
retained for the portions derived from them — the symmetry core and the
`spacegroup_*` / `msg_*` / `sitesym_*` / `hall_generators*` data tables from
spglib, and the rod-group tables from PyXtal. The full verbatim notices are in
the [THIRD-PARTY NOTICES](LICENSE) section of the license file. The element
covalent-radius data is from Cordero et al. (2008), cited there.
