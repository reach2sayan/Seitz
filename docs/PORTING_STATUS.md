# CppCrystal — porting status & resume roadmap

A modern-C++ port of spglib's C core (Eigen + Boost + STL; `boost::leaf::result<T>`
errors; no Python/Ruby/Fortran). Reference being ported/validated against:
**spglib v2.7.0**, cloned at `~/.cache/spglib-ref`.

## Build & test

```bash
cmake --preset default && ctest --test-dir cmake-build-debug          # unit tests
cmake --preset oracle  && ctest --test-dir cmake-build-debug          # + diff vs reference spglib (v2.7.0)
```
`cmake --preset oracle` FetchContent-builds the reference spglib (from the local
clone) as `Spglib::symspg` and links it only into `cppcrystal_oracle_tests`.
Always check **both** Debug+oracle and a Release build (`-O3` once exposed an
Eigen link bug Debug masked).

## Status (as of this checkpoint — 91 Debug+oracle / 53 Release, all green)

| Phase | State | Modules | Oracle check |
|---|---|---|---|
| 0 Foundations | ✅ done | `core/{types,error,tolerance,version}`, `math/{integer_matrix,fractional}` | unit |
| 1 Data model | ✅ done | `core/{cell,symmetry_operation,point_group,overlap}` | `spg_get_symmetry` ops overlap |
| 2 Reductions | ✅ done | `reduce/{delaunay,niggli}` (incl. 2D `delaunay_reduce_2d`) | match reference to 1e-8 |
| 3 Symmetry | ✅ done | `symmetry/{find_symmetry,pointgroup,primitive}` (incl. `reduce_symmetry`) | `get_symmetry`/`get_pointgroup`/`find_primitive` |
| 4 Database | ✅ done | `data/spg_database` (+ generated `data/spacegroup_{metadata,operation}_tables.hpp`) | all 530 Hall #s = `spg_get_symmetry_from_database` + `spg_get_spacegroup_type` |
| 5 Determination | ✅ **done** | `spacegroup/{hall_symbol,spacegroup}`, `dataset` (`get_dataset`) | `get_dataset` SG#/Hall#/intl symbol/hall symbol/choice + op count = `spg_get_dataset` over a 17-structure corpus (all systems + centerings) |
| 6 Standardize + Wyckoff | ✅ **done** | `data/sitesym_database`, `refine/{refinement,site_symmetry,standardize,operations}`, `dataset` | `get_dataset` std_lattice / std_rotation_matrix / transformation_matrix / std cell (set) / wyckoffs / equivalent_atoms / site_symmetry_symbols / **operations** = `spg_get_dataset` over the 17-structure corpus |
| 7 Magnetic | ✅ **done** | `data/msg_database`, `spin`, symmetry→SG pipeline, `magnetic` (MSG id + transform), `magnetic_dataset` (`get_magnetic_dataset`) | end-to-end `get_magnetic_dataset` (UNI#/type/Hall#, transformation matrix / origin shift / std rotation, std lattice + positions + tensors, operations, equiv atoms, primitive lattice) = `spg_get_magnetic_dataset` over FM / AFM / non-collinear / type-II grey |
| 8 K-points / BZ | ✅ **done** | `kpoint/{grid, reciprocal_mesh, brillouin_zone}` | grid / IR mesh / stabilized / grid-pts-by-rot / BZ relocation / BZ-grid-pts-by-rot all vs `spg_get_dense_*` (cubic/bcc/fcc/hexagonal, ±shift/±time-reversal) |

Phase 6 is complete, including the exact database-derived refined operations
(`refine/operations.cpp` = `get_refined_symmetry_operations` +
`recover_symmetry_in_original_cell`): `Dataset::operations` are the exact
operations recovered in the input cell, matching `spg_get_dataset`.

Code conventions enforced throughout (see also the memory file
`cppcrystal-port-plan`): no magic-number sentinels — `std::optional` for
"absent" (never `-1`), enums/`variant` for modes, lookup tables over
`switch`-on-convention. Big generated data tables are transcribed by
`tools/transcribe_*.py`, never hand-typed, and split into a public metadata
header + a private operation header (the encoded-op packing never reaches a
public header). Each DB module decodes its metadata once at compile time into a
`constexpr` catalog (direct index by the contiguous key — Hall # 1..530, UNI #
1..1651 — plus any sorted secondary index), and decodes the Eigen-valued
operations on demand (lazily-built per-key cache for the 1D non-magnetic DB;
by-value per call for the sparse 2D magnetic (UNI, Hall) DB). This replaced the
earlier `boost::multi_index` + `flat_map` design. The Hall-symbol generator /
VSpU tables likewise live in generated `data/hall_generators.hpp` (Eigen views
in `hall_generators_view.hpp`), and the crystal-system / rhombohedral
classification of a Hall # is a compile-time table in `data/hall_classification.hpp`.

## Phase 5 — DONE (the determination → `get_dataset`)

`spglib::get_dataset(cell, symprec)` is implemented and oracle-verified. The
ported pieces (all 3D path, layer/aperiodic branches dropped):
- `reduce/delaunay.cpp` → `delaunay_reduce_2d` (rank-2, unique axis b).
- `symmetry/find_symmetry.cpp` → `reduce_symmetry` (= `sym_reduce_operation`).
- `symmetry/primitive.hpp` → `Primitive` grew `orig_lattice` + `tolerance` +
  `angle_tolerance` (set in `find_primitive`).
- `spacegroup/spacegroup.{hpp,cpp}` → the whole `spacegroup.c` orchestration:
  `get_conventional_symmetry`, `get_centering`/`get_base_center`/centering
  shifts, `is_equivalent_lattice`, `change_basis_tricli` (Niggli) /
  `change_basis_monocli` (2D Delaunay), `match_hall_symbol_db` + per-system
  matchers (`monocli`/`ortho`/`cubic`/`rhombo`/`others`/`change_of_basis_loop`),
  `search_hall_number` → `iterative_search` → `search_spacegroup` →
  `Spacegroup{type, bravais_lattice, origin_shift}`. All the small
  change-of-basis / centering tables are hand-written `static` Eigen arrays in
  the .cpp (they are algorithm constants, not the big generated DB tables).
- `dataset.{hpp,cpp}` → `Dataset` value type + `get_dataset` (the
  `determination.c` outer loop + `spg_get_dataset` assembly, minus refinement).
- Oracle: `tests/test_oracle_spacegroup.cpp` diffs SG#/Hall#/intl/hall-symbol/
  choice + op count vs `spg_get_dataset` over 17 structures (`reference_dataset`
  added to `tests/oracle.hpp`).

## Phase 6 — DONE (standardization + Wyckoff)

`get_dataset` now returns the full standardized dataset. Ported (3D path):
- `tools/transcribe_sitesym_database.py` → `data/sitesym_tables.hpp`
  (`position_wyckoff` by Hall #, `coordinates_first`, `multiplicities`,
  `site_symmetry_symbols`); decoder `data/sitesym_database.{hpp,cpp}`
  (`wyckoff_coordinate`, `wyckoff_indices`, `site_symmetry_symbol`).
- `refine/refinement.{hpp,cpp}`: `conventional_lattice`
  (`ref_get_conventional_lattice` + the 7 `set_*` metric→lattice setters),
  `find_similar_bravais_lattice`, `measure_rigid_rotation`/orthonormal basis.
- `refine/site_symmetry.{hpp,cpp}`: `exact_positions` (= `ssm_get_exact_positions`):
  exact special positions, equivalent atoms, Wyckoff letter + site-symmetry
  symbol lookup.
- `refine/standardize.{hpp,cpp}`: `get_wyckoff_positions` (conventional
  primitive → expand across pure translations → bravais cell + per-atom Wyckoff
  data + crystallographic orbits + supercell broken-symmetry path).
- `refine/operations.{hpp,cpp}`: `refined_operations` (= `get_refined_symmetry_
  operations`): conventional DB ops → origin shift → primitive setting →
  recovered + replicated across pure translations in the input cell.
- `dataset.{hpp,cpp}`: `Dataset` grew the std cell, transformation /
  std_rotation matrices, wyckoffs, site_symmetry_symbols, equivalent_atoms,
  crystallographic_orbits, std/primitive mappings; `get_dataset` wires it all.

## Phase 7 — IN PROGRESS (magnetic)

**Magnetic database — DONE** (the Phase-4 analogue, the foundation for the rest
of Phase 7). Ported (3D path):
- `core/magnetic_symmetry_operation.hpp`: `MagneticSymmetryOperation`
  {`Matrix3i rotation`, `Vector3d translation`, `bool time_reversal`} (anti =
  true) + `MagneticSymmetryOperations`. Replaces spglib's parallel
  `MagneticSymmetry.rot/.trans/.timerev` arrays.
- `tools/transcribe_msg_database.py` → generated `data/magnetic_spacegroup_
  metadata_tables.hpp` (`kMagneticSpacegroupTypes` 1652, `kMagneticHallMapping`
  531, `kMagneticUniMapping` 1652) + `data/magnetic_spacegroup_operation_
  tables.hpp` (`kMagneticOperations` 76683 encoded, `kMagneticOperationIndex`
  [1652][18][2], `kAlternativeTransformations` [1652][18][7]). The nested
  `[18][2]`/`[18][7]` C tables are parsed by a recursive brace parser and
  zero-padded to rectangular.
- `data/msg_database.{hpp,cpp}`: `constexpr MagneticSpacegroupCatalog` (direct
  index by UNI 1..1651) + `magnetic_spacegroup_type` / `uni_candidates`
  (constexpr), and runtime decoders `magnetic_operations_from_database(uni,
  hall=0)` / `magnetic_std_transformations(uni, hall=0)` (port of
  `msgdb_get_*`). Op encoding: `enc = timerev*34012224 + spgdb_enc`
  (34012224 = 3⁹·12³); the spgdb base-3/base-12 decode is replicated in the
  .cpp (spg_database's copy is in an anon namespace). Hall→offset via
  `get_hall_number_offset` returning `std::optional`. The (UNI, Hall) key space
  is 2D and sparse, so ops are decoded by value per call (no per-key cache,
  unlike the 1D non-magnetic DB).
- Oracle: `tests/test_oracle_msg_database.cpp` — type metadata + decoded ops
  (default setting AND every Hall setting via `uni_candidates`) match
  `spg_get_magnetic_spacegroup_type` / `spg_get_magnetic_symmetry_from_database`
  for all 1651 UNI numbers. **67 Debug+oracle / 50 Release, all green.**

**Magnetic symmetry search (`spin.c`) — DONE.** Ported (3D path):
- `core/magnetic_cell.hpp`: `MagneticCell` = `Cell` + `SiteTensors`
  (`std::variant<vector<double>, vector<Vector3d>>` — the variant alternative
  *is* the tensor rank, so no separate `tensor_rank` field). `rank()` /
  `scalar(i)` / `vector(i)` accessors.
- `spin/spin.{hpp,cpp}` (namespace `spglib::spin`):
  `operations_with_site_tensors(sym_nonspin, mcell, with_time_reversal,
  is_axial, symprec, angle_tol, mag_symprec=nullopt)` → `MagneticSymmetrySearch`
  {operations (input-cell basis), equivalent_atoms, permutations,
  primitive_lattice}. Ports get_operations (filter spatial ops by moment
  consistency → spin-flip sign; undetermined-with-time-reversal ops split into
  ordinary+anti), get_symmetry_permutations, get_orbits,
  `collect_pure_translations`. `mag_symprec` is `std::optional<double>` (nullopt
  → symprec; replaces spglib's `< 0` sentinel). time_reversal is a `bool` (anti
  = true). `spn_get_idealized_cell` is DEFERRED to the magnetic_spacegroup
  milestone (where it is actually used + oracle-testable).
- Exposed `symmetry::primitive_lattice_vectors` (was the anon-namespace
  `primitive_lattice` helper in primitive.cpp) = `prm_get_primitive_lattice_
  vectors`; the magnetic search rebuilds the primitive cell from the magnetic
  pure translations through it.
- New error tag `e_magnetic_symmetry_search_failed`.
- Oracle: `tests/test_oracle_spin.cpp` — ops compared as a SET (the reference's
  internal spatial-symmetry order differs from `find_symmetry`'s), equiv-atoms
  element-wise, primitive lattice via `same_lattice` (L_a⁻¹·L_b integer
  unimodular, since the basis choice can differ). Reference time_reversal from
  `spin_flips == -1` (spin_flips = 1 − 2·timerev).

**Symmetry → space-group pipeline (foundation for the MSG determination) —
DONE.** The magnetic determination identifies its family/maximal space groups by
running a space-group search on a *set of operations* (no positions); ported
that path first, independently oracle-verified:
- `symmetry::primitive_symmetry(operations, symprec)` → `{prim_sym, t_mat}`
  (port of primitive.c `prm_get_primitive_symmetry`): pure translations of the
  op set → primitive lattice in translation space (reuses `find_primitive` on a
  unit cell of those translations) → distinct rotations transformed to the
  primitive setting. Plus the earlier `primitive_lattice_vectors`.
- `spacegroup::search_spacegroup_with_symmetry(operations, prim_lattice,
  symprec)` (port of `spa_search_spacegroup_with_symmetry`): wraps the existing
  internal matcher with a one-atom primitive at the given lattice.
- `spacegroup::spacegroup_type_from_symmetry<LatticeSetting>(operations,
  lattice, symprec)` (port of spglib.c `get_hall_number_from_symmetry`):
  primitive symmetry → Niggli-reduce → Hall match. **The C `transform_lattice_
  by_tmat` int flag became a compile-time `enum class LatticeSetting
  {conventional, primitive}` template parameter (`if constexpr`, explicit
  instantiation), per the no-magic-flag rule.** `conventional` mirrors
  `spg_get_spacegroup_type_from_symmetry`; `primitive` mirrors
  `spg_get_hall_number_from_symmetry`.
- Oracle `test_oracle_spacegroup_from_symmetry.cpp`: same op set fed to both
  sides; SG#/Hall#/intl-short match `spg_get_spacegroup_type_from_symmetry`
  across P/I/F-cubic + hexagonal + tetragonal, and the primitive branch matches
  `spg_get_hall_number_from_symmetry`.

**Magnetic space-group identification (`msg_identify_magnetic_space_group_type`)
— DONE.** `magnetic/magnetic_spacegroup.{hpp,cpp}` (namespace `spglib::magnetic`):
`identify_magnetic_spacegroup_type(lattice, magnetic_symmetry, symprec)` →
`MagneticDataset` {uni_number, `MagneticType msg_type`, hall_number,
transformation_matrix, origin_shift, std_rotation_matrix}. Ports
get_reference_space_group (FSG via family / XSG via maximal — the C
`ignore_time_reversal` int became `enum class SpaceGroupKind {family, maximal}`;
both use the new `primitive_symmetry` + `search_spacegroup_with_symmetry` +
`find_similar_bravais_lattice`), get_magnetic_space_group_type (+ MagneticType
enum I–IV replacing the bare int), get_representative, get_changed_magnetic_
symmetry, get_distinct_changed_magnetic_symmetry, get_changed_pure_translations
(denominator search for non-unimodular tmat), is_equal, get_rigid_rotation, and
the DB-matching loop (`uni_candidates` → `magnetic_operations_from_database` +
`magnetic_std_transformations` → `is_equal`). Oracle `test_oracle_magnetic.cpp`:
UNI number + MSG type match `spg_get_magnetic_spacegroup_type_from_symmetry`
(which feeds the same magnetic symmetry straight into msg_identify) for
collinear FM / AFM (±time-reversal) / non-collinear / type-II grey group. NOTE:
`std_rotation_matrix` is computed but NOT yet oracle-checked (the v2.7.0
get_rigid_rotation collapses to ≈identity — see the verbatim port); it is
verified end-to-end in the std-cell milestone.

**Magnetic standardized cell + public dataset — DONE.** Phase 7 is complete.
- `spin::idealized_cell` (= `spn_get_idealized_cell`): positions + site tensors
  averaged over the magnetic ops (each atom mapped back by the op's permutation).
- `magnetic::transform_cell` (= `msg_get_transformed_cell`): primitive via
  `find_primitive_with_pure_translations` (new, = `prm_get_primitive_with_pure_
  trans`) → `tmat_prm` → replicate across the transformed pure translations →
  rotate rank-1 tensors by the rigid rotation. The identification struct was
  renamed `magnetic::MagneticTypeIdentification` to free `MagneticDataset` for
  the public type.
- Public `spglib::MagneticDataset` + `get_magnetic_dataset(MagneticCell,
  is_axial, symprec, …)` (`magnetic_dataset.{hpp,cpp}`): wires find_symmetry →
  `spin::operations_with_site_tensors` (time reversal on) →
  `identify_magnetic_spacegroup_type` → `idealized_cell` → `transform_cell`.
- Oracle `test_oracle_magnetic_dataset.cpp` vs `spg_get_magnetic_dataset`:
  UNI#/type/Hall#, **exact** transformation_matrix / origin_shift /
  std_rotation_matrix / std_lattice, std cell atoms (type+pos+tensor, set),
  operations (set), equivalent_atoms, primitive_lattice — over collinear FM/AFM,
  non-collinear AFM, and the type-II grey group. (The exact std-cell match
  retroactively confirms the verbatim `get_rigid_rotation` port.)

## Phase 8 — IN PROGRESS (k-points / BZ)

Plan: `~/.claude/plans/get-the-plan-for-lovely-gray.md`. Unified API (one set,
collapsing spglib's legacy `int` + dense `size_t` variants): `std::size_t`
indices, `std::vector` tables, value-type result structs. `mesh`/`is_shift`/grid
addresses are `Vector3i` (is_shift is a 0/1-per-axis offset used arithmetically —
user-confirmed over `std::array<bool,3>`). Default spglib macros assumed (no
`GRID_ORDER_XYZ` / `GRID_BOUNDARY_AS_NEGATIVE`). Three milestones:

- **8a — grid + reciprocal point group — DONE.** `kpoint/grid.{hpp,cpp}` (ns
  `spglib::kpoint`): `grid_point_single_mesh` (index = a2·m0·m1 + a1·m0 + a0),
  `grid_point_double_mesh`, `address_double_mesh`, `all_grid_addresses`,
  `grid_point_from_address` + file-local `modulo_i3`/`reduce_grid_address(_double)`
  folds. `kpoint/reciprocal_mesh.{hpp,cpp}`: `point_group_reciprocal`
  (transpose + optional −transpose for time reversal + dedup) and
  `point_group_reciprocal_with_q` (q-stabilizer). Oracle
  `test_oracle_kpoint_grid.cpp` vs `spg_get_dense_grid_point_from_address`
  (all addresses incl. boundary fold, meshes 4×4×4 / 3×3×3 / 2×3×4) + a unit
  `test_kpoint.cpp`.
- **8b — IR mesh + stabilized — DONE.** `ir_reciprocal_mesh` (rotation-set form
  + structure form via `get_dataset` → `Result`), `mesh_has_conventional_
  symmetry` (= check_mesh_symmetry, the duplicate eq[0]/eq[1] quirk ported
  verbatim) → normal/distortion split, `num_ir`, `stabilized_reciprocal_mesh`,
  `grid_points_by_rotations`. **The IR mapping is computed as the orbit-minimum
  representative directly** (`mapping[i] = min over rotations of the rotated grid
  point`): the reciprocal group is closed so one rotation reaches the whole
  orbit, making this branch-independent of spglib's serial-vs-OpenMP code (both
  yield the same canonical table) — so the flagged OpenMP question is moot. The
  distortion path keeps `std::array<std::int64_t,3>` divisor/long-address (the
  ONE non-`Vector3i` spot — int32 would overflow on dense meshes, which is the
  whole reason the path exists). Oracle `test_oracle_kpoint_mesh.cpp` vs
  `spg_get_dense_ir_reciprocal_mesh` / `_stabilized_` /
  `_grid_points_by_rotations` (cubic+bcc normal, hexagonal distortion,
  ±shift/±time-reversal). NB: a linter had inverted the `views::filter`
  predicate in `magnetic_spacegroup.cpp::space_group_of_magnetic_symmetry`
  (keeping the dropped ops) — fixed (predicate negated to a "keep" form).
- **8c — BZ relocation + BZ grid points by rotations — DONE.**
  `kpoint/brillouin_zone.{hpp,cpp}`: `kBzSearchSpace` (constexpr 125 `std::array`,
  Eigen-free), `relocate_BZ_grid_address` → `BzGrid{bz_grid_address,
  bz_map=vector<optional<size_t>>}` (`std::optional` replaces spglib's
  `num_bzmesh` sentinel; min image at slot i, boundary copies appended; serial,
  matching the reference's boundary-counter dependency), `bz_tolerance`
  (= get_tolerance_for_BZ_reduction), `BZ_grid_points_by_rotations`. Oracle
  `test_oracle_kpoint_bz.cpp` vs `spg_relocate_dense_BZ_grid_address` /
  `_BZ_grid_points_by_rotations` (cubic/fcc/bcc/hexagonal, ±shift). NB:
  `grid_address[i] + offset.cwiseProduct(mesh)` — NOT `offset * mesh` (Eigen
  reads `Vector3i * Vector3i` as an invalid matrix product); ported as a
  component-wise loop.

**Phase 8 complete.** Only Phase 9 (hardening + docs) remains.

## Phase 9 — Hardening (IN PROGRESS; docs deferred)

Plan: `~/.claude/plans/get-the-plan-for-lovely-gray.md`. 9a sanitizer+warnings;
9b input-validation + error-paths; 9c broad corpus (230 ref cells + jitter).

- **9a — sanitizer + warning baseline — DONE.** Added the `oracle-asan` preset
  (oracle + `SPGLIB_ENABLE_SANITIZERS` = `-fsanitize=address,undefined`). Added
  `-Wshadow -Wcast-align -Wdouble-promotion -Wnull-dereference
  -Wimplicit-fallthrough -Wconversion -Wsign-conversion` to the **library**
  (PRIVATE — the oracle tests bridge the reference C API, so the conversion
  warnings stay off them). Lib is **zero-warning** after 5 trivial
  signed-index→`size_t` casts (from the `std::views::enumerate` refactors in
  overlap / delaunay / pointgroup / hall_symbol). Full suite **91/91 green under
  ASan/UBSan**, zero warnings; Release still green.
  - **Found & fixed a real determination bug** the rebuild surfaced:
    `spacegroup.cpp::change_of_centering_monocli()` had been turned `consteval`
    (immediate function) while still returning a reference to a
    `static constexpr` table indexed by a **runtime** index — a consteval
    function's static local is never emitted to runtime storage, so the read
    returned `0` (`Centering::error`) and **C-centered monoclinic determination
    failed for every Hall number**. Fix: drop `consteval`. (valgrind/ASan don't
    flag this — compile-time data read at runtime, not heap/stack UB.)

- **9b — input validation + UB elimination + error-path tests — DONE.**
  `core/error.hpp`: added `e_invalid_lattice {double determinant;}` and
  `e_invalid_mesh {}`; deleted the orphans `e_reciprocal_mesh_failed` /
  `e_array_size_shortage`. New `core/validation.hpp::validate_cell(Cell) ->
  Result<void>` (empty → `e_empty_cell`; `|det(lattice)| < kZeroPrec` →
  `e_invalid_lattice`), called first in `get_dataset`, `find_symmetry`,
  `get_magnetic_dataset` — one up-front guard makes the downstream internal
  `Matrix3d::inverse()` calls safe. (Used an explicit `if (auto v = …; !v)
  return v.error();` rather than `BOOST_LEAF_CHECK`, which trips
  `-Wgnu-statement-expression`.) K-points: non-positive mesh → `e_invalid_mesh`
  in the Result-returning `ir_reciprocal_mesh(cell, …)`; the value-returning
  grid functions degrade to empty/0 via a `valid_mesh` guard at the two roots
  (`all_grid_addresses`, `grid_point_double_mesh`) — no more mod-by-zero /
  over-allocation. The magnetic factor-group `static_vector<…,48>` became a
  `small_vector<…,48>` (spills to heap instead of overflowing). New unit test
  `tests/test_error_paths.cpp` (no reference needed): empty / singular-lattice /
  non-positive-mesh → the exact error tag via `leaf::try_handle_all`, and the
  value functions degrade without UB. 95 Debug+oracle / 53 Release green, zero
  lib warnings.

- **9c — broad oracle corpus — DONE.** `tests/corpus.hpp` (a small hand-rolled
  parser; no YAML lib) loads spglib's reference corpus —
  `${spglib_reference_SOURCE_DIR}/test/functional/python/data/<system>/
  unitcell_<N>.yaml`, one input cell per space group (~230), via the
  `SPGLIB_REF_DATA_DIR` compile-def set in `tests/CMakeLists.txt`. Lattice rows →
  Matrix3d columns; the embedded `space_group.number` is a free orientation
  cross-check. `tests/test_oracle_corpus.cpp`: `get_dataset` vs `spg_get_dataset`
  over all 230 (SG#/Hall#/intl-symbol/op-count/n_std_atoms) **and** vs the YAML
  SG number — **all 230 match, no latent determination bug surfaced.**
  `tests/test_robustness.cpp`: SG number stable under sub-symprec position jitter
  and invariant for a 2×2×2 supercell, over a sample of the corpus. 98 Debug+
  oracle / 53 Release green (the corpus adds ~70 s of determination — full suite
  ~98 s).

### Phase 9 complete except docs. Remaining: README + public-API docs + the
### convention/footgun catalog (deferred Phase 9-docs).

## Public C-API parity — audited COMPLETE

`spg_standardize_cell` / `spg_refine_cell` / `spg_find_primitive` are ported as
`standardize.{hpp,cpp}`: `standardize_cell(Cell, StandardizeOptions{to_primitive,
no_idealize}, …)` (all four flag combos, derived from `get_dataset` — idealized
cases read the dataset std_* cell, no_idealize cases transform the input via the
dataset transformation matrix + `symmetry::trim_to_lattice` = `cell.c trim_cell`),
with inline `refine_cell` (default opts) and `find_primitive` (`{.to_primitive =
true}`) convenience wrappers. Oracle `test_oracle_standardize.cpp` vs
`spg_standardize_cell` over all four combos (metric tensor / composition / atom
set; F/I/C-centered + primitive + primitive→conventional expansion + strained
no_idealize).

All 59 reference public functions in `spglib.h` audited against the port: full
coverage. The `spgat_*` / `spgms_*` variants collapse into optional
`angle_tolerance` / `mag_symprec` params; the legacy `int` + dense `size_t`
k-point pairs collapse into one `size_t` set. By-design omissions (not gaps):
`spg_get_error_code/message` (→ `boost::leaf::result`), `spg_free_*dataset` (→
RAII value types), deprecated accessors (→ `Dataset` fields / `data::
spacegroup_type`), `spg_get_version_full`/`_commit` (no vendored git metadata).

## Then Phase 9b/9c
