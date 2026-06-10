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

## Status (as of this checkpoint — 64 Debug+oracle / 50 Release, all green)

| Phase | State | Modules | Oracle check |
|---|---|---|---|
| 0 Foundations | ✅ done | `core/{types,error,tolerance,version}`, `math/{integer_matrix,fractional}` | unit |
| 1 Data model | ✅ done | `core/{cell,symmetry_operation,point_group,overlap}` | `spg_get_symmetry` ops overlap |
| 2 Reductions | ✅ done | `reduce/{delaunay,niggli}` (incl. 2D `delaunay_reduce_2d`) | match reference to 1e-8 |
| 3 Symmetry | ✅ done | `symmetry/{find_symmetry,pointgroup,primitive}` (incl. `reduce_symmetry`) | `get_symmetry`/`get_pointgroup`/`find_primitive` |
| 4 Database | ✅ done | `data/spg_database` (+ generated `data/spacegroup_{metadata,operation}_tables.hpp`) | all 530 Hall #s = `spg_get_symmetry_from_database` + `spg_get_spacegroup_type` |
| 5 Determination | ✅ **done** | `spacegroup/{hall_symbol,spacegroup}`, `dataset` (`get_dataset`) | `get_dataset` SG#/Hall#/intl symbol/hall symbol/choice + op count = `spg_get_dataset` over a 17-structure corpus (all systems + centerings) |
| 6 Standardize + Wyckoff | ✅ **done** | `data/sitesym_database`, `refine/{refinement,site_symmetry,standardize,operations}`, `dataset` | `get_dataset` std_lattice / std_rotation_matrix / transformation_matrix / std cell (set) / wyckoffs / equivalent_atoms / site_symmetry_symbols / **operations** = `spg_get_dataset` over the 17-structure corpus |
| 7 Magnetic | ⬜ todo | — | — |
| 8 K-points / BZ | ⬜ todo | — | — |

Phase 6 is complete, including the exact database-derived refined operations
(`refine/operations.cpp` = `get_refined_symmetry_operations` +
`recover_symmetry_in_original_cell`): `Dataset::operations` are the exact
operations recovered in the input cell, matching `spg_get_dataset`.

Code conventions enforced throughout (see also the memory file
`cppcrystal-port-plan`): no magic-number sentinels — `std::optional` for
"absent" (never `-1`), enums/`variant` for modes, lookup tables / Boost maps
over `switch`-on-convention. The space-group DB uses a `boost::multi_index`
catalog (ByHall/ByNumber/ByPointgroup) + `flat_map<int,SymmetryOperations>`.
Big generated data tables are transcribed by `tools/transcribe_*.py`, never
hand-typed.

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

## Then Phases 7–8 (each its own milestone)

- **7 — Magnetic:** `spin.c`, `magnetic_spacegroup.c`; transcribe
  `msg_database.c` (1,651 UNI entries). Oracle pinned to v2.7.0
  (`spg_get_magnetic_dataset`, `spg_get_magnetic_symmetry_from_database`).
- **8 — K-points / BZ:** `kgrid.c`, `kpoint.c` (ir reciprocal mesh, stabilized
  mesh, grid points by rotation, BZ relocation). Collapse the legacy `int` and
  dense `size_t` variants into one overloaded API. Oracle: ir mapping tables.

These ~10k lines are best parallelized via a multi-agent workflow (the
transcriptions and the independent module ports fan out cleanly), then
integrated and oracle-tested one module at a time.
