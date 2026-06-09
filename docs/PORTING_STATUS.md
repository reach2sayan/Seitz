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

## Status (as of this checkpoint — 55 Debug+oracle / 46 Release, all green)

| Phase | State | Modules | Oracle check |
|---|---|---|---|
| 0 Foundations | ✅ done | `core/{types,error,tolerance,version}`, `math/{integer_matrix,fractional}` | unit |
| 1 Data model | ✅ done | `core/{cell,symmetry_operation,point_group,overlap}` | `spg_get_symmetry` ops overlap |
| 2 Reductions | ✅ done | `reduce/{delaunay,niggli}` | match reference to 1e-8 |
| 3 Symmetry | ✅ done | `symmetry/{find_symmetry,pointgroup,primitive}` | `get_symmetry`/`get_pointgroup`/`find_primitive` |
| 4 Database | ✅ done | `data/spg_database` (+ generated `data/spacegroup_tables.hpp`) | all 530 Hall #s = `spg_get_symmetry_from_database` + `spg_get_spacegroup_type` |
| 5 Determination | 🟡 **matcher done; orchestration TODO** | `spacegroup/hall_symbol` (+ generated `data/hall_generators.hpp`) | self-consistency: all 530 settings match own ops, zero shift |
| 6 Standardize + Wyckoff | ⬜ todo | — | — |
| 7 Magnetic | ⬜ todo | — | — |
| 8 K-points / BZ | ⬜ todo | — | — |

Code conventions enforced throughout (see also the memory file
`cppcrystal-port-plan`): no magic-number sentinels — `std::optional` for
"absent" (never `-1`), enums/`variant` for modes, lookup tables / Boost maps
over `switch`-on-convention. The space-group DB uses a `boost::multi_index`
catalog (ByHall/ByNumber/ByPointgroup) + `flat_map<int,SymmetryOperations>`.
Big generated data tables are transcribed by `tools/transcribe_*.py`, never
hand-typed.

## Resume: finish Phase 5 (the determination → `get_dataset`)

Goal: `spglib::get_dataset(cell, symprec)` producing `spacegroup_number`,
`hall_number`, international symbol, operations, transformation/origin_shift —
matching `spg_get_dataset`. Port `spacegroup.c` (3D path) wiring the existing
pieces (`find_primitive` → `find_symmetry` → `get_pointgroup` →
`match_hall_symbol`). Concretely, port from `~/.cache/spglib-ref/src/spacegroup.c`:

1. **Data tables** (transcribe, like Phase 4): `spacegroup_to_hall_number[230]`
   (first Hall # per space group — the candidate list) and the fixed centering
   matrices `A_mat/C_mat/F_mat/I_mat/R_mat` + `change_of_basis_monocli`.
2. **`del_layer_delaunay_reduce_2D`** (delaunay.c) — NOT yet ported; needed by
   `change_basis_monocli`. The 2D Delaunay path.
3. **`search_hall_number`** (spacegroup.c:789): `ptg_get_transformation_matrix`
   (have it) → for LAUE1/LAUE2M correct the basis (`change_basis_tricli` via
   `niggli_reduce`, have it; `change_basis_monocli` via 2D Delaunay) →
   `get_centering` (tmat+laue → Centering + correction_mat) → build
   `conv_lattice` + `get_initial_conventional_symmetry`
   (`get_conventional_symmetry`: transform primitive ops to the conventional
   setting) → loop candidates calling `match_hall_symbol_db`.
4. **`match_hall_symbol_db`** + per-system variants
   (`_monocli`/`_ortho`/`_cubic`/`_rhombo` and their `_in_loop`): these try the
   axis/setting permutations, each ultimately calling our existing
   `spacegroup::match_hall_symbol` (= `hal_match_hall_symbol_db`, done).
5. **`spa_search_spacegroup`** (spacegroup.c:451) + `get_spacegroup` → a
   `Spacegroup{number, hall_number, pointgroup_number, symbols, bravais_lattice,
   origin_shift}` value type.
6. **`determination.c`** (138 lines) + `spglib.c`'s `spg_get_dataset` assembly →
   build the `Dataset` (operations from the primitive symmetry expanded to the
   input cell; wyckoffs/equiv-atoms come in Phase 6).
   Oracle: compare `spacegroup_number`/`hall_number`/`international_short` to
   `spg_get_dataset` over the structure corpus (set/gauge-aware comparators in
   `tests/oracle.hpp`).

## Then Phases 6–8 (each its own milestone)

- **6 — Standardize + Wyckoff:** `refinement.c` (`standardize_cell`,
  `find_primitive`, `refine_cell`, std cell + transformation/origin) and
  `site_symmetry.c` (Wyckoff letters, site-symmetry symbols, equivalent atoms,
  crystallographic orbits). Transcribe `sitesym_database.c` (1,645 lines) →
  `data/sitesym_database.hpp`. Oracle: std cell + wyckoffs + equivalent_atoms.
- **7 — Magnetic:** `spin.c`, `magnetic_spacegroup.c`; transcribe
  `msg_database.c` (1,651 UNI entries). Oracle pinned to v2.7.0
  (`spg_get_magnetic_dataset`, `spg_get_magnetic_symmetry_from_database`).
- **8 — K-points / BZ:** `kgrid.c`, `kpoint.c` (ir reciprocal mesh, stabilized
  mesh, grid points by rotation, BZ relocation). Collapse the legacy `int` and
  dense `size_t` variants into one overloaded API. Oracle: ir mapping tables.

These ~10k lines are best parallelized via a multi-agent workflow (the
transcriptions and the independent module ports fan out cleanly), then
integrated and oracle-tested one module at a time.
