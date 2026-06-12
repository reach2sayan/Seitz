# Workflow

How data flows through CppCrystal's pipelines. Each diagram is a Mermaid
flowchart; GitHub renders them inline. Boxes name the module or function that
does the work.

## Top-level pipelines

```mermaid
flowchart TD
    Cell["Cell<br/>(lattice, positions, types, periodicity)"]

    Cell --> Det["Space-group determination<br/>get_dataset"]
    Cell --> Mag["Magnetic analysis<br/>get_magnetic_dataset"]
    Cell --> KP["Reciprocal-space reduction<br/>kpoint::ir_reciprocal_mesh"]

    Group["Group catalog<br/>SpaceGroup / PointGroup / RodGroup"]
    Group --> Gen["Structure generation<br/>random_crystal / random_cluster"]

    Det --> Dataset["Dataset"]
    Mag --> MDataset["MagneticDataset"]
    KP --> Mesh["Ir mesh + mapping"]
    Gen --> Xtal["GeneratedCrystal / GeneratedCluster"]
```

## Space-group determination — `get_dataset`

The core pipeline. A cell goes in; a fully classified, standardized dataset
comes out.

```mermaid
flowchart TD
    A["Cell"] --> V{"validate_cell"}
    V -->|empty / singular| Err1["e_empty_cell /<br/>e_invalid_lattice"]
    V -->|ok| B["find_symmetry<br/>(symmetry operations of the cell)"]
    B --> C["find_primitive<br/>(primitive lattice + ops)"]
    C --> D["reduce lattice<br/>Niggli (triclinic) / 2D Delaunay (monoclinic)"]
    D --> E["search_spacegroup<br/>match conventional ops to a Hall setting"]
    E -->|no match at any tolerance| Err2["e_spacegroup_search_failed"]
    E --> F["Spacegroup<br/>(type, bravais lattice, origin shift)"]
    F --> G["refine: conventional lattice<br/>+ exact Wyckoff positions"]
    G --> H["recover operations in the input cell"]
    F --> H
    G --> I["Dataset<br/>identity, operations, std cell,<br/>wyckoffs, equivalent_atoms"]
    H --> I
```

The Wyckoff data comes from the site-symmetry database keyed by the matched Hall
number; the standardized cell is the idealized conventional setting, and the
recovered operations are expressed back in the original input basis.

`get_layer_dataset` runs the same pipeline with one aperiodic axis set on the
cell; matching then draws from the layer-group settings (negative Hall numbers).

### The analyzer facade

`SymmetryAnalyzer` wraps this pipeline and memoizes it. Each public getter is a
projection of one cached dataset (or of an independently cached intermediate):

```mermaid
flowchart LR
    SA["SymmetryAnalyzer<br/>(owns Cell + tolerances)"]
    SA -.cached.-> DS["dataset()"]
    DS --> N["spacegroup_number()"]
    DS --> W["wyckoffs()"]
    DS --> ST["standardized_cell()"]
    SA -.cached.-> PR["primitive()"]
    PR --> PC["primitive_cell()"]
    SA -.cached.-> SG["spacegroup()"]
```

## Structure generation — `random_crystal`

Generation is the determination pipeline run in reverse: pick a group, seat a
composition on its Wyckoff sites, and materialize coordinates.

```mermaid
flowchart TD
    A["SpaceGroup + Composition"] --> B["list_wyckoff_combinations<br/>(valid seatings of the atoms<br/>on Wyckoff positions)"]
    B -->|none fit| Err["error: composition cannot be placed"]
    B --> C["for each combination"]
    C --> D["sample a random lattice<br/>(compatible metric, volume_factor)"]
    D --> E["sample free coordinates<br/>of each occupied Wyckoff site"]
    E --> F{"distance_check<br/>(minimum-distance criterion)"}
    F -->|fail| D
    F -->|pass| G["GeneratedCrystal<br/>(Cell + assignment)"]
```

`random_layer_crystal` and the rod-group overload share this loop with 2D/1D
lattices; `random_cluster` replaces the lattice with a point-group metric and
runs an all-pairs (non-periodic) distance check.

## Magnetic analysis — `get_magnetic_dataset`

```mermaid
flowchart TD
    A["MagneticCell<br/>(Cell + per-site tensors)"] --> B["find_symmetry<br/>(spatial operations)"]
    B --> C["spin::operations_with_site_tensors<br/>(keep ops consistent with the moments;<br/>assign time-reversal sign)"]
    C --> D["identify_magnetic_spacegroup_type<br/>(reference space group via family / maximal,<br/>then match the UNI database)"]
    D --> E["idealized_cell<br/>(average positions + tensors over the ops)"]
    E --> F["transform_cell<br/>(standardized magnetic cell)"]
    D --> G["MagneticDataset<br/>(UNI #, type, transform, std cell + tensors,<br/>operations, equivalent_atoms)"]
    F --> G
```

## Reciprocal-space reduction — `kpoint::ir_reciprocal_mesh`

```mermaid
flowchart TD
    A["Cell + mesh + shift + time_reversal"] --> V{"valid mesh?"}
    V -->|non-positive| Err["e_invalid_mesh"]
    V -->|ok| B["get_dataset<br/>(rotations of the cell)"]
    B --> C["point_group_reciprocal<br/>(transpose ops; add −op for time reversal)"]
    C --> D["all_grid_addresses<br/>(enumerate the mesh)"]
    D --> E["map each grid point to its<br/>orbit-minimum representative"]
    E --> F["Ir mesh: num_ir + mapping<br/>(+ stabilized / BZ variants)"]
```

The irreducible representative is taken as the orbit minimum directly: the
reciprocal group is closed, so one pass over the rotations reaches every point in
an orbit, giving a canonical mapping independent of iteration order.

## Where the symmetry data comes from

Every pipeline above reads from compiled-in catalogs rather than parsing files
at runtime.

```mermaid
flowchart LR
    Ref["spglib data tables (C source)"] --> Py["tools/transcribe_*.py"]
    Py --> Meta["generated metadata header<br/>(constexpr, public)"]
    Py --> Ops["generated operation header<br/>(encoded, private)"]
    Meta --> Cat["constexpr catalog<br/>(direct index by Hall # / UNI #)"]
    Ops --> Dec["lazy runtime decode<br/>(Eigen operation matrices, cached)"]
    Cat --> Use["determination · generation ·<br/>magnetic · k-points"]
    Dec --> Use
```

Metadata (numbers, symbols, centering, multiplicities) is decoded once at
compile time and indexed directly by the group's key. Operation matrices, which
cannot be `constexpr` because they are `Eigen` types, are decoded on first use
and cached. `spglib::warmup()` primes these caches up front so they are race-free
for concurrent reads.
