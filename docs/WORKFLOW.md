# Workflow

How data flows through CppCrystal's pipelines. Each diagram is a Mermaid
flowchart; GitHub renders them inline. Boxes name the module or function that
does the work.

## Top-level pipelines

```mermaid
flowchart TD
    Cell["Cell<br/>(lattice, positions, types, periodicity)"]

    Cell --> Det["analysis::SymmetryAnalyzer"]
    Cell --> Mag["analysis::MagneticSymmetryAnalyzer"]
    Det --> KP["SymmetryAnalyzer::reciprocal_mesh<br/>(Mesh, TimeReversal)"]

    Group["Group catalog<br/>SpaceGroup / PointGroup / RodGroup"]
    Group --> Gen["generate::Generator&lt;G&gt;"]

    Det --> Dataset["Dataset"]
    Mag --> MDataset["MagneticDataset"]
    KP --> Mesh["ReciprocalMesh<br/>(+ BrillouinZone)"]
    Gen --> Xtal["Generated<br/>(Cell + assignment)"]
```

## Space-group determination — `SymmetryAnalyzer`

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
key; the standardized cell is the idealized conventional setting, and the
recovered operations are expressed back in the original input basis.

A layer cell — one whose periodicity has an aperiodic axis — runs the same
pipeline, instantiated with `GroupFamily::layer` instead of `GroupFamily::space`.
The family is a compile-time template parameter, so every stage's family
difference is an `if constexpr`, not a runtime test on a Hall number's sign;
`dispatch_family()` at the top is the single runtime branch.

### The analyzer facade

`SymmetryAnalyzer` wraps this pipeline and memoizes it. Each public getter is a
projection of one cached dataset (or of an independently cached intermediate):

```mermaid
flowchart LR
    SA["Analyzer&lt;Derived, Traits&gt;<br/>(owns Cell + tolerances)"]
    SA -.cached.-> DS["dataset()"]
    DS --> N["hall()"]
    DS --> W["sites()"]
    DS --> ST["standardized_cell&lt;Setting, Idealize&gt;()"]
    DS --> RM["reciprocal_mesh(Mesh, TimeReversal)"]
    SA -.cached.-> PR["primitive()"]
    PR --> PC["primitive_cell()"]
```

`SymmetryAnalyzer` and `MagneticSymmetryAnalyzer` are the two `Derived` types:
the memo, the projections and the `const &`-qualified accessors live once in
`Analyzer`, and each `Traits` supplies the cell type, the dataset type and the
one pipeline call that fills it. The projections return references into the memo,
so the rvalue overloads are deleted — the analyzer must outlive what it hands
back.

## Structure generation — `Generator<G>`

Generation is the determination pipeline run in reverse: pick a group, seat a
composition on its Wyckoff sites, and materialize coordinates.

```mermaid
flowchart TD
    A["Generator&lt;G&gt;<br/>(group + GenerateOptions)"] --> B["enumerate_assignments<br/>(valid seatings of the atoms<br/>on Wyckoff positions)"]
    B -->|none fit| Err["error: composition cannot be placed"]
    B --> R{"Placement::general_only?"}
    R --> C["shuffle, then for each assignment"]
    C --> D["GroupTraits&lt;G&gt;::lattice<br/>(the family's random metric)"]
    D --> E["seed from GroupTraits&lt;G&gt;::seed_box,<br/>project onto the Wyckoff locus,<br/>expand the orbit"]
    E --> F{"distances_valid<br/>(minimum-distance criterion)"}
    F -->|fail| D
    F -->|pass| G["Generated<br/>(Cell + assignment)"]
```

One generator covers every family. `GroupTraits<G>` answers the only four
questions that differ — what to call the family in an error, how its cell is
periodic, where in that cell a seed coordinate belongs, and the random metric one
attempt is drawn in — and a layer group is not a separate specialisation but a
`SpaceGroup` whose Hall key names the layer family. Because the cell carries its
own periodicity, the orbit expansion and the distance check need no family branch
at all: a cluster is a cell with no periodic axis, and the same
`distances_valid(Cell)` is its all-pairs Euclidean check.

## Magnetic analysis — `MagneticSymmetryAnalyzer`

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

## Reciprocal-space reduction — `kpoint::ReciprocalMesh`

```mermaid
flowchart TD
    A["Mesh::of(divisions, shift)"] -->|non-positive| Err["nullopt — an invalid<br/>mesh is unrepresentable"]
    A --> B["SymmetryAnalyzer::reciprocal_mesh<br/>(the cell's rotations)"]
    B --> C["ReciprocalMesh::from_rotations<br/>(transpose ops; add −op for time reversal)"]
    C --> D["Mesh::addresses<br/>(enumerate the mesh)"]
    D --> E["map each grid point to its<br/>orbit-minimum representative"]
    E --> F["ReciprocalMesh: num_irreducible + mapping"]
    F --> G["stabilized(qpoints)<br/>brillouin_zone(lattice)"]
```

The irreducible representative is taken as the orbit minimum directly: the
reciprocal group is closed, so one pass over the rotations reaches every point in
an orbit, giving a canonical mapping independent of iteration order. `Mesh` is
constexpr throughout — its index/address round trip is `static_assert`ed in the
header rather than tested at runtime — and it carries its own shift, so the
`(mesh, is_shift)` pair no longer travels separately.

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
and cached. `cppcrystal::warmup()` builds every Hall setting of both families up
front — about 30 ms — so the caches are race-free for concurrent reads and the
cost is off the query path.
