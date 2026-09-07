"""Modern C++23 crystallography and symmetry analysis.

A structure in, its symmetry out::

    import numpy as np
    import seitz as cc

    cell = cc.Cell(
        cc.Lattice(3.0 * np.eye(3)),
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        [26, 26],
    )
    analyzer = cc.analyze(cell)
    print(analyzer.spacegroup_type.international_short)  # Im-3m

The analyzer memoizes, so it is the object you keep rather than a call you
repeat, and every query on it is thread-safe.  Errors are never sentinels: a
fallible call raises a :class:`~seitz.errors.SeitzError` subclass, and
"absent" is always ``None``.  :mod:`seitz.results` is the same surface with the
failure returned as an ``Ok``/``Err`` value instead, for callers who want the
C++ layer's ``Result<T>`` discipline.

Layer groups are not a separate entry point.  A cell with one aperiodic axis --
``periodicity=cc.aperiodic_along(2)`` -- goes through the same analyzer, with
the family carried by the Hall key it resolves to.
"""

from __future__ import annotations

import os
from pathlib import Path

from . import _core, errors, records
from ._core import (
    AxisKind,
    BrillouinZone,
    Cell,
    CellSetting,
    Centering,
    CIF_SYMPREC,
    CifBlock,
    CifStructure,
    CrystalClass,
    Dataset,
    GroupFamily,
    HallNumber,
    Holohedry,
    K_DEFAULT_SYMPREC,
    K_ZERO_PREC,
    Lattice,
    LatticeSetting,
    Laue,
    MagneticCell,
    MagneticDataset,
    MagneticMatch,
    MagneticOperations,
    MagneticSpacegroupType,
    MagneticSymmetryAnalyzer,
    MagneticSymmetryOperation,
    MagneticType,
    Mesh,
    OccupancyCollapse,
    Operations,
    PointGroupMatch,
    PointGroupType,
    ReciprocalMesh,
    Setting,
    Site,
    SiteTensor,
    SpaceGroup,
    SpacegroupMatch,
    SpacegroupType,
    SubgroupEdge,
    SubgroupKind,
    SymmetryAnalyzer,
    SymmetryOperation,
    TensorKind,
    TimeReversal,
    UniNumber,
    Warm,
    Wyckoff,
    all_periodic,
    aperiodic_along,
    aperiodic_axis,
    conjugated_by,
    default_hall,
    default_halls_with_pointgroup,
    elements,
    family_of,
    halls_with_number,
    magnetic_operations_from_database,
    magnetic_spacegroup_type,
    magnetic_std_transformations,
    minimal_image,
    none_periodic,
    operations_from_database,
    parse_cif,
    periodic_along,
    pointgroup_by_number,
    same_operation,
    spacegroup_type,
    subgroups,
    to_positions,
    uni_candidates,
    version_string,
    warmup,
    wrap,
)
from .errors import SeitzError
from .options import MagneticTolerance, Tolerance
from .records import CellRecord, DatasetRecord, MagneticDatasetRecord

__version__ = _core.__version__

__all__ = [
    "AxisKind",
    "BrillouinZone",
    "Cell",
    "CellRecord",
    "CellSetting",
    "Centering",
    "CIF_SYMPREC",
    "CifBlock",
    "CifStructure",
    "CrystalClass",
    "Dataset",
    "DatasetRecord",
    "GroupFamily",
    "HallNumber",
    "Holohedry",
    "K_DEFAULT_SYMPREC",
    "K_ZERO_PREC",
    "Lattice",
    "LatticeSetting",
    "Laue",
    "MagneticCell",
    "MagneticDataset",
    "MagneticDatasetRecord",
    "MagneticMatch",
    "MagneticOperations",
    "MagneticSpacegroupType",
    "MagneticSymmetryAnalyzer",
    "MagneticSymmetryOperation",
    "MagneticTolerance",
    "MagneticType",
    "Mesh",
    "OccupancyCollapse",
    "Operations",
    "PointGroupMatch",
    "PointGroupType",
    "ReciprocalMesh",
    "SeitzError",
    "Setting",
    "Site",
    "SiteTensor",
    "SpaceGroup",
    "SpacegroupMatch",
    "SpacegroupType",
    "SubgroupEdge",
    "SubgroupKind",
    "SymmetryAnalyzer",
    "SymmetryOperation",
    "TensorKind",
    "TimeReversal",
    "Tolerance",
    "UniNumber",
    "Warm",
    "Wyckoff",
    "all_periodic",
    "analyze",
    "analyze_magnetic",
    "aperiodic_along",
    "aperiodic_axis",
    "conjugated_by",
    "default_hall",
    "default_halls_with_pointgroup",
    "elements",
    "errors",
    "family_of",
    "halls_with_number",
    "magnetic_operations_from_database",
    "magnetic_spacegroup_type",
    "magnetic_std_transformations",
    "minimal_image",
    "none_periodic",
    "operations_from_database",
    "parse_cif",
    "periodic_along",
    "pointgroup_by_number",
    "read_cif",
    "records",
    "results",
    "same_operation",
    "spacegroup_type",
    "subgroups",
    "to_positions",
    "uni_candidates",
    "version_string",
    "warmup",
    "wrap",
    "write_cif",
]


def analyze(cell: Cell, tolerance: Tolerance | dict[str, float | None] | None = None, *, setting: HallNumber | None = None,) \
        -> SymmetryAnalyzer:
    """Determine the symmetry of ``cell``.

    Takes a :class:`~seitz.options.Tolerance`, a dict, or nothing.  An unset
    ``setting`` searches every Hall setting of the cell's family; a set one
    fixes it.

    Nothing is computed here: the analyzer computes on first query and caches.
    """
    return SymmetryAnalyzer.from_cell(cell, _tolerance(tolerance).to_core(), setting)


def analyze_magnetic(cell: MagneticCell,
                     tolerance: MagneticTolerance | dict[str, float | None] | None = None) \
        -> MagneticSymmetryAnalyzer:
    """Determine the magnetic symmetry of ``cell``.

    The magnetic counterpart of :func:`analyze`, and the same bargain: nothing
    is computed until a query, and every query is memoized and thread-safe.  A
    :class:`MagneticCell` carries the moments, so the rank -- collinear or
    non-collinear -- is the shape of the array handed to it, never a flag.
    """
    return MagneticSymmetryAnalyzer.from_cell(cell, _magnetic_tolerance(tolerance).to_core())


def _magnetic_tolerance(tolerance: MagneticTolerance | dict[str, float | None] | None) -> MagneticTolerance:
    if tolerance is None: return MagneticTolerance()
    if isinstance(tolerance, MagneticTolerance): return tolerance
    return MagneticTolerance.model_validate(tolerance)


def _tolerance(tolerance: Tolerance | dict[str, float | None] | None) -> Tolerance:
    """A tolerance argument as a validated model; shared by analyze/read_cif."""
    if tolerance is None: return Tolerance()
    if isinstance(tolerance, Tolerance): return tolerance
    return Tolerance.model_validate(tolerance)


def read_cif(source: str | os.PathLike[str],
             tolerance: Tolerance | dict[str, float | None] | None = None) -> list[CifStructure]:
    """Read every structure in a CIF document.

    A :class:`os.PathLike` ``source`` is read from disk, a ``str`` is the text
    itself -- by type, so a one-line CIF and a filename cannot be confused.

    Defaults to :data:`CIF_SYMPREC` (1e-3 A), not the search default:
    five-decimal coordinates put two images of a site ~1e-4 apart.
    """
    text = Path(source).read_text() if isinstance(source, os.PathLike) else source
    validated = _tolerance(tolerance) if tolerance is not None else Tolerance(symprec=CIF_SYMPREC)
    return _core.read_cif(text, validated.to_core())


def write_cif(obj: Cell | SymmetryAnalyzer, *, name: str = "seitz",
              path: str | os.PathLike[str] | None = None) -> str:
    """Render ``obj`` as a CIF document, and optionally write it to ``path``.

    A :class:`Cell` is written in P1.  A :class:`SymmetryAnalyzer` is written
    symmetrized: standardized cell, its setting's database operations, and one
    atom per orbit with Wyckoff letter and multiplicity.

    The text is returned either way, so ``print(seitz.write_cif(cell))`` is the
    printable form.
    """
    if isinstance(obj, Cell): text = _core.write_cif_cell(obj, name)
    elif isinstance(obj, SymmetryAnalyzer): text = _core.write_cif_analyzer(obj, name)
    else: raise TypeError(f"write_cif takes a Cell or a SymmetryAnalyzer, not {type(obj).__name__}")
    if path is not None: Path(path).write_text(text)
    return text


# Imported last, not with the others at the top: seitz.results wraps read_cif
# and write_cif, which are defined above, so a top-of-file import would reach
# for names this module has not bound yet.
from . import results  # noqa: E402
