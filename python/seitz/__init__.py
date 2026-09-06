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
"absent" is always ``None``.

Layer groups are not a separate entry point.  A cell with one aperiodic axis --
``periodicity=cc.aperiodic_along(2)`` -- goes through the same analyzer, with
the family carried by the Hall key it resolves to.
"""

from __future__ import annotations

import os
from pathlib import Path

from . import _core, errors, records
from ._core import (
    CIF_SYMPREC,
    K_DEFAULT_SYMPREC,
    K_ZERO_PREC,
    AxisKind,
    Cell,
    CifBlock,
    CifStructure,
    Centering,
    CellSetting,
    CrystalClass,
    Dataset,
    GroupFamily,
    HallNumber,
    Holohedry,
    Lattice,
    LatticeSetting,
    Laue,
    OccupancyCollapse,
    Operations,
    PointGroupType,
    Setting,
    Site,
    SpaceGroup,
    SpacegroupMatch,
    SpacegroupType,
    SymmetryAnalyzer,
    SymmetryOperation,
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
    minimal_image,
    none_periodic,
    parse_cif,
    operations_from_database,
    periodic_along,
    pointgroup_by_number,
    same_operation,
    spacegroup_type,
    SubgroupEdge,
    SubgroupKind,
    subgroups,
    to_positions,
    version_string,
    warmup,
    wrap,
)
from .errors import SeitzError
from .options import MagneticTolerance, Tolerance
from .records import CellRecord, DatasetRecord

__version__ = _core.__version__

__all__ = [
    "CIF_SYMPREC",
    "K_DEFAULT_SYMPREC",
    "K_ZERO_PREC",
    "AxisKind",
    "Cell",
    "CellRecord",
    "CifBlock",
    "CifStructure",
    "CellSetting",
    "Centering",
    "SeitzError",
    "CrystalClass",
    "Dataset",
    "DatasetRecord",
    "GroupFamily",
    "HallNumber",
    "Holohedry",
    "Lattice",
    "LatticeSetting",
    "Laue",
    "MagneticTolerance",
    "OccupancyCollapse",
    "Operations",
    "PointGroupType",
    "Setting",
    "Site",
    "SpaceGroup",
    "SpacegroupMatch",
    "SpacegroupType",
    "SymmetryAnalyzer",
    "SymmetryOperation",
    "TimeReversal",
    "Tolerance",
    "UniNumber",
    "Warm",
    "Wyckoff",
    "all_periodic",
    "analyze",
    "aperiodic_along",
    "aperiodic_axis",
    "conjugated_by",
    "default_hall",
    "default_halls_with_pointgroup",
    "elements",
    "errors",
    "family_of",
    "halls_with_number",
    "minimal_image",
    "none_periodic",
    "parse_cif",
    "read_cif",
    "operations_from_database",
    "periodic_along",
    "pointgroup_by_number",
    "records",
    "same_operation",
    "spacegroup_type",
    "SubgroupEdge",
    "SubgroupKind",
    "subgroups",
    "to_positions",
    "version_string",
    "warmup",
    "wrap",
    "write_cif",
]


def analyze(cell: Cell, tolerance: Tolerance | dict[str, float | None] | None = None, *, setting: HallNumber | None = None,) \
        -> SymmetryAnalyzer:
    """Determine the symmetry of ``cell``.

    Accepts a validated :class:`~seitz.options.Tolerance`, a plain dict, or
    nothing at all; :func:`_tolerance` is where any of the three becomes the
    struct the extension takes.

    An unset ``setting`` searches every Hall setting of the cell's family; a set
    one fixes it.

    Nothing is computed here.  The returned analyzer computes on first query and
    caches, so asking it several questions costs one determination.
    """
    return SymmetryAnalyzer.from_cell(cell, _tolerance(tolerance).to_core(), setting)


def _tolerance(tolerance: Tolerance | dict[str, float | None] | None) -> Tolerance:
    """The one place a tolerance argument becomes a validated model.

    Shared by :func:`analyze` and :func:`read_cif`, which differ only in the
    default they fall back to when the argument is ``None``.
    """
    if tolerance is None: return Tolerance()
    if isinstance(tolerance, Tolerance): return tolerance
    return Tolerance.model_validate(tolerance)


def read_cif(source: str | os.PathLike[str],
             tolerance: Tolerance | dict[str, float | None] | None = None) -> list[CifStructure]:
    """Read every structure in a CIF document.

    ``source`` is either the CIF text itself or a path to a file holding it; a
    :class:`os.PathLike` is always read, a ``str`` is always the text.  The
    split is by type rather than by guessing, so a one-line CIF and a filename
    can never be confused.

    The default tolerance is :data:`CIF_SYMPREC` (1e-3 A), not the search
    default: CIF coordinates are written to five decimals, which puts two
    images of one site as much as 1e-4 apart.
    """
    text = Path(source).read_text() if isinstance(source, os.PathLike) else source
    validated = _tolerance(tolerance) if tolerance is not None else Tolerance(symprec=CIF_SYMPREC)
    return _core.read_cif(text, validated.to_core())


def write_cif(obj: Cell | SymmetryAnalyzer, *, name: str = "seitz",
              path: str | os.PathLike[str] | None = None) -> str:
    """Render ``obj`` as a CIF document, and optionally write it to ``path``.

    A :class:`Cell` is written in P1, every atom spelled out.  A
    :class:`SymmetryAnalyzer` is written symmetrized: its standardized cell, the
    database operations of the setting it determined, and one representative
    atom per orbit with its Wyckoff letter and multiplicity.

    The text is returned either way, so this is also the printable form::

        print(seitz.write_cif(cell))
    """
    if isinstance(obj, Cell): text = _core.write_cif_cell(obj, name)
    elif isinstance(obj, SymmetryAnalyzer): text = _core.write_cif_analyzer(obj, name)
    else: raise TypeError(f"write_cif takes a Cell or a SymmetryAnalyzer, not {type(obj).__name__}")
    if path is not None: Path(path).write_text(text)
    return text
