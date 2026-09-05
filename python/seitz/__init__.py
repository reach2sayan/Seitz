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

from . import _core, errors, records
from ._core import (
    K_DEFAULT_SYMPREC,
    K_ZERO_PREC,
    AxisKind,
    Cell,
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
    operations_from_database,
    periodic_along,
    pointgroup_by_number,
    same_operation,
    spacegroup_type,
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
    "K_DEFAULT_SYMPREC",
    "K_ZERO_PREC",
    "AxisKind",
    "Cell",
    "CellRecord",
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
    "operations_from_database",
    "periodic_along",
    "pointgroup_by_number",
    "records",
    "same_operation",
    "spacegroup_type",
    "subgroups",
    "to_positions",
    "version_string",
    "warmup",
    "wrap",
]


def analyze(cell: Cell, tolerance: Tolerance | dict[str, float | None] | None = None, *, setting: HallNumber | None = None,) \
        -> SymmetryAnalyzer:
    """Determine the symmetry of ``cell``.

    This is the single place a validated :class:`~seitz.options.Tolerance`
    becomes the plain struct the extension takes, which is why it accepts a
    model, a plain dict, or nothing at all.

    An unset ``setting`` searches every Hall setting of the cell's family; a set
    one fixes it.

    Nothing is computed here.  The returned analyzer computes on first query and
    caches, so asking it several questions costs one determination.
    """
    if tolerance is None: validated = Tolerance()
    elif isinstance(tolerance, Tolerance): validated = tolerance
    else: validated = Tolerance.model_validate(tolerance)
    return SymmetryAnalyzer.from_cell(cell, validated.to_core(), setting)
