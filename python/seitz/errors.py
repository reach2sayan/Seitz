"""The exception hierarchy, re-exported from the extension module.

Every error the library reports arrives as one of these.  The C++ side returns
``Result<T>`` and never throws; the binding layer performs the one translation,
mapping each typed error tag to the class below and attaching its payload as an
attribute.
"""

from __future__ import annotations

from . import _core

SeitzError = _core.SeitzError
SpacegroupSearchFailedError = _core.SpacegroupSearchFailedError
CellStandardizationFailedError = _core.CellStandardizationFailedError
SymmetryOperationSearchFailedError = _core.SymmetryOperationSearchFailedError
MagneticSymmetrySearchFailedError = _core.MagneticSymmetrySearchFailedError
PointgroupNotFoundError = _core.PointgroupNotFoundError
NiggliFailedError = _core.NiggliFailedError
DelaunayFailedError = _core.DelaunayFailedError
EmptyCellError = _core.EmptyCellError

#: Carries ``.determinant``, the value that made the basis singular.
InvalidLatticeError = _core.InvalidLatticeError

#: Carries ``.distance``, how close the offending pair actually was.
AtomsTooCloseError = _core.AtomsTooCloseError

__all__ = [
    "AtomsTooCloseError",
    "CellStandardizationFailedError",
    "SeitzError",
    "DelaunayFailedError",
    "EmptyCellError",
    "InvalidLatticeError",
    "MagneticSymmetrySearchFailedError",
    "NiggliFailedError",
    "PointgroupNotFoundError",
    "SpacegroupSearchFailedError",
    "SymmetryOperationSearchFailedError",
]
