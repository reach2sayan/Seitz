"""Serializable records.

The plain-data half of a determination, as frozen pydantic models: what a caller
wants to compare, print, store, or hand to something else as JSON.  The
heavyweight handles -- :class:`~seitz.SymmetryAnalyzer`,
:class:`~seitz.SpaceGroup`, :class:`~seitz.Cell` -- stay extension
types, because their lifetime rules live on the C++ object and a model that
merely held one would not inherit them.

Building a record copies.  :class:`Site` in particular is per atom, so a large
cell should reach for ``SymmetryAnalyzer.site_arrays()`` instead, which returns
the same answers as parallel NumPy arrays with no per-site object at all.
"""

from __future__ import annotations

from typing import Any, Self

import numpy as np
from pydantic import BaseModel, ConfigDict

from . import _core
from ._arrays import Basis, Positions, Vector3

__all__ = ["CellRecord", "DatasetRecord", "Setting", "Site", "SpacegroupType"]


class _Frozen(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")

    def __eq__(self, other: object) -> bool:
        """Field-by-field equality that copes with the array fields.

        Pydantic's own ``__eq__`` compares ``self.__dict__ == other.__dict__``,
        and a dict comparison over a NumPy array raises "truth value of an array
        is ambiguous" rather than answering.  Every record here holds at least
        one array, so comparing two of them would be an error instead of a
        result -- which shows up the moment anyone round-trips one through JSON
        and checks it survived.
        """
        if other.__class__ is not self.__class__:
            return NotImplemented
        for name in type(self).model_fields:
            mine, theirs = getattr(self, name), getattr(other, name)
            if isinstance(mine, np.ndarray) or isinstance(theirs, np.ndarray):
                if not np.array_equal(mine, theirs):
                    return False
            elif mine != theirs:
                return False
        return True

    # Defining __eq__ in a class body sets __hash__ to None; these are frozen
    # models, so restore it -- by identity, since the array fields are not
    # hashable and a value hash would have to copy them on every lookup.
    __hash__ = object.__hash__

    @classmethod
    def from_core(cls, obj: Any) -> Self:
        """Copy out of the extension type of the same shape, field by name.

        Every field below is named after the attribute it reads, so the model
        is the only place the list of them is written down: adding a field to
        one of these records is a one-line change, and a field the C++ side
        does not have fails here rather than silently arriving as None.
        Pydantic does the rest -- an IntEnum into an ``int``, an array through
        the annotations in ``_arrays``.
        """
        return cls(**{name: getattr(obj, name) for name in cls.model_fields})

class Site(_Frozen):
    """The per-atom result of a determination."""

    wyckoff: int
    site_symmetry: str
    equivalent_atom: int
    orbit: int
    primitive_atom: int

class Setting(_Frozen):
    """How the input cell maps onto the standardized setting."""

    transformation: Basis
    origin_shift: Vector3
    rigid_rotation: Basis

class SpacegroupType(_Frozen):
    """One Hall setting's metadata."""

    number: int
    schoenflies: str
    hall_symbol: str
    international: str
    international_full: str
    international_short: str
    choice: str
    centering: int
    pointgroup_number: int


class CellRecord(_Frozen):
    """A cell as data, for storage and provenance.
    The explicit adapter rather than making :class:`~seitz.Cell` itself a
    model: a Cell is on the hot path in both directions -- out of a
    determination and straight back into another one -- and a validation layer
    in the middle of that would be paid on every round trip for nothing.
    """

    basis: Basis
    positions: Positions
    types: list[int]
    periodicity: tuple[int, int, int]

    @classmethod
    def from_cell(cls, cell: _core.Cell) -> CellRecord:
        return cls(
            basis=cell.lattice.matrix,
            positions=cell.positions,
            types=[int(t) for t in cell.types],
            periodicity=tuple(int(axis) for axis in cell.periodicity),
        )

    def to_cell(self) -> _core.Cell:
        return _core.Cell(
            _core.Lattice(np.asarray(self.basis)),
            np.asarray(self.positions),
            list(self.types),
            tuple(_core.AxisKind(axis) for axis in self.periodicity),
        )


class DatasetRecord(_Frozen):
    """A whole determination as data.

    Eager, including ``sites``: this type exists to be serialized, so anything
    it left lazy would be built by the first ``model_dump()`` anyway.  The lazy
    path is the analyzer itself -- its projections are memoized, and
    ``site_arrays()`` answers the same per-atom questions as parallel NumPy
    arrays with no per-site object built at all, which is the route for a cell
    large enough that thousands of :class:`Site` models would show up in a
    profile.
    """

    hall_family: int
    hall_index: int
    number: int
    international_short: str
    bravais: Basis
    setting: Setting
    primitive: Basis
    standardized: CellRecord
    std_mapping_to_primitive: list[int]
    num_operations: int
    sites: tuple[Site, ...]

    @classmethod
    def from_analyzer(cls, analyzer: _core.SymmetryAnalyzer) -> DatasetRecord:
        dataset = analyzer.dataset
        spacegroup_type = _core.spacegroup_type(dataset.hall)
        return cls(
            hall_family=int(dataset.hall.family),
            hall_index=dataset.hall.index,
            number=spacegroup_type.number,
            international_short=spacegroup_type.international_short,
            bravais=dataset.bravais.matrix,
            setting=Setting.from_core(dataset.setting),
            primitive=dataset.primitive.matrix,
            standardized=CellRecord.from_cell(dataset.standardized),
            std_mapping_to_primitive=list(dataset.std_mapping_to_primitive),
            num_operations=len(dataset.operations),
            sites=tuple(Site.from_core(s) for s in dataset.sites),
        )
