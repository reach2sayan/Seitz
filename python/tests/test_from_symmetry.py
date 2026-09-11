"""Determination from a set of operations, with no atomic positions."""

from __future__ import annotations

import numpy as np
import pytest

import pyseitz as sz


def test_the_operations_of_a_cell_find_its_own_space_group(bcc_fe: sz.Cell) -> None:
    """The round trip: determine a cell, then hand the operations back alone."""
    analyzer = sz.analyze(bcc_fe)
    match = analyzer.operations.spacegroup(bcc_fe.lattice)

    assert match.type.number == 229
    assert match.hall.family == sz.GroupFamily.space


def test_a_database_setting_identifies_itself() -> None:
    """Operations straight out of the database, matched in their own basis."""
    hall = sz.default_hall(sz.GroupFamily.space, 225)
    assert hall is not None
    operations = sz.operations_from_database(hall)

    match = operations.spacegroup(sz.Lattice(4.0 * np.eye(3)))
    assert match.type.number == 225


def test_the_point_group_of_a_rotation_set(bcc_fe: sz.Cell) -> None:
    match = sz.analyze(bcc_fe).operations.point_group()

    assert match.type.symbol == "m-3m"
    assert match.type.number == 32
    assert match.transformation.shape == (3, 3)
    assert round(np.linalg.det(match.transformation)) != 0


def test_a_point_group_search_that_finds_nothing_raises() -> None:
    """A rotation set that is not a crystallographic group is an error, not a 0."""
    operations = sz.Operations([sz.SymmetryOperation(2 * np.eye(3, dtype=np.int32))])
    with pytest.raises(sz.errors.PointgroupNotFoundError):
        operations.point_group()


def test_to_primitive_recovers_the_centering(bcc_fe: sz.Cell) -> None:
    """Im-3m is body-centred, so the transformation has |det| = 2."""
    primitive = sz.analyze(bcc_fe).cell_operations.to_primitive()
    assert primitive is not None
    operations, transformation = primitive

    assert len(operations) == 48
    assert round(abs(np.linalg.det(transformation))) == 2
