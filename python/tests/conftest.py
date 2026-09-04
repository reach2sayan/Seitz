"""Fixtures mirroring tests/helpers.hpp, so the two suites test the same cells."""

from __future__ import annotations

import numpy as np
import pytest

import cppcrystal as cc


@pytest.fixture
def bcc_fe() -> cc.Cell:
    """The demo cell from main.cpp: 2 Fe in a 3 A cube -> Im-3m, No. 229."""
    return cc.Cell(
        cc.Lattice(3.0 * np.eye(3)),
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        [26, 26],
    )


@pytest.fixture
def simple_cubic() -> cc.Cell:
    """One atom in a cube -> Pm-3m, No. 221."""
    return cc.Cell(cc.Lattice(np.eye(3) * 4.0), [[0.0, 0.0, 0.0]], [11])


@pytest.fixture
def layer_cell() -> cc.Cell:
    """A 2D-periodic cell, which routes through the same analyzer."""
    return cc.Cell(
        cc.Lattice(np.diag([3.0, 3.0, 20.0])),
        [[0.0, 0.0, 0.5]],
        [6],
        cc.aperiodic_along(2),
    )


def space_hall(index: int) -> cc.HallNumber:
    return cc.HallNumber(cc.GroupFamily.space, index)


def layer_hall(index: int) -> cc.HallNumber:
    return cc.HallNumber(cc.GroupFamily.layer, index)
