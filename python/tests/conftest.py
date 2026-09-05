"""Fixtures mirroring tests/helpers.hpp, so the two suites test the same cells."""

from __future__ import annotations

import os
import pathlib

import numpy as np
import pytest

import seitz as sz


def pytest_configure(config: pytest.Config) -> None:
    """Fail loudly if we imported a different _core than the one under test.

    CTest points PYTHONPATH at the build tree it just compiled -- which for the
    sanitizer preset is the whole point of test_lifetimes.py.  But PYTHONPATH is
    only consulted by the path finder, and anything that installs a meta-path
    finder for this package (a scikit-build-core editable install is the one
    that bit us) overrides it silently: the suite passes, having exercised some
    other build entirely.  Compare the two and say so instead.
    """
    expected = os.environ.get("PYTHONPATH", "").split(os.pathsep)[0]
    if not expected:
        return
    loaded = pathlib.Path(sz._core.__file__).resolve()
    if not loaded.is_relative_to(pathlib.Path(expected).resolve()):
        raise pytest.UsageError(
            f"seitz._core was imported from {loaded}, but PYTHONPATH points "
            f"at {expected}. Something is shadowing the build tree -- most "
            f"likely an editable install of seitz in this interpreter. "
            f"Remove it (`uv pip uninstall seitz`); this project is reached "
            f"through the CMake build tree, not an editable install."
        )


@pytest.fixture
def bcc_fe() -> sz.Cell:
    """The demo cell from main.cpp: 2 Fe in a 3 A cube -> Im-3m, No. 229."""
    return sz.Cell(
        sz.Lattice(3.0 * np.eye(3)),
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        [26, 26],
    )


@pytest.fixture
def simple_cubic() -> sz.Cell:
    """One atom in a cube -> Pm-3m, No. 221."""
    return sz.Cell(sz.Lattice(np.eye(3) * 4.0), [[0.0, 0.0, 0.0]], [11])


@pytest.fixture
def layer_cell() -> sz.Cell:
    """A 2D-periodic cell, which routes through the same analyzer."""
    return sz.Cell(
        sz.Lattice(np.diag([3.0, 3.0, 20.0])),
        [[0.0, 0.0, 0.5]],
        [6],
        sz.aperiodic_along(2),
    )


def space_hall(index: int) -> sz.HallNumber:
    return sz.HallNumber(sz.GroupFamily.space, index)


def layer_hall(index: int) -> sz.HallNumber:
    return sz.HallNumber(sz.GroupFamily.layer, index)
