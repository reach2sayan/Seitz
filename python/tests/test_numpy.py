"""The NumPy contract, stated as assertions so it cannot drift silently."""

from __future__ import annotations

import numpy as np

import seitz as sz


def test_positions_are_c_contiguous_float64(bcc_fe: sz.Cell) -> None:
    positions = bcc_fe.positions
    assert positions.dtype == np.float64
    assert positions.shape == (2, 3)
    assert positions.flags.c_contiguous


def test_types_are_an_int32_array(bcc_fe: sz.Cell) -> None:
    types = bcc_fe.types
    assert types.dtype == np.int32
    assert types.shape == (2,)
    assert np.array_equal(types, [26, 26])


def test_the_basis_columns_are_the_basis_vectors() -> None:
    """Matrix3d is column-major, so this holds with no transpose anywhere."""
    basis = np.array([[1.0, 0.5, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]])
    lattice = sz.Lattice(basis)

    assert lattice.matrix.flags.f_contiguous
    assert np.allclose(lattice.matrix[:, 0], lattice.a)
    assert np.allclose(lattice.matrix[:, 1], lattice.b)
    assert np.allclose(lattice.matrix[:, 2], lattice.c)
    assert np.allclose(lattice.matrix, basis)


def test_the_matrix_property_hands_back_a_copy() -> None:
    """A Lattice is immutable in C++; Python must not get a writable alias."""
    lattice = sz.Lattice(np.eye(3) * 2.0)
    first = lattice.matrix
    first[0, 0] = 99.0
    assert lattice.matrix[0, 0] == 2.0


def test_positions_accept_a_plain_nested_list(bcc_fe: sz.Cell) -> None:
    from_list = sz.Cell(sz.Lattice(np.eye(3)), [[0.0, 0.0, 0.0]], [1])
    assert from_list.positions.shape == (1, 3)


def test_an_orbit_is_an_n_by_3_array() -> None:
    group = sz.SpaceGroup.of(sz.HallNumber(sz.GroupFamily.space, 1))
    general = group.wyckoffs[-1]
    orbit = general.orbit(np.array([0.1, 0.2, 0.3]))
    assert orbit.shape == (general.multiplicity, 3)
    assert orbit.dtype == np.float64
