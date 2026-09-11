"""Reciprocal-space sampling: mesh geometry, symmetry reduction, the zone."""

from __future__ import annotations

import numpy as np
import pytest

import pyseitz as sz


def test_a_mesh_needs_positive_divisions() -> None:
    """Both idioms, as everywhere else: `of` answers None, the ctor raises."""
    assert sz.Mesh.of([0, 4, 4]) is None
    with pytest.raises(ValueError):
        sz.Mesh([4, -2, 4])


def test_addresses_and_indices_round_trip() -> None:
    mesh = sz.Mesh([4, 4, 4])
    assert len(mesh) == 64
    assert mesh.index_of([1, 2, 3]) == 1 + 2 * 4 + 3 * 16
    assert mesh.index_of([-1, 0, 0]) == 3  # folds into the grid
    assert mesh.index_of(mesh.address_of(57)) == 57


def test_addresses_are_an_n_by_3_int_array() -> None:
    mesh = sz.Mesh([2, 3, 4])
    addresses = mesh.addresses

    assert addresses.shape == (len(mesh), 3)
    assert addresses.dtype == np.int32
    assert list(addresses[7]) == list(mesh.address_of(7))


def test_the_doubled_mesh_carries_the_shift() -> None:
    mesh = sz.Mesh([4, 4, 4], [True, False, False])
    assert mesh.doubled().divisions == [8, 8, 8]
    assert mesh.doubled_address([1, 0, 0]) == [3, 0, 0]


def test_a_cubic_mesh_reduces_to_its_irreducible_wedge(bcc_fe: sz.Cell) -> None:
    reciprocal = sz.analyze(bcc_fe).reciprocal_mesh(sz.Mesh([4, 4, 4]))

    assert len(reciprocal.rotations) == 48
    assert reciprocal.num_irreducible == len(set(reciprocal.mapping))
    assert reciprocal.num_irreducible < len(reciprocal.mesh)
    # A representative maps to itself, and every point maps to a representative.
    mapping = reciprocal.mapping
    assert all(mapping[m] == m for m in set(mapping))


def test_time_reversal_can_be_switched_off(bcc_fe: sz.Cell) -> None:
    analyzer = sz.analyze(bcc_fe)
    with_tr = analyzer.reciprocal_mesh(sz.Mesh([4, 4, 4]), sz.TimeReversal.on)
    without = analyzer.reciprocal_mesh(sz.Mesh([4, 4, 4]), sz.TimeReversal.off)

    assert without.num_irreducible >= with_tr.num_irreducible


def test_stabilizing_on_a_q_point_reduces_less(bcc_fe: sz.Cell) -> None:
    reciprocal = sz.analyze(bcc_fe).reciprocal_mesh(sz.Mesh([4, 4, 4]))
    stabilized = reciprocal.stabilized([np.array([0.25, 0.0, 0.0])])

    assert len(stabilized.rotations) <= len(reciprocal.rotations)
    assert stabilized.num_irreducible >= reciprocal.num_irreducible


def test_images_of_an_address_are_grid_points(bcc_fe: sz.Cell) -> None:
    reciprocal = sz.analyze(bcc_fe).reciprocal_mesh(sz.Mesh([4, 4, 4]))
    images = reciprocal.images_of([1, 0, 0])

    assert len(images) == len(reciprocal.rotations)
    assert all(0 <= i < len(reciprocal.mesh) for i in images)


def test_the_brillouin_zone_keeps_every_grid_point(bcc_fe: sz.Cell) -> None:
    """Boundary points are duplicated, so the zone holds at least the mesh."""
    analyzer = sz.analyze(bcc_fe)
    reciprocal = analyzer.reciprocal_mesh(sz.Mesh([4, 4, 4]))
    zone = reciprocal.brillouin_zone(sz.Lattice(np.linalg.inv(bcc_fe.lattice.matrix).T))

    assert zone.addresses.shape[0] >= len(reciprocal.mesh)
    assert zone.addresses.dtype == np.int32
    assert zone.map(0) is not None
    assert len(zone.images_of([1, 0, 0])) == len(reciprocal.rotations)
