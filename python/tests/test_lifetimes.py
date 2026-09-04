"""Borrowed references and keep_alive, which are what the bindings promise.

Every test here would segfault, not merely fail, if the corresponding lifetime
annotation went missing -- so this file is worth running under
-DCPPCRYSTAL_ENABLE_SANITIZERS=ON at least once per change to the bindings.
"""

from __future__ import annotations

import gc

import numpy as np

import cppcrystal as cc


def test_a_wyckoff_outlives_the_name_that_reached_it() -> None:
    """wyckoffs borrows into the group's storage, parented on the group."""
    group = cc.SpaceGroup.from_number(cc.GroupFamily.space, 225)
    position = group.wyckoffs[0]
    del group
    gc.collect()
    assert position.multiplicity > 0
    assert position.letter == "a"


def test_a_wyckoff_from_lookup_also_keeps_its_group_alive() -> None:
    group = cc.SpaceGroup.from_number(cc.GroupFamily.space, 225)
    position = group.wyckoff("a")
    del group
    gc.collect()
    assert position.site_symmetry is not None


def test_cell_atoms_survives_the_cell(bcc_fe: cc.Cell) -> None:
    """Cell.atoms() is a lazy view in C++ and is materialised for exactly this."""
    atoms = bcc_fe.atoms
    listed = list(bcc_fe)
    del bcc_fe
    gc.collect()
    assert len(atoms) == 2
    assert len(listed) == 2
    assert np.allclose(atoms[0][0], [0.0, 0.0, 0.0])


def test_a_dataset_outlives_its_analyzer(bcc_fe: cc.Cell) -> None:
    """Projections copy out of the memo rather than borrowing from it."""
    analyzer = cc.analyze(bcc_fe)
    dataset = analyzer.dataset
    operations = analyzer.operations
    del analyzer
    gc.collect()
    assert len(operations) == 96
    assert dataset.hall.index > 0


def test_operations_iterate_safely(bcc_fe: cc.Cell) -> None:
    operations = cc.analyze(bcc_fe).operations
    gc.collect()
    assert len(list(operations)) == len(operations)
    assert operations[0].rotation.shape == (3, 3)
