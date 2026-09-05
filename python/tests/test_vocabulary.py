"""The value types: cells, lattices, operations, periodicity, pydantic models."""

from __future__ import annotations

import copy
import pickle

import numpy as np
import pytest
from pydantic import ValidationError

import seitz as sz


def test_a_cell_is_a_sequence_of_atoms(bcc_fe: sz.Cell) -> None:
    assert len(bcc_fe) == 2
    position, atom_type = bcc_fe[1]
    assert atom_type == 26
    assert np.allclose(position, [0.5, 0.5, 0.5])
    assert np.allclose(bcc_fe[-1][0], [0.5, 0.5, 0.5])

    with pytest.raises(IndexError):
        bcc_fe[2]


def test_cells_and_lattices_round_trip_through_pickle(bcc_fe: sz.Cell) -> None:
    restored = pickle.loads(pickle.dumps(bcc_fe))
    assert len(restored) == len(bcc_fe)
    assert np.allclose(restored.positions, bcc_fe.positions)
    assert np.allclose(copy.deepcopy(bcc_fe).positions, bcc_fe.positions)


def test_lattice_geometry() -> None:
    lattice = sz.Lattice(np.diag([2.0, 3.0, 4.0]))
    assert lattice.volume == pytest.approx(24.0)
    assert np.allclose(lattice.metric, np.diag([4.0, 9.0, 16.0]))
    assert np.allclose(lattice.to_cartesian([1.0, 0.0, 0.0]), [2.0, 0.0, 0.0])
    assert np.allclose(lattice.to_fractional([2.0, 0.0, 0.0]), [1.0, 0.0, 0.0])


def test_lattice_reduction_returns_a_lattice() -> None:
    skewed = sz.Lattice(np.array([[4.0, 3.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]]))
    assert skewed.niggli().volume == pytest.approx(skewed.volume)
    assert skewed.delaunay().volume == pytest.approx(skewed.volume)


def test_periodicity_helpers_describe_the_families() -> None:
    assert sz.family_of(sz.all_periodic()) == sz.GroupFamily.space
    assert sz.family_of(sz.aperiodic_along(2)) == sz.GroupFamily.layer
    assert sz.aperiodic_axis(sz.aperiodic_along(2)) == 2
    assert sz.aperiodic_axis(sz.all_periodic()) is None
    assert sz.aperiodic_axis(sz.none_periodic()) is None


def test_symmetry_operations_compose_and_invert() -> None:
    op = sz.SymmetryOperation(np.eye(3, dtype=np.int32), [0.5, 0.0, 0.0])
    assert np.allclose(op.apply([0.1, 0.0, 0.0]), [0.6, 0.0, 0.0])
    assert op.is_identity_rotation

    inverse = op.inverse
    assert inverse is not None
    assert np.allclose((op * inverse).translation, [0.0, 0.0, 0.0])


def test_operations_compare_through_same_operation_not_eq() -> None:
    """SymmetryOperation has no __eq__ on purpose -- see the docstring."""
    a = sz.SymmetryOperation(np.eye(3, dtype=np.int32), [0.0, 0.0, 0.0])
    b = sz.SymmetryOperation(np.eye(3, dtype=np.int32), [1e-9, 0.0, 0.0])
    assert sz.same_operation(a, b, 1e-5)
    assert not sz.same_operation(a, b, 1e-12)


def test_tolerance_validates_its_inputs() -> None:
    assert sz.Tolerance().symprec == sz.K_DEFAULT_SYMPREC
    assert sz.Tolerance(symprec=1e-3).angle_tolerance is None

    with pytest.raises(ValidationError):
        sz.Tolerance(symprec=-1.0)
    with pytest.raises(ValidationError):
        sz.Tolerance(symprec=1e-5, nonsense=1)


def test_tolerance_round_trips_through_json() -> None:
    tolerance = sz.Tolerance(symprec=1e-4, angle_tolerance=5.0)
    restored = sz.Tolerance.model_validate_json(tolerance.model_dump_json())
    assert restored == tolerance


def test_cell_records_round_trip(bcc_fe: sz.Cell) -> None:
    record = sz.CellRecord.from_cell(bcc_fe)
    restored = sz.CellRecord.model_validate_json(record.model_dump_json())
    assert restored == record

    cell = restored.to_cell()
    assert len(cell) == len(bcc_fe)
    assert np.allclose(cell.positions, bcc_fe.positions)
    assert np.allclose(cell.lattice.matrix, bcc_fe.lattice.matrix)


def test_dataset_records_serialize(bcc_fe: sz.Cell) -> None:
    record = sz.DatasetRecord.from_analyzer(sz.analyze(bcc_fe))
    assert record.number == 229
    assert record.international_short == "Im-3m"
    assert record.num_operations == 96
    assert len(record.sites) == 2

    restored = sz.DatasetRecord.model_validate_json(record.model_dump_json())
    assert restored.number == record.number
    assert restored.sites == record.sites


def test_positions_validation_rejects_the_wrong_shape() -> None:
    with pytest.raises(ValidationError):
        sz.CellRecord(
            basis=np.eye(3), positions=[[0.0, 0.0]], types=[1],
            periodicity=(0, 0, 0),
        )


def test_the_version_is_reported() -> None:
    assert sz.version_string().count(".") == 2
    assert sz.__version__ == sz.version_string()
