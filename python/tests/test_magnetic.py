"""The magnetic layer: moments in, a UNI number out."""

from __future__ import annotations

import numpy as np
import pytest

import pyseitz as sz


@pytest.fixture
def collinear(bcc_fe: sz.Cell) -> sz.MagneticCell:
    """Antiferromagnetic bcc iron: one moment up, one down."""
    return sz.MagneticCell(bcc_fe, np.array([1.0, -1.0]), sz.TensorKind.axial)


@pytest.fixture
def noncollinear(bcc_fe: sz.Cell) -> sz.MagneticCell:
    moments = np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
    return sz.MagneticCell(bcc_fe, moments, sz.TensorKind.axial)


def test_the_rank_is_the_shape_of_the_moments(collinear, noncollinear) -> None:
    """No rank flag anywhere: the array's shape decides, once."""
    assert collinear.rank == sz.SiteTensor.collinear
    assert noncollinear.rank == sz.SiteTensor.noncollinear
    assert collinear.moments.shape == (2,)
    assert noncollinear.moments.shape == (2, 3)


def test_moments_of_the_wrong_shape_are_refused(bcc_fe: sz.Cell) -> None:
    with pytest.raises(ValueError):
        sz.MagneticCell(bcc_fe, np.array([1.0, -1.0, 1.0]))
    with pytest.raises(ValueError):
        sz.MagneticCell(bcc_fe, np.zeros((2, 2)))


def test_a_collinear_antiferromagnet_determines(collinear: sz.MagneticCell) -> None:
    analyzer = sz.analyze_magnetic(collinear)

    assert 1 <= analyzer.uni.value <= 1651
    assert analyzer.hall.family == sz.GroupFamily.space
    assert len(analyzer.operations) > 0
    assert len(analyzer.equivalent_atoms) == len(collinear)
    assert analyzer.spacegroup_type.uni_number == analyzer.uni.value
    assert analyzer.spacegroup_type.bns_number


def test_a_noncollinear_structure_determines(noncollinear: sz.MagneticCell) -> None:
    assert 1 <= sz.analyze_magnetic(noncollinear).uni.value <= 1651


def test_the_dataset_carries_the_match_and_the_cell(collinear: sz.MagneticCell) -> None:
    dataset = sz.analyze_magnetic(collinear).dataset

    assert dataset.uni.value >= 1
    assert dataset.type in tuple(sz.MagneticType)
    assert dataset.setting.transformation.shape == (3, 3)
    assert len(dataset.standardized) >= len(collinear)


def test_the_operations_identify_the_same_group(collinear: sz.MagneticCell) -> None:
    """The public door from a bare operation set, as the non-magnetic one has."""
    analyzer = sz.analyze_magnetic(collinear)
    match = analyzer.operations.spacegroup(collinear.cell.lattice)

    assert match.uni == analyzer.uni
    assert match.hall == analyzer.hall


def test_the_magnetic_database_answers_by_uni_number() -> None:
    uni = sz.UniNumber(1651)
    metadata = sz.magnetic_spacegroup_type(uni)

    assert metadata.uni_number == 1651
    assert metadata.number == 230
    assert 1 <= metadata.type <= 4
    assert len(sz.magnetic_operations_from_database(uni)) > 0
    assert len(sz.magnetic_std_transformations(uni)) > 0


def test_uni_candidates_bracket_a_hall_setting() -> None:
    hall = sz.default_hall(sz.GroupFamily.space, 225)
    assert hall is not None
    first, last = sz.uni_candidates(hall)
    assert first <= last


def test_a_magnetic_analyzer_refuses_to_be_copied(collinear: sz.MagneticCell) -> None:
    import copy

    with pytest.raises(TypeError, match="not copyable"):
        copy.deepcopy(sz.analyze_magnetic(collinear))


def test_the_record_round_trips_through_json(collinear: sz.MagneticCell) -> None:
    record = sz.MagneticDatasetRecord.from_analyzer(sz.analyze_magnetic(collinear))
    assert record == sz.MagneticDatasetRecord.model_validate_json(record.model_dump_json())
