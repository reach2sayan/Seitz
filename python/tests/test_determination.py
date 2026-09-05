"""The determination path, end to end."""

from __future__ import annotations

import numpy as np
import pytest

import seitz as sz


def test_bcc_iron_is_im3m(bcc_fe: sz.Cell) -> None:
    """The round trip from main.cpp, which is the canonical smoke test."""
    analyzer = sz.analyze(bcc_fe)
    spacegroup_type = analyzer.spacegroup_type

    assert spacegroup_type.number == 229
    assert spacegroup_type.international_short == "Im-3m"
    assert len(analyzer.operations) == 96


def test_the_analyzer_memoizes_rather_than_recomputing(bcc_fe: sz.Cell) -> None:
    analyzer = sz.analyze(bcc_fe)
    assert analyzer.hall == analyzer.hall
    assert analyzer.dataset.hall == analyzer.hall


def test_simple_cubic_is_pm3m(simple_cubic: sz.Cell) -> None:
    assert sz.analyze(simple_cubic).spacegroup_type.number == 221


def test_a_layer_cell_uses_the_same_analyzer(layer_cell: sz.Cell) -> None:
    """Layer groups are not a separate entry point -- the periodicity decides."""
    analyzer = sz.analyze(layer_cell)
    assert analyzer.hall.family == sz.GroupFamily.layer
    assert 1 <= analyzer.spacegroup_type.number <= 80


def test_tolerance_accepts_a_model_a_dict_or_nothing(bcc_fe: sz.Cell) -> None:
    loose = sz.Tolerance(symprec=1e-3)
    assert sz.analyze(bcc_fe, loose).spacegroup_type.number == 229
    assert sz.analyze(bcc_fe, {"symprec": 1e-3}).spacegroup_type.number == 229
    assert sz.analyze(bcc_fe).spacegroup_type.number == 229


def test_a_fixed_setting_is_honoured(bcc_fe: sz.Cell) -> None:
    hall = sz.default_hall(sz.GroupFamily.space, 229)
    assert hall is not None
    assert sz.analyze(bcc_fe, setting=hall).hall == hall


def test_sites_and_site_arrays_agree(bcc_fe: sz.Cell) -> None:
    analyzer = sz.analyze(bcc_fe)
    sites = analyzer.sites
    arrays = analyzer.site_arrays()

    assert len(sites) == len(bcc_fe)
    assert np.array_equal(arrays["wyckoff"], [s.wyckoff for s in sites])
    assert arrays["site_symmetry"] == [s.site_symmetry for s in sites]


def test_standardized_cell_in_each_setting(bcc_fe: sz.Cell) -> None:
    analyzer = sz.analyze(bcc_fe)
    conventional = analyzer.standardized_cell_in(sz.CellSetting.conventional)
    primitive = analyzer.standardized_cell_in(sz.CellSetting.primitive)

    # Im-3m is body-centred: the primitive cell holds half the atoms.
    assert len(conventional) == 2 * len(primitive)
    assert len(analyzer.standardized_cell) == len(conventional)


def test_an_analyzer_refuses_to_be_copied(bcc_fe: sz.Cell) -> None:
    """It owns a memo guarded by a mutex; sharing it is the supported move."""
    import copy

    analyzer = sz.analyze(bcc_fe)
    with pytest.raises(TypeError, match="not copyable"):
        copy.deepcopy(analyzer)
