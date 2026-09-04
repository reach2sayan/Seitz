"""The catalogs: groups, Wyckoff positions, metadata, elements, subgroups."""

from __future__ import annotations

import pytest

import cppcrystal as cc


def test_every_space_hall_setting_has_a_symbol() -> None:
    for index in range(1, cc._core.K_SPACE_HALL_SETTINGS + 1):
        hall = cc.HallNumber(cc.GroupFamily.space, index)
        assert cc.spacegroup_type(hall).international_short


def test_every_layer_hall_setting_has_a_symbol() -> None:
    for index in range(1, cc._core.K_LAYER_HALL_SETTINGS + 1):
        hall = cc.HallNumber(cc.GroupFamily.layer, index)
        assert cc.spacegroup_type(hall).international_short


def test_a_key_out_of_range_is_none_or_raises() -> None:
    over = cc._core.K_SPACE_HALL_SETTINGS + 1
    assert cc.HallNumber.of(cc.GroupFamily.space, over) is None
    with pytest.raises(cc.CppCrystalError):
        cc.HallNumber(cc.GroupFamily.space, over)


def test_hall_numbers_compare_hash_and_pickle() -> None:
    import pickle

    hall = cc.HallNumber(cc.GroupFamily.space, 229)
    assert hall == cc.HallNumber(cc.GroupFamily.space, 229)
    assert hall != cc.HallNumber(cc.GroupFamily.space, 228)
    assert len({hall, cc.HallNumber(cc.GroupFamily.space, 229)}) == 1
    assert pickle.loads(pickle.dumps(hall)) == hall


def test_a_group_exposes_its_wyckoff_positions() -> None:
    group = cc.SpaceGroup.from_number(cc.GroupFamily.space, 229)
    assert group.number == 229
    assert group.order == len(group.operations)
    assert len(group.wyckoffs) > 0
    # Ordered by ascending letter, the last being the general position.
    assert group.wyckoffs[0].letter == "a"
    assert group.wyckoffs[-1].multiplicity == max(
        w.multiplicity for w in group.wyckoffs
    )


def test_wyckoff_lookup_by_letter_returns_the_same_object() -> None:
    group = cc.SpaceGroup.of(cc.HallNumber(cc.GroupFamily.space, 1))
    assert group.wyckoff("a").letter == "a"


def test_spacegroup_of_is_a_shared_flyweight() -> None:
    hall = cc.HallNumber(cc.GroupFamily.space, 229)
    assert cc.SpaceGroup.of(hall) is cc.SpaceGroup.of(hall)


def test_element_symbols_round_trip() -> None:
    for z in range(1, cc.elements.K_NUM_ELEMENTS + 1):
        symbol = cc.elements.element_symbol(z)
        assert symbol is not None
        assert cc.elements.atomic_number(symbol) == z


def test_unknown_elements_answer_none() -> None:
    assert cc.elements.element_symbol(0) is None
    assert cc.elements.covalent_radius(0) is None
    assert cc.elements.atomic_number("Xx") is None
    assert cc.elements.is_known_element(26)


def test_the_subgroup_graph_is_consistent() -> None:
    # P1 is a t-subgroup of everything and has no maximal subgroups of its own.
    assert cc.subgroups.maximal_subgroups(1) == []
    assert cc.subgroups.is_subgroup(1, 229)
    assert cc.subgroups.is_subgroup(229, 229)
    assert not cc.subgroups.is_subgroup(229, 1)

    path = cc.subgroups.path(229, 1)
    assert path is not None
    assert path[0] == 229 and path[-1] == 1
    assert cc.subgroups.path(1, 229) is None


def test_point_group_metadata_is_aligned_to_its_numbering() -> None:
    assert cc.pointgroup_by_number(1).symbol == "1"
    assert cc.pointgroup_by_number(32).symbol == "m-3m"
    assert cc.pointgroup_by_number(32).crystal_class == cc.CrystalClass.oh
    assert cc.pointgroup_by_number(33).number == 0
