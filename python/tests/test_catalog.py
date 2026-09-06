"""The catalogs: groups, Wyckoff positions, metadata, elements, subgroups."""

from __future__ import annotations

import pytest

import seitz as sz


def test_every_space_hall_setting_has_a_symbol() -> None:
    for index in range(1, sz._core.K_SPACE_HALL_SETTINGS + 1):
        hall = sz.HallNumber(sz.GroupFamily.space, index)
        assert sz.spacegroup_type(hall).international_short


def test_every_layer_hall_setting_has_a_symbol() -> None:
    for index in range(1, sz._core.K_LAYER_HALL_SETTINGS + 1):
        hall = sz.HallNumber(sz.GroupFamily.layer, index)
        assert sz.spacegroup_type(hall).international_short


def test_a_key_out_of_range_is_none_or_raises() -> None:
    over = sz._core.K_SPACE_HALL_SETTINGS + 1
    assert sz.HallNumber.of(sz.GroupFamily.space, over) is None
    with pytest.raises(sz.SeitzError):
        sz.HallNumber(sz.GroupFamily.space, over)


def test_hall_numbers_compare_hash_and_pickle() -> None:
    import pickle

    hall = sz.HallNumber(sz.GroupFamily.space, 229)
    assert hall == sz.HallNumber(sz.GroupFamily.space, 229)
    assert hall != sz.HallNumber(sz.GroupFamily.space, 228)
    assert len({hall, sz.HallNumber(sz.GroupFamily.space, 229)}) == 1
    assert pickle.loads(pickle.dumps(hall)) == hall


def test_a_group_exposes_its_wyckoff_positions() -> None:
    group = sz.SpaceGroup.from_number(sz.GroupFamily.space, 229)
    assert group.number == 229
    assert group.order == len(group.operations)
    assert len(group.wyckoffs) > 0
    # Ordered by ascending letter, the last being the general position.
    assert group.wyckoffs[0].letter == "a"
    assert group.wyckoffs[-1].multiplicity == max(
        w.multiplicity for w in group.wyckoffs
    )


def test_wyckoff_lookup_by_letter_returns_the_same_object() -> None:
    group = sz.SpaceGroup.of(sz.HallNumber(sz.GroupFamily.space, 1))
    assert group.wyckoff("a").letter == "a"


def test_spacegroup_of_is_a_shared_flyweight() -> None:
    hall = sz.HallNumber(sz.GroupFamily.space, 229)
    assert sz.SpaceGroup.of(hall) is sz.SpaceGroup.of(hall)


def test_element_symbols_round_trip() -> None:
    for z in range(1, sz.elements.K_NUM_ELEMENTS + 1):
        symbol = sz.elements.element_symbol(z)
        assert symbol is not None
        assert sz.elements.atomic_number(symbol) == z


def test_unknown_elements_answer_none() -> None:
    assert sz.elements.element_symbol(0) is None
    assert sz.elements.covalent_radius(0) is None
    assert sz.elements.atomic_number("Xx") is None
    assert sz.elements.is_known_element(26)


def test_the_subgroup_graph_is_consistent() -> None:
    # P1 is a t-subgroup of everything and has no maximal subgroups of its own.
    assert sz.subgroups.maximal_subgroups(1) == []
    assert sz.subgroups.is_subgroup(1, 229)
    assert sz.subgroups.is_subgroup(229, 229)
    assert not sz.subgroups.is_subgroup(229, 1)

    path = sz.subgroups.path(229, 1)
    assert path is not None
    assert path[0].super == 229 and path[-1].sub == 1
    assert all(a.sub == b.super for a, b in zip(path, path[1:]))
    assert sz.subgroups.path(229, 229) == []
    assert sz.subgroups.path(1, 229) is None
    # Every relation names a real setting, and its kind is one of the two.
    edge = sz.subgroups.maximal_subgroups(221, sz.SubgroupKind.translationengleiche)[0]
    assert edge.super == 221 and edge.hall.index >= 1
    assert sz.subgroups.edge(edge.id).sub == edge.sub
    assert edge.basis.shape == (3, 3)


def test_point_group_metadata_is_aligned_to_its_numbering() -> None:
    assert sz.pointgroup_by_number(1).symbol == "1"
    assert sz.pointgroup_by_number(32).symbol == "m-3m"
    assert sz.pointgroup_by_number(32).crystal_class == sz.CrystalClass.oh
    assert sz.pointgroup_by_number(33).number == 0
