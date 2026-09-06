"""The CIF layer: reading, writing, and the round trip between them."""

from __future__ import annotations

import numpy as np
import pytest

import seitz


NACL = """data_nacl
_cell_length_a 5.64
_cell_length_b 5.64
_cell_length_c 5.64
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_space_group_name_H-M_alt 'F m -3 m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Na1 Na 0.0 0.0 0.0
Cl1 Cl 0.5 0.5 0.5
"""


def test_parse_cif_is_uninterpreted() -> None:
    blocks = seitz.parse_cif(NACL)
    assert len(blocks) == 1
    block = blocks[0]
    assert block.name == "nacl"
    assert block.value("_cell_length_a") == "5.64"
    assert block.column("_atom_site_label") == ["Na1", "Cl1"]
    assert block.value("_atom_site_label") == "Na1"
    assert block.column("_nope") is None
    assert block.loops == [
        [
            "_atom_site_label",
            "_atom_site_type_symbol",
            "_atom_site_fract_x",
            "_atom_site_fract_y",
            "_atom_site_fract_z",
        ]
    ]
    assert block.columns["_cell_length_b"] == ["5.64"]
    assert "CifBlock(name=nacl" in repr(block)


def test_read_cif_expands_the_asymmetric_unit() -> None:
    (structure,) = seitz.read_cif(NACL)
    assert structure.name == "nacl"
    assert len(structure.cell.types) == 8
    assert structure.hall is not None
    assert seitz.spacegroup_type(structure.hall).number == 225
    assert structure.collapsed == []
    assert set(structure.labels) == {"Na1", "Cl1"}


def test_read_cif_takes_a_path(tmp_path) -> None:
    path = tmp_path / "nacl.cif"
    path.write_text(NACL)
    from_path = seitz.read_cif(path)
    from_text = seitz.read_cif(NACL)
    assert len(from_path) == len(from_text) == 1
    assert np.allclose(from_path[0].cell.positions, from_text[0].cell.positions)


def test_read_cif_defaults_to_the_cif_tolerance() -> None:
    # Not the search default: five-decimal coordinates need the looser one.
    assert seitz.CIF_SYMPREC == pytest.approx(1e-3)
    assert seitz.CIF_SYMPREC > seitz.K_DEFAULT_SYMPREC
    # An explicit tolerance is accepted as a model or a plain dict.
    assert len(seitz.read_cif(NACL, {"symprec": 1e-3})) == 1
    assert len(seitz.read_cif(NACL, seitz.Tolerance(symprec=1e-3))) == 1


def test_a_shared_site_is_reported_not_hidden() -> None:
    text = NACL.replace(
        "_atom_site_fract_z\nNa1 Na 0.0 0.0 0.0",
        "_atom_site_fract_z\n_atom_site_occupancy\nNa1 Na 0.0 0.0 0.0 0.6\n"
        "K1 K 0.0 0.0 0.0 0.4",
    ).replace("Cl1 Cl 0.5 0.5 0.5\n", "Cl1 Cl 0.5 0.5 0.5 1.0\n")
    (structure,) = seitz.read_cif(text)
    assert len(structure.collapsed) == 1
    collapse = structure.collapsed[0]
    assert collapse.kept == "Na1"
    assert collapse.dropped == ["K1"]
    assert collapse.occupancy == pytest.approx(0.6)
    assert "OccupancyCollapse(kept=Na1" in repr(collapse)


def test_write_cif_p1_round_trips(bcc_fe) -> None:
    text = seitz.write_cif(bcc_fe, name="written")
    (back,) = seitz.read_cif(text)
    assert back.name == "written"
    assert len(back.cell.types) == len(bcc_fe.types)
    assert np.allclose(back.cell.lattice.matrix, bcc_fe.lattice.matrix, atol=1e-5)


def test_write_cif_symmetrized_names_its_group(bcc_fe) -> None:
    analyzer = seitz.analyze(bcc_fe)
    text = seitz.write_cif(analyzer, name="determined")
    assert "_space_group_symop_operation_xyz" in text
    assert "_atom_site_Wyckoff_symbol" in text
    (back,) = seitz.read_cif(text)
    assert back.hall is not None
    assert (
        seitz.spacegroup_type(back.hall).number
        == analyzer.spacegroup_type.number
    )


def test_write_cif_writes_a_file(bcc_fe, tmp_path) -> None:
    path = tmp_path / "out.cif"
    text = seitz.write_cif(bcc_fe, path=path)
    assert path.read_text() == text


def test_write_cif_rejects_anything_else() -> None:
    with pytest.raises(TypeError, match="Cell or a SymmetryAnalyzer"):
        seitz.write_cif("not a cell")


def test_symmetry_operation_xyz_round_trips() -> None:
    for op in seitz.operations_from_database(seitz.default_hall(seitz.GroupFamily.space, 227)):
        back = seitz.SymmetryOperation.from_xyz(op.to_xyz())
        assert seitz.same_operation(back, op, 1e-9)


def test_symmetry_operation_from_xyz_reads_what_files_write() -> None:
    op = seitz.SymmetryOperation.from_xyz("1/2+x, -y, z")
    assert np.array_equal(op.rotation, np.diag([1, -1, 1]))
    assert np.allclose(op.translation, [0.5, 0.0, 0.0])
    assert seitz.SymmetryOperation.from_xyz("x,y,z").to_xyz() == "x,y,z"
