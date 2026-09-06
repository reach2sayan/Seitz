"""One adapter per package. ``prepare`` runs outside the timed region; a missing
workload is ``None`` and the case is skipped. Missing packages yield no backend."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from structures import Structure

SYMPREC = 1e-3  # the relaxed-structure regime; the jitter in structures.py sits below it


@dataclass(frozen=True)
class Backend:
    prepare: Callable[[Structure], Any]
    number: Callable[[Any], int] | None = None  # space-group number of a dataset result
    dataset: Callable[[Any], Any] | None = None
    standardize: Callable[[Any], Any] | None = None
    primitive: Callable[[Any], Any] | None = None
    read_cif: Callable[[Path], Any] | None = None
    max_atoms: int = 1 << 30
    note: str = ""


def _seitz() -> Backend:
    import seitz as sz
    tol = sz.Tolerance(symprec=SYMPREC)
    sz.warmup()
    return Backend(
        prepare=lambda s: sz.Cell(sz.Lattice(s.basis.T), s.positions, s.numbers.tolist()),
        number=lambda d: sz.spacegroup_type(d.hall).number,
        dataset=lambda c: sz.analyze(c, tol).dataset,
        standardize=lambda c: sz.analyze(c, tol).standardized_cell,
        # Not .primitive_cell: that is a pure-translation search with no
        # space-group determination, so it is not what find_primitive does.
        primitive=lambda c: sz.analyze(c, tol).standardized_cell_in(sz.CellSetting.primitive),
        read_cif=sz.read_cif,
    )


def _spglib() -> Backend:
    import spglib
    return Backend(
        prepare=lambda s: (s.basis, s.positions, s.numbers),
        number=lambda d: d.number,
        dataset=lambda c: spglib.get_symmetry_dataset(c, symprec=SYMPREC),
        standardize=lambda c: spglib.standardize_cell(c, symprec=SYMPREC),
        primitive=lambda c: spglib.find_primitive(c, symprec=SYMPREC),
    )


def _moyopy() -> Backend:
    import moyopy
    run = lambda c: moyopy.MoyoDataset(c, symprec=SYMPREC, angle_tolerance=None)
    return Backend(
        prepare=lambda s: moyopy.Cell(s.basis.tolist(), s.positions.tolist(), s.numbers.tolist()),
        number=lambda d: d.number,
        dataset=run,
        standardize=lambda c: run(c).std_cell,
        primitive=lambda c: run(c).prim_std_cell,
        note="one call computes dataset, standardized and primitive cells together",
    )


def _pymatgen_structure(s: Structure):
    from pymatgen.core import Lattice, Structure as PmgStructure
    return PmgStructure(Lattice(s.basis), s.numbers, s.positions)


def _pymatgen() -> Backend:
    from pymatgen.io.cif import CifParser
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    sga = lambda st: SpacegroupAnalyzer(st, symprec=SYMPREC, angle_tolerance=-1)  # derived, like the rest
    return Backend(
        prepare=_pymatgen_structure,
        number=lambda d: d.number,
        dataset=lambda st: sga(st).get_symmetry_dataset(),
        standardize=lambda st: sga(st).get_refined_structure(),
        primitive=lambda st: sga(st).find_primitive(),
        read_cif=lambda p: CifParser(p).parse_structures(primitive=False),
    )


def _ase() -> Backend:
    import ase.io
    from ase import Atoms
    from ase.spacegroup.symmetrize import check_symmetry
    return Backend(
        prepare=lambda s: Atoms(numbers=s.numbers, scaled_positions=s.positions, cell=s.basis, pbc=True),
        number=lambda d: d.number,
        dataset=lambda a: check_symmetry(a, symprec=SYMPREC),
        read_cif=ase.io.read,
    )


def _pyxtal() -> Backend:
    from pyxtal import pyxtal
    def dataset(st):
        c = pyxtal(); c.from_seed(st, tol=SYMPREC); return c
    def read_cif(p):
        c = pyxtal(); c.from_seed(str(p)); return c
    return Backend(
        prepare=_pymatgen_structure,
        number=lambda c: c.group.number,
        dataset=dataset,
        read_cif=read_cif,
        max_atoms=432,
        note="from_seed also assigns Wyckoff positions",
    )


def _gemmi() -> Backend:
    import gemmi
    return Backend(
        prepare=lambda s: None,
        read_cif=lambda p: gemmi.read_small_structure(str(p)).get_all_unit_cell_sites(),
        note="parse and expand only, no symmetry search",
    )


def _load(make: Callable[[], Backend]) -> Backend | None:
    try: return make()
    except ImportError: return None


BACKENDS: dict[str, Backend] = {
    name: b for name, make in {
        "seitz": _seitz, "spglib": _spglib, "moyopy": _moyopy,
        "pymatgen": _pymatgen, "ase": _ase, "pyxtal": _pyxtal, "gemmi": _gemmi,
    }.items() if (b := _load(make)) is not None
}
