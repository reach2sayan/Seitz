"""The cells under test, backend-neutral, mirroring tests/helpers.hpp."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from pathlib import Path

import numpy as np

NA, CL = 11, 17


@dataclass(frozen=True)
class Structure:
    """``basis`` rows are the lattice vectors (the spglib convention)."""

    name: str
    basis: np.ndarray
    positions: np.ndarray
    numbers: np.ndarray

    def __len__(self) -> int: return len(self.numbers)


def rocksalt_supercell(k: int, jitter: float = 1e-4, seed: int = 7) -> Structure:
    """k x k x k NaCl in a 4 A cube, each coordinate jittered by up to ``jitter`` angstrom.

    Below SYMPREC, so the search still resolves Fm-3m but pays for real numerics."""
    rng = np.random.default_rng(seed)
    basis = np.eye(3) * 4.0 * k
    offsets = np.indices((k, k, k)).reshape(3, -1).T / k
    motif = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]) / k
    positions = (offsets[:, None, :] + motif).reshape(-1, 3)
    positions += rng.uniform(-jitter, jitter, positions.shape) @ np.linalg.inv(basis)
    numbers = np.tile([NA, CL], k**3)
    return Structure(f"rocksalt-{len(numbers)}", basis, positions, numbers)


def triclinic(n: int = 48, seed: int = 11) -> Structure:
    """The P1 path from tests/test_benchmark.cpp: every candidate setting is tried."""
    rng = np.random.default_rng(seed)
    basis = np.array([[5.1, 0.0, 0.0], [1.3, 6.4, 0.0], [0.7, 1.1, 7.9]])
    numbers = np.resize([8, 14, 26], n)
    return Structure(f"triclinic-{n}", basis, rng.uniform(size=(n, 3)), numbers)


CASES: dict[str, Structure] = {
    s.name: s for s in [*(rocksalt_supercell(k) for k in (2, 4, 6, 8)), triclinic()]
}


def cif_corpus() -> list[Path]:
    """PyXtal's 91 real-structure CIFs, the same set tests/test_cif_corpus.cpp reads."""
    root = Path(str(resources.files("pyxtal") / "database" / "cifs"))
    return sorted(root.glob("*.cif"))
