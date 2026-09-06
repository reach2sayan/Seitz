"""pytest-benchmark cases: engine call only, native objects prebuilt."""

from __future__ import annotations

import pytest

from backends import BACKENDS, Backend
from structures import CASES, cif_corpus

WORKLOADS = ("dataset", "standardize", "primitive")
CIF_BACKENDS = [n for n, b in BACKENDS.items() if b.read_cif]


@pytest.fixture(scope="session")
def reference_numbers() -> dict[str, int]:
    """spglib's answer per case, so a backend that disagrees fails rather than times."""
    ref = BACKENDS.get("spglib")
    if ref is None: return {}
    return {name: ref.number(ref.dataset(ref.prepare(s))) for name, s in CASES.items()}


@pytest.fixture(scope="session")
def cif_files() -> list:
    """Only the files every CIF backend parses, so all of them time the same set."""
    def readable(p):
        try: return all(BACKENDS[n].read_cif(p) is not None for n in CIF_BACKENDS)
        except Exception: return False
    files = [p for p in cif_corpus() if readable(p)]
    print(f"\ncif corpus: {len(files)} of {len(cif_corpus())} files readable by all backends")
    return files


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("case", CASES)
@pytest.mark.parametrize("workload", WORKLOADS)
def test_symmetry(benchmark, workload, case, backend, reference_numbers) -> None:
    b: Backend = BACKENDS[backend]
    fn = getattr(b, workload)
    if fn is None: pytest.skip(f"{backend} has no {workload}")
    if len(CASES[case]) > b.max_atoms: pytest.skip(f"{backend} capped at {b.max_atoms} atoms")
    benchmark.group = f"{workload}:{case}"
    benchmark.extra_info["backend"] = backend
    result = benchmark(fn, b.prepare(CASES[case]))
    if workload == "dataset" and case in reference_numbers:
        assert b.number(result) == reference_numbers[case]


@pytest.mark.parametrize("backend", CIF_BACKENDS)
def test_read_cif(benchmark, backend, cif_files) -> None:
    read = BACKENDS[backend].read_cif
    benchmark.group = "read_cif"
    benchmark.extra_info["backend"] = backend
    benchmark.extra_info["files"] = len(cif_files)
    benchmark(lambda: [read(p) for p in cif_files])
