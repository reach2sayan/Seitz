"""pytest-benchmark cases: engine call only, native objects prebuilt."""

from __future__ import annotations

import pytest

from backends import BACKENDS, MESH, Backend
from structures import CASES, cif_corpus

WORKLOADS = ("dataset", "standardize", "primitive")
CIF_BACKENDS = [n for n, b in BACKENDS.items() if b.read_cif]
MESH_BACKENDS = [n for n, b in BACKENDS.items() if b.reciprocal_mesh]
MAGNETIC_BACKENDS = [n for n, b in BACKENDS.items() if b.magnetic_dataset]
# One small and one large cell each: the mesh reduction scales with the grid
# and the point group, the magnetic search with the atom count.
MESH_CASES = ("rocksalt-16", "rocksalt-128")
MAGNETIC_CASES = ("rocksalt-16", "rocksalt-128")


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
    _say(f"cif corpus: {len(files)} of {len(cif_corpus())} files readable by every backend ({', '.join(CIF_BACKENDS)})")
    return files


def _say(text: str, end: str = "\n") -> None:
    print(text, end=end, flush=True)


def _done(benchmark, extra: str = "") -> None:
    st = benchmark.stats.stats
    _say(f"median {st.median * 1e3:10.3f} ms   min {st.min * 1e3:10.3f} ms   {st.rounds:5d} rounds{extra}")


@pytest.mark.parametrize("backend", BACKENDS)
@pytest.mark.parametrize("case", CASES)
@pytest.mark.parametrize("workload", WORKLOADS)
def test_symmetry(benchmark, workload, case, backend, reference_numbers) -> None:
    b: Backend = BACKENDS[backend]
    fn = getattr(b, workload)
    label = f"{workload:11s} {case:14s} {backend:9s}"
    if fn is None:
        _say(f"{label} skipped: no {workload}"); pytest.skip(f"{backend} has no {workload}")
    if len(CASES[case]) > b.max_atoms:
        _say(f"{label} skipped: capped at {b.max_atoms} atoms"); pytest.skip(f"{backend} capped at {b.max_atoms} atoms")
    _say(f"{label} {len(CASES[case]):5d} atoms ... ", end="")
    benchmark.group = f"{workload}:{case}"
    benchmark.extra_info["backend"] = backend
    result = benchmark(fn, b.prepare(CASES[case]))
    found = b.number(result) if workload == "dataset" else None
    _done(benchmark, f"   space group {found}" if found else "")
    if found is not None and case in reference_numbers:
        assert found == reference_numbers[case], f"spglib says {reference_numbers[case]}"


@pytest.mark.parametrize("backend", CIF_BACKENDS)
def test_read_cif(benchmark, backend, cif_files) -> None:
    read = BACKENDS[backend].read_cif
    _say(f"{'read_cif':11s} {len(cif_files):3d} files       {backend:9s} ... ", end="")
    benchmark.group = "read_cif"
    benchmark.extra_info["backend"] = backend
    benchmark.extra_info["files"] = len(cif_files)
    benchmark(lambda: [read(p) for p in cif_files])
    _done(benchmark)


@pytest.mark.parametrize("backend", MESH_BACKENDS)
@pytest.mark.parametrize("case", MESH_CASES)
def test_reciprocal_mesh(benchmark, case, backend) -> None:
    """Irreducible-wedge reduction of a {0}^3 grid, the determination prebuilt."""
    b: Backend = BACKENDS[backend]
    grid = "x".join(str(d) for d in MESH)
    _say(f"{'mesh':11s} {case:14s} {backend:9s} {grid:>11s} ... ", end="")
    benchmark.group = f"mesh:{case}"
    benchmark.extra_info["backend"] = backend
    benchmark.extra_info["mesh"] = grid
    result = benchmark(b.reciprocal_mesh, b.prepare_mesh(CASES[case]))
    _done(benchmark, f"   {b.num_irreducible(result)} irreducible")


@pytest.mark.parametrize("backend", MAGNETIC_BACKENDS)
@pytest.mark.parametrize("case", MAGNETIC_CASES)
def test_magnetic_dataset(benchmark, case, backend) -> None:
    """Magnetic determination of a collinear antiferromagnet."""
    b: Backend = BACKENDS[backend]
    _say(f"{'magnetic':11s} {case:14s} {backend:9s} {len(CASES[case]):5d} atoms ... ", end="")
    benchmark.group = f"magnetic:{case}"
    benchmark.extra_info["backend"] = backend
    result = benchmark(b.magnetic_dataset, b.prepare_magnetic(CASES[case]))
    _done(benchmark, f"   UNI {b.uni_number(result)}")
