"""Shared by python/tests and python/benchmarks: the build-tree guard.

At the repo root, not under python/: pytest puts a conftest's directory on
sys.path, and python/ holds the *source* seitz package, whose _core is stubs only.
"""

from __future__ import annotations

import os
import pathlib

import pytest

import seitz as sz


def pytest_configure(config: pytest.Config) -> None:
    """Fail loudly if we imported a different _core than the one under test.

    CTest points PYTHONPATH at the build tree it just compiled -- which for the
    sanitizer preset is the whole point of test_lifetimes.py.  But PYTHONPATH is
    only consulted by the path finder, and anything that installs a meta-path
    finder for this package (a scikit-build-core editable install is the one
    that bit us) overrides it silently: the suite passes, having exercised some
    other build entirely.  Compare the two and say so instead.
    """
    expected = os.environ.get("PYTHONPATH", "").split(os.pathsep)[0]
    if not expected:
        return
    loaded = pathlib.Path(sz._core.__file__).resolve()
    if not loaded.is_relative_to(pathlib.Path(expected).resolve()):
        raise pytest.UsageError(
            f"seitz._core was imported from {loaded}, but PYTHONPATH points "
            f"at {expected}. Something is shadowing the build tree -- most "
            f"likely an editable install of seitz in this interpreter. "
            f"Remove it (`uv pip uninstall seitz`); this project is reached "
            f"through the CMake build tree, not an editable install."
        )
