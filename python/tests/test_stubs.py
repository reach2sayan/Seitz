"""The committed stubs, and the one way they could break the package."""

from __future__ import annotations

import importlib.machinery
import pathlib

import seitz as sz


def test_the_stub_package_does_not_shadow_the_extension() -> None:
    """`_core` stubs are a directory sitting next to the built `_core` module.

    Python resolves a real module in a directory ahead of a same-named
    subdirectory with no ``__init__.py``, so the extension wins -- but that is
    a rule worth pinning rather than trusting, because if it ever went the other
    way ``import seitz`` would start returning an empty namespace package
    and every symbol would vanish at once.

    The suffix is asked of the interpreter rather than spelled out: it is
    ``.so`` on Linux and macOS but ``.pyd`` on Windows, and what the assertion
    means is "a compiled extension won", not "this one platform's suffix".
    """
    assert sz._core.__file__.endswith(
        tuple(importlib.machinery.EXTENSION_SUFFIXES)
    )
    assert sz.version_string()


def test_the_stubs_are_shipped_beside_the_module() -> None:
    package = pathlib.Path(sz.__file__).parent
    stubs = package / "_core"
    if not stubs.is_dir():
        # A build tree that has never had stubs generated is still usable.
        return
    assert (stubs / "__init__.pyi").is_file()
    assert (package / "py.typed").is_file()
