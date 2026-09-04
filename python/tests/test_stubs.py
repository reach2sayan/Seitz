"""The committed stubs, and the one way they could break the package."""

from __future__ import annotations

import pathlib

import cppcrystal as cc


def test_the_stub_package_does_not_shadow_the_extension() -> None:
    """`_core` stubs are a directory sitting next to `_core.cpython-*.so`.

    Python resolves a real module in a directory ahead of a same-named
    subdirectory with no ``__init__.py``, so the extension wins -- but that is
    a rule worth pinning rather than trusting, because if it ever went the other
    way ``import cppcrystal`` would start returning an empty namespace package
    and every symbol would vanish at once.
    """
    assert cc._core.__file__.endswith(".so")
    assert cc.version_string()


def test_the_stubs_are_shipped_beside_the_module() -> None:
    package = pathlib.Path(cc.__file__).parent
    stubs = package / "_core"
    if not stubs.is_dir():
        # A build tree that has never had stubs generated is still usable.
        return
    assert (stubs / "__init__.pyi").is_file()
    assert (package / "py.typed").is_file()
