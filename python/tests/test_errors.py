"""The Result<T> -> exception bridge.

``test_an_empty_cell_raises_the_specific_class`` is load-bearing infrastructure,
not an API test.  Boost.LEAF keeps an error's payload in a thread-local slot,
one exported symbol per tag type, and the library's ``new_error`` and the
extension's ``try_handle_all`` have to resolve those to the SAME object.  They
do today -- ``libseitz.so`` exports them ``STB_GNU_UNIQUE``, one instance
process-wide -- but a stray ``-Bsymbolic`` or ``--exclude-libs`` on the
extension target would bind them locally, and then every typed error would
silently arrive as the plain base class with its payload gone.

That test is the one that catches it, because ``e_empty_cell`` is raised inside
the library (``src/core/validation.hpp``, compiled into ``libseitz.so``) and
handled in the extension, so it genuinely crosses the boundary.  By contrast the
``InvalidLatticeError`` assertions below never touch LEAF -- the binding layer
rejects a singular basis itself, before the pipeline sees it -- so they are API
tests, not infrastructure ones.
"""

from __future__ import annotations

import numpy as np
import pytest

import seitz as sz
from seitz import errors


def test_a_singular_basis_raises_and_carries_its_determinant() -> None:
    """The payload reaches Python as an attribute, not buried in a sentence."""
    singular = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    with pytest.raises(errors.InvalidLatticeError) as caught:
        sz.Lattice(singular)

    assert caught.value.determinant == pytest.approx(0.0)


def test_from_basis_answers_none_where_the_constructor_raises() -> None:
    singular = np.zeros((3, 3))
    assert sz.Lattice.from_basis(singular) is None
    with pytest.raises(errors.InvalidLatticeError):
        sz.Lattice(singular)


def test_an_empty_cell_raises_the_specific_class() -> None:
    """The cross-library tag test -- see this module's docstring.

    ``e_empty_cell`` is created inside libseitz.so and matched by a handler
    compiled into the extension.  Asserting the *specific* class is the point:
    if the two ever stopped sharing LEAF's thread-local slots, this would arrive
    as a bare SeitzError and only a specific assertion would notice.
    """
    empty = sz.Cell(sz.Lattice(np.eye(3)), np.zeros((0, 3)), [])
    with pytest.raises(errors.EmptyCellError):
        sz.analyze(empty).hall


def test_an_unknown_wyckoff_letter_raises_with_a_message() -> None:
    """An e_message with no tag beside it keeps its sentence."""
    group = sz.SpaceGroup.of(sz.HallNumber(sz.GroupFamily.space, 1))
    with pytest.raises(errors.SeitzError, match="Wyckoff"):
        group.wyckoff("z")


def test_the_hierarchy_is_rooted_at_seitzerror() -> None:
    for name in errors.__all__:
        cls = getattr(errors, name)
        assert issubclass(cls, errors.SeitzError)
    assert issubclass(errors.SeitzError, Exception)


def test_an_out_of_range_group_number_raises() -> None:
    with pytest.raises(errors.SeitzError):
        sz.SpaceGroup.from_number(sz.GroupFamily.space, 231)
