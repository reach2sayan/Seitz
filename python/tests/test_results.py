"""The exception -> Result inverse.

``test_errors.py`` pins that a failure arrives as the right exception class with
its payload attached.  This module pins that ``pyseitz.results`` hands that same
object back as a value instead, and -- the half that is easy to get wrong -- that
it converts SeitzError and nothing else.
"""

from __future__ import annotations

import numpy as np
import pydantic
import pytest
from result import Err, Ok

import pyseitz as sz
from pyseitz import errors, results

# The same malformed documents test_errors.py uses, so the two suites cannot
# drift about what "a CIF syntax error" is.
BAD_SYNTAX = "data_x\n_a 1\nloop_ _p _q\n1 2 3\n"
MISSING_TAG = (
    "data_x\n_cell_length_a 4.0\n_cell_length_c 4.0\n"
    "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
    "loop_ _atom_site_label _atom_site_fract_x _atom_site_fract_y "
    "_atom_site_fract_z\nNa1 0 0 0\n"
)


def test_success_is_ok_and_holds_what_the_raising_call_returns(bcc_fe: sz.Cell) -> None:
    text = sz.write_cif(bcc_fe, name="fe")
    wrapped = results.write_cif(bcc_fe, name="fe")

    assert isinstance(wrapped, Ok)
    assert wrapped.unwrap() == text


def test_failure_is_err_carrying_the_exception_itself() -> None:
    """Not a description of the error -- the very object `raise` would throw."""
    wrapped = results.parse_cif(BAD_SYNTAX)

    assert isinstance(wrapped, Err)
    assert isinstance(wrapped.err_value, errors.CifSyntaxError)
    assert wrapped.err_value.line == 3
    assert wrapped.err_value.column == 1


def test_the_payload_survives_on_the_other_entry_points() -> None:
    wrapped = results.read_cif(MISSING_TAG)

    assert isinstance(wrapped, Err)
    assert wrapped.err_value.tag == "_cell_length_b"


def test_match_binds_the_class_and_its_payload() -> None:
    """Structural matching is half the reason for a Result type at all."""
    match results.parse_cif(BAD_SYNTAX):
        case Ok(_):
            pytest.fail("that document does not parse")
        case Err(errors.CifSyntaxError() as error):
            assert (error.line, error.column) == (3, 1)
        case Err(other):
            pytest.fail(f"wrong class: {type(other).__name__}")


def test_attempt_covers_the_memoized_properties() -> None:
    """The analyzer's projections are properties, so only attempt() reaches them.

    This is also the one assertion here that crosses the libseitz.so ->
    extension LEAF boundary that test_errors.py's docstring describes: a
    specific class, not a bare SeitzError.
    """
    empty = sz.Cell(sz.Lattice(np.eye(3)), np.zeros((0, 3)), [])
    wrapped = results.attempt(lambda: sz.analyze(empty).hall)

    assert isinstance(wrapped, Err)
    assert isinstance(wrapped.err_value, errors.EmptyCellError)


def test_attempt_covers_constructors_and_carries_their_payload() -> None:
    singular = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    wrapped = results.attempt(lambda: sz.Lattice(singular))

    assert isinstance(wrapped, Err)
    assert wrapped.err_value.determinant == pytest.approx(0.0)


def test_attempt_is_transparent_on_success(bcc_fe: sz.Cell) -> None:
    analyzer = sz.analyze(bcc_fe)
    assert results.attempt(lambda: analyzer.hall).unwrap() == analyzer.hall


def test_checked_wraps_a_bound_method_and_keeps_its_signature() -> None:
    group = sz.SpaceGroup.of(sz.HallNumber(sz.GroupFamily.space, 1))
    wyckoff = results.checked(group.wyckoff)

    assert isinstance(wyckoff("a"), Ok)
    assert isinstance(wyckoff("z"), Err)
    assert isinstance(wyckoff("z").err_value, errors.SeitzError)


def test_err_round_trips_back_to_a_raise() -> None:
    """`raise res.err_value` is the bridge back, because the value is the error."""
    wrapped = results.parse_cif(BAD_SYNTAX)
    assert isinstance(wrapped, Err)

    with pytest.raises(errors.CifSyntaxError) as caught:
        raise wrapped.err_value
    assert caught.value.line == 3


# Only SeitzError becomes an Err. A bug in the caller stays an exception.


def test_a_type_error_still_raises_through_checked() -> None:
    with pytest.raises(TypeError):
        results.write_cif(object(), name="x")  # type: ignore[arg-type]


def test_a_validation_error_still_raises_through_checked() -> None:
    with pytest.raises(pydantic.ValidationError):
        results.read_cif("data_x\n", {"symprec": -1.0})


def test_a_non_seitz_exception_still_raises_through_attempt() -> None:
    def boom() -> int:
        raise ZeroDivisionError("not a crystallographic failure")

    with pytest.raises(ZeroDivisionError):
        results.attempt(boom)


def test_the_module_exports_the_result_vocabulary() -> None:
    """A caller should not need to import `result` itself to spell a match arm."""
    for name in results.__all__:
        assert hasattr(results, name)
    assert results.Ok is Ok
    assert results.Err is Err
