"""The other half of the error model: failure as a value, not as a jump.

The C++ library never throws.  A fallible call returns ``Result<T>``
(``boost::leaf::result<T>``) and the caller decides, at the call site, what a
failure means.  ``python/src/errors.hpp`` performs the one translation that
turns such a result into a raised :class:`~pyseitz.errors.SeitzError` subclass,
because raising is what Python callers expect by default.

This module is that translation's named inverse.  Every helper here returns a
``Result[T, SeitzError]`` -- ``Ok`` or ``Err`` from the `result
<https://pypi.org/project/result/>`_ package -- so the C++ discipline is
available to Python code that wants it::

    from result import Ok, Err
    from pyseitz import errors, results

    match results.read_cif(path):
        case Ok(structures):
            print(len(structures))
        case Err(errors.CifSyntaxError() as error):
            print(error.line, error.column)

Nothing about the raising surface changes; ``pyseitz.read_cif`` still raises.  The
two share one hierarchy, and share it by identity rather than by parallel
definition: an ``Err`` **contains** the very exception object that would have
been raised, payload attribute and all, so ``raise res.err_value`` is the round
trip back.

Only :class:`~pyseitz.errors.SeitzError` becomes an ``Err``.  A ``TypeError``, an
``IndexError``, pydantic's ``ValidationError`` from a malformed tolerance, the
``OSError`` from reading a CIF off disk -- all still raise.  A bug in the caller
is not a crystallographic failure, and a Result API that flattens the two stops
carrying information.
"""

from __future__ import annotations

from typing import Callable, ParamSpec, TypeAlias, TypeVar

from result import Err, Ok, Result, as_result

from . import _core
from . import read_cif as _raising_read_cif
from . import write_cif as _raising_write_cif
from .errors import SeitzError

_P = ParamSpec("_P")
_T = TypeVar("_T")

#: ``Result`` with the error half already spelled: ``SeitzResult[Cell]``.
SeitzResult: TypeAlias = Result[_T, SeitzError]

__all__ = ["Err", "Ok", "Result", "SeitzResult", "attempt", "checked",
           "parse_cif", "read_cif", "write_cif"]


def attempt(call: Callable[[], _T]) -> Result[_T, SeitzError]:
    """Run ``call``, keeping a :class:`~pyseitz.errors.SeitzError` as a value.

    Takes the *call*, not its result -- the same shape, for the same reason, as
    the C++ ``unwrap()`` this inverts.  That is what lets it cover what a
    decorator cannot: an analyzer's memoized projections are properties, so
    there is no function to wrap, only an access to defer::

        results.attempt(lambda: analyzer.hall)
        results.attempt(lambda: lattice.niggli())
        results.attempt(lambda: pyseitz.Lattice(basis))
    """
    try:
        return Ok(call())
    except SeitzError as error:
        return Err(error)


def checked(fn: Callable[_P, _T]) -> Callable[_P, Result[_T, SeitzError]]:
    """``fn``, rewritten to return a ``Result`` instead of raising.

    The decorator form, for a callable used more than once.  The signature
    survives, statically and at runtime, so the wrapper is a drop-in whose only
    difference is the return type::

        wyckoff = results.checked(group.wyckoff)
        for letter in "abc":
            match wyckoff(letter): ...
    """
    return as_result(SeitzError)(fn)


# The fallible entry points of `pyseitz`, pre-wrapped. These three are the whole
# set at module level: everything else that can fail is a method, a static
# factory or a memoized property, and reaches this module through attempt().
parse_cif = checked(_core.parse_cif)
read_cif = checked(_raising_read_cif)
write_cif = checked(_raising_write_cif)
