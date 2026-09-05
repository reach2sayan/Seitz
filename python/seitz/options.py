"""Validated inputs.

These are the small configuration values a caller assembles by hand, so they get
pydantic: keyword arguments, defaults that mirror the C++ ones, range checks,
construction from a plain dict, and a JSON round-trip for provenance.  Each
carries a ``to_core()`` that produces the plain struct the extension expects --
validation happens once, at the boundary, never inside a loop.
"""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field

from . import _core

__all__ = ["MagneticTolerance", "Tolerance"]

#: The library's own default, read from the extension so there is one home.
K_DEFAULT_SYMPREC: float = _core.K_DEFAULT_SYMPREC


class Tolerance(BaseModel):
    """Tolerances threaded through the symmetry search.

    ``angle_tolerance`` left unset means "derive an effective value from
    symprec", which is what ``std::nullopt`` means on the C++ side -- not zero,
    and not a sentinel.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    symprec: float = Field(default=K_DEFAULT_SYMPREC, gt=0.0)
    angle_tolerance: float | None = Field(default=None, gt=0.0)

    def to_core(self) -> _core.Tolerance:
        return _core.Tolerance(self.symprec, self.angle_tolerance)


class MagneticTolerance(Tolerance):
    """Tolerance plus the magnetic search's moment comparison.

    ``moment`` unset falls back to ``symprec``, as it does in C++.
    """

    moment: float | None = Field(default=None, gt=0.0)

    def to_core(self) -> _core.MagneticTolerance:
        return _core.MagneticTolerance(
            self.symprec, self.angle_tolerance, self.moment
        )
