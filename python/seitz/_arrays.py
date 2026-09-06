"""Pydantic schemas for the NumPy shapes the library speaks.

Pydantic v2 cannot build a core schema for ``NDArray`` on its own, and
``arbitrary_types_allowed`` would only weaken it to an isinstance check with no
coercion and no JSON.  So each alias below replaces the schema outright with
pydantic's own annotation types -- ``PlainValidator`` (which tolerates a source
type it cannot introspect), ``PlainSerializer``, ``WithJsonSchema`` -- and the
models stay strict and still round-trip through JSON.
"""

from __future__ import annotations

from typing import Annotated, Any

import numpy as np
from numpy.typing import NDArray
from pydantic import PlainSerializer, PlainValidator, WithJsonSchema

__all__ = ["Basis", "Positions", "Vector3"]

_ROW = {"type": "array", "items": {"type": "number"}, "minItems": 3, "maxItems": 3}
_MATRIX = _ROW | {"items": _ROW}


def _block(shape: tuple[int, ...], order: str, described: str) -> Any:
    """Validator coercing to a read-only float64 array of ``shape`` in ``order``.

    ``-1`` is the free dimension.  ``order`` is the memory layout the C++ side
    wants: ``"C"`` for Positions (``Matrix<double, Dynamic, 3, RowMajor>``) and
    ``"F"`` for Matrix3d, whose columns are the basis vectors -- so
    ``basis[:, 0]`` is **a**, with no transpose anywhere.  Either way pybind11
    can memcpy the block instead of converting it element by element.

    ``.view()`` before ``setflags``: ``asarray`` hands back the caller's own
    array when it already has the right dtype and layout, and freezing that
    would freeze a variable the caller still owns.
    """

    def validate(value: Any) -> NDArray[np.float64]:
        array = np.asarray(value, dtype=np.float64, order=order).view()
        if array.size == 0 and -1 in shape:
            array = array.reshape(shape)
        if array.ndim != len(shape) or any(n not in (-1, m) for n, m in zip(shape, array.shape)):
            raise ValueError(f"expected {described}, got shape {array.shape}")
        array.setflags(write=False)
        return array

    return validate


def _array(shape: tuple[int, ...], order: str, described: str, json_schema: dict[str, Any]) -> Any:
    return Annotated[NDArray[np.float64],
                     PlainValidator(_block(shape, order, described)),
                     PlainSerializer(np.ndarray.tolist, when_used="json"),
                     WithJsonSchema(json_schema)]


Positions = _array((-1, 3), "C", "an (N, 3) array", {"type": "array", "items": _ROW})
Basis = _array((3, 3), "F", "a (3, 3) matrix", _MATRIX)
Vector3 = _array((3,), "C", "a 3-vector", _ROW)
