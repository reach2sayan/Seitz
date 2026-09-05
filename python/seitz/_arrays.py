"""Pydantic schemas for the NumPy shapes the library speaks.

Pydantic v2 cannot build a core schema for ``NDArray`` on its own, and
``arbitrary_types_allowed`` would only weaken it to an isinstance check with no
coercion and no JSON.  The annotation classes here replace the schema outright
instead, so the models stay strict and still round-trip through JSON.
"""

from __future__ import annotations

from typing import Annotated, Any

import numpy as np
from numpy.typing import NDArray
from pydantic import GetCoreSchemaHandler, GetJsonSchemaHandler
from pydantic.json_schema import JsonSchemaValue
from pydantic_core import CoreSchema, core_schema

__all__ = ["Basis", "Positions", "Vector3"]

def _as_positions(value: Any) -> NDArray[np.float64]:
    """Coerce to the C-contiguous (N, 3) float64 block ``Positions`` is in C++.

    ``ascontiguousarray`` rather than ``asarray``: Positions is
    ``Matrix<double, Dynamic, 3, RowMajor>``, so a C-contiguous float64 input is
    the one shape pybind11 can memcpy instead of converting element by element.
    """
    array = np.ascontiguousarray(value, dtype=np.float64)
    if array.size == 0: array = array.reshape(0, 3)
    if array.ndim != 2 or array.shape[1] != 3:
        raise ValueError(f"expected an (N, 3) array, got shape {array.shape}")
    array.setflags(write=False)
    return array

def _as_basis(value: Any) -> NDArray[np.float64]:
    """Coerce to the (3, 3) F-ordered block ``Matrix3d`` is in C++.

    Fortran order because Eigen's Matrix3d is column-major and the columns are
    the basis vectors -- ``basis[:, 0]`` is **a**, with no transpose anywhere.
    """
    array = np.asfortranarray(value, dtype=np.float64)
    if array.shape != (3, 3):
        raise ValueError(f"expected a (3, 3) matrix, got shape {array.shape}")
    array.setflags(write=False)
    return array

def _as_vector3(value: Any) -> NDArray[np.float64]:
    array = np.ascontiguousarray(value, dtype=np.float64)
    if array.shape != (3,):
        raise ValueError(f"expected a 3-vector, got shape {array.shape}")
    array.setflags(write=False)
    return array


class _ArrayAnnotation:
    """Gives pydantic a schema for an array type it cannot introspect.

    ``no_info_plain_validator_function`` replaces the schema outright, so
    pydantic never tries to build one for ``NDArray`` and the model needs no
    ``arbitrary_types_allowed`` escape hatch.
    """

    def __init__(self, validator: Any, json_schema: JsonSchemaValue) -> None:
        self._validator = validator
        self._json_schema = json_schema

    def __get_pydantic_core_schema__(self, source: Any, handler: GetCoreSchemaHandler) -> CoreSchema:
        return core_schema.no_info_plain_validator_function(self._validator,
            serialization=core_schema.plain_serializer_function_ser_schema(
                lambda array: array.tolist(), when_used="json"),
        )

    def __get_pydantic_json_schema__(self, schema: CoreSchema, handler: GetJsonSchemaHandler) -> JsonSchemaValue:
        return dict(self._json_schema)

_ROW = {"type": "array", "items": {"type": "number"}, "minItems": 3, "maxItems": 3}

Positions = Annotated[NDArray[np.float64], _ArrayAnnotation(_as_positions, {"type": "array", "items": _ROW}),]
Basis = Annotated[NDArray[np.float64], _ArrayAnnotation(_as_basis,{"type": "array", "items": _ROW, "minItems": 3, "maxItems": 3},),]
Vector3 = Annotated[NDArray[np.float64], _ArrayAnnotation(_as_vector3, _ROW)]
