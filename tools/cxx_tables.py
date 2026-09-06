"""Reading C initializers and emitting constexpr C++ tables, for tools/transcribe_*.py."""
import functools
import numbers
import re

import numpy as np

def c_initializer(text, name, symbols=None):
    """The value of `... name[...] = {...};` in C source, as nested lists.
    Brace initializers are Python list literals once braces become brackets;
    identifiers (enumerators) resolve through `symbols`."""
    text = re.sub(r"/\*.*?\*/|//[^\n]*", "", text, flags=re.S)
    body = re.search(rf"\b{re.escape(name)}\s*(?:\[[^\]]*\])*\s*=\s*(\{{.*?\}})\s*;", text, re.S)
    return eval(body.group(1).replace("{", "[").replace("}", "]"), {"__builtins__": {}}, symbols or {})

class code(str):
    """A string emitted verbatim (an enumerator, an expression), not as a literal."""

@functools.singledispatch
def literal(v):
    raise TypeError(f"no C++ literal for {v!r}")

@literal.register(type(None))
def _(v): return "std::nullopt"

@literal.register(code)
def _(v): return str(v)

@literal.register(str)
def _(v): return '"' + v.replace("\\", "\\\\").replace('"', '\\"') + '"'

@literal.register(float)
@literal.register(np.floating)
def _(v): return repr(float(v))

@literal.register(int)
@literal.register(np.integer)
def _(v): return str(int(v))

def scalar(v):
    return v is None or isinstance(v, (str, numbers.Number, np.generic))


def aggregate(x):
    """Nested std::array initializer: scalars as literals, sequences double-braced."""
    return literal(x) if scalar(x) else "{{" + ", ".join(map(aggregate, x)) + "}}"

def struct(fields):
    """Struct initializer with brace elision: nested sequences single-braced."""
    return "{" + ", ".join(literal(f) if scalar(f) else struct(f) for f in fields) + "}"

def lines(items, per_line, indent="    "):
    items = list(items)
    return "".join(indent + ", ".join(items[i:i + per_line]) + ",\n" for i in range(0, len(items), per_line))

def array(name, ctype, items, per_line=1, comment="", size=None, indent="    "):
    """`inline constexpr std::array<ctype, N> name = {{...}};`, items already rendered."""
    items = list(items)
    head = f"inline constexpr std::array<{ctype}, {size or len(items)}> {name} = {{{{\n"
    return (f"// {comment}\n" if comment else "") + head + lines(items, per_line, indent) + "}};\n\n"

def preamble(banner, includes=("array",), public=False):
    out = "#pragma once\n\n" + "".join(f"// {b}\n" for b in banner) + "\n"
    out += "".join(f"#include <{i}>\n" for i in includes) + "\n"
    return out + ("#pragma GCC visibility push(default)\n\n" if public else "") + "namespace seitz::data {\n\n"

def epilogue(public=False):
    return "} // namespace seitz::data\n" + ("\n#pragma GCC visibility pop\n" if public else "")
