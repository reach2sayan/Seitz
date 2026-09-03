#!/usr/bin/env python3
"""Transcribe PyXtal's rod-group (1D-periodic) symmetry operations into a
constexpr header.

Rod groups are the 75 subperiodic groups periodic along ONE axis. They are NOT
in spglib (unlike layer groups, which CppCrystal reaches through spglib's
negative-Hall datasets), so the source of truth here is PyXtal — the same
crystal-generation library this layer mirrors — read through its stable public
API rather than its on-disk database format:

    from pyxtal.symmetry import Group
    g = Group(n, dim=1)            # rod group n (1..75)

For each rod group we emit only its symmetry OPERATIONS (and HM symbol). That is
deliberately the minimal sufficient data: the C++ layer derives the rod Wyckoff
positions IN-HOUSE from these operations (the same approach group::PointGroup
uses for the 0D point groups — stabiliser/orbit of a generic point), so the
brittle, less-documented PyXtal Wyckoff attribute names never enter the build.
Wyckoff letters are then assigned by ascending multiplicity in C++.

Operation convention (matches the spglib operation tables): each entry holds one
period's worth of factor-group representatives; the pure (0,0,1) lattice
translation along the periodic axis is implicit. PyXtal places the periodic axis
along c (index 2); a, b are the two aperiodic directions. The orbit expansion in
the generator must therefore fold ONLY the c component (the rod analogue of the
layer-group aperiodic-axis fix: a rod-flipping op sends y -> -y in Cartesian, and
that must not be wrapped to 1 - y).

Usage:
    pip install pyxtal            # provides pymatgen (SymmOp.affine_matrix) too
    transcribe_rod_groups.py [out.hpp]

`out.hpp` defaults to src/data/rod_group_tables.hpp relative to the repo root.
"""
import os
import sys

NUM_ROD_GROUPS = 75
TRANSLATION_DENOMINATOR = 12
PERIODIC_AXIS = 2  # c; PyXtal rod-group convention

# Default output, resolved relative to this file (tools/ -> repo root).
DEFAULT_OUT = os.path.normpath(
    os.path.join(os.path.dirname(__file__), os.pardir,
                 "src", "data", "rod_group_tables.hpp"))


def require_pyxtal():
    try:
        from pyxtal.symmetry import Group  # noqa: F401
    except ImportError:
        sys.exit("error: PyXtal is not installed. Install it (and its pymatgen "
                 "dependency) with:\n    pip install pyxtal\n"
                 "then re-run this transcriber.")
    return Group


def rationalize_translation(value):
    """A translation component (a float fraction) as an integer numerator over
    TRANSLATION_DENOMINATOR, folded into [0, denominator)."""
    frac = value - int(value)            # strip whole periods
    if frac < 0:
        frac += 1.0
    num = round(frac * TRANSLATION_DENOMINATOR)
    num %= TRANSLATION_DENOMINATOR
    return num


def integerize_rotation(matrix):
    """A 3x3 rotation (row-major list of 9 ints), asserting it is integral."""
    out = []
    for r in range(3):
        for c in range(3):
            v = matrix[r][c]
            i = round(v)
            assert abs(v - i) < 1e-6, ("non-integer rotation entry", v)
            out.append(int(i))
    return out


def operations_of(group):
    """The full symmetry-operation set of a rod group = the operations of its
    general position (maximum multiplicity, trivial site symmetry). Returned as a
    list of (rotation[9], translation[3]) tuples.

    Uses only pymatgen's SymmOp.affine_matrix (a 4x4 numpy array), which is
    stable across PyXtal/pymatgen versions, rather than any PyXtal-specific
    Wyckoff attribute.
    """
    wyckoffs = group.Wyckoff_positions
    general = max(wyckoffs, key=lambda w: w.multiplicity)

    ops = []
    seen = set()
    for symm_op in general.ops:
        affine = symm_op.affine_matrix  # 4x4
        rotation = integerize_rotation(
            [[affine[r][c] for c in range(3)] for r in range(3)])
        translation = [rationalize_translation(affine[a][3]) for a in range(3)]
        key = (tuple(rotation), tuple(translation))
        if key not in seen:                 # distinct factor-group reps only
            seen.add(key)
            ops.append((rotation, translation))
    return ops


def collect():
    """Return (symbols, ops_per_group) for rod groups 1..75."""
    Group = require_pyxtal()
    symbols = []
    ops_per_group = []
    for number in range(1, NUM_ROD_GROUPS + 1):
        g = Group(number, dim=1)
        symbols.append(str(g.symbol))
        ops_per_group.append(operations_of(g))
    return symbols, ops_per_group


def fmt_rotation(rotation):
    return "{" + ", ".join("%d" % v for v in rotation) + "}"


def fmt_translation(translation):
    return "{" + ", ".join("%d" % v for v in translation) + "}"


def emit(symbols, ops_per_group, out_path):
    flat = [op for group_ops in ops_per_group for op in group_ops]
    offsets = [0]
    for group_ops in ops_per_group:
        offsets.append(offsets[-1] + len(group_ops))

    with open(out_path, "w") as f:
        w = f.write
        w("#pragma once\n\n")
        w("// GENERATED by tools/transcribe_rod_groups.py from PyXtal "
          "(Group(n, dim=1)).\n// Do not edit by hand.\n\n")
        w("#include <array>\n")
        w("#include <string_view>\n\n")
        w("namespace cppcrystal::data {\n\n")

        w("// Number of rod (1D-periodic subperiodic) groups.\n")
        w("inline constexpr int kNumRodGroups = %d;\n\n" % NUM_ROD_GROUPS)

        w("// Rod groups are 1D-periodic along the c axis (index 2); a, b "
          "(indices 0, 1)\n// are the two aperiodic directions.\n")
        w("inline constexpr int kRodPeriodicAxis = %d;\n\n" % PERIODIC_AXIS)

        w("// Translation components are integer numerators over this "
          "denominator.\n")
        w("inline constexpr int kRodTranslationDenominator = %d;\n\n"
          % TRANSLATION_DENOMINATOR)

        w("// HM symbol indexed by (rod_number - 1).\n")
        w("inline constexpr std::array<std::string_view, kNumRodGroups> "
          "kRodGroupSymbols = {{\n")
        for i in range(0, len(symbols), 5):
            w("    " + ", ".join('"%s"' % s for s in symbols[i:i + 5]) + ",\n")
        w("}};\n\n")

        w("// One symmetry operation: 9 row-major rotation entries + 3 "
          "translation\n// numerators over kRodTranslationDenominator. These are "
          "one period's worth\n// of factor-group representatives (the pure "
          "(0,0,1) lattice translation is\n// implicit), as in the space-group "
          "operation tables.\n")
        w("struct RodOperation {\n")
        w("  std::array<int, 9> rotation;\n")
        w("  std::array<int, 3> translation;\n")
        w("};\n\n")

        w("inline constexpr std::array<RodOperation, %d> kRodOperations = "
          "{{\n" % len(flat))
        for rotation, translation in flat:
            w("    {%s, %s},\n"
              % (fmt_rotation(rotation), fmt_translation(translation)))
        w("}};\n\n")

        w("// Half-open range [kRodOperationOffset[g-1], kRodOperationOffset[g])"
          "\n// into kRodOperations for rod group g (1..%d).\n"
          % NUM_ROD_GROUPS)
        w("inline constexpr std::array<int, kNumRodGroups + 1> "
          "kRodOperationOffset = {{\n")
        for i in range(0, len(offsets), 12):
            w("    " + ", ".join("%d" % o for o in offsets[i:i + 12]) + ",\n")
        w("}};\n\n")

        w("} // namespace cppcrystal::data\n")

    return len(flat)


def main():
    out_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT
    symbols, ops_per_group = collect()
    total = emit(symbols, ops_per_group, out_path)
    counts = [len(o) for o in ops_per_group]
    print("rod groups: %d, operations: %d (min %d, max %d per group)"
          % (NUM_ROD_GROUPS, total, min(counts), max(counts)))
    print("written: %s" % out_path)


if __name__ == "__main__":
    main()
