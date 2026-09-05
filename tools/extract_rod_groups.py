#!/usr/bin/env python3
"""Extract PyXtal's rod-group (1D-periodic) symmetry operations into
tools/rod_group_data.csv.

The 75 rod groups are not in spglib (unlike layer groups, which Seitz reaches
through spglib's negative-Hall datasets), so their source of truth is PyXtal,
read through its stable public API rather than its on-disk database format:

    from pyxtal.symmetry import Group
    g = Group(n, dim=1)            # rod group n (1..75)

Only the OPERATIONS (and the HM symbol) are extracted: the C++ layer derives the
rod Wyckoff positions from them in-house, as group::PointGroup does in 0D, so
PyXtal's less-documented Wyckoff attributes never enter the build.

This script is offline and run by hand -- it needs PyXtal, which a C++ build has
no business requiring. Its CSV is the checked-in input; the build runs
tools/transcribe_rod_groups.py over it.

    pip install pyxtal            # provides pymatgen (SymmOp.affine_matrix)
    extract_rod_groups.py [out.csv]
"""
import csv
import os
import sys

NUM_ROD_GROUPS = 75
TRANSLATION_DENOMINATOR = 12

DEFAULT_OUT = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "rod_group_data.csv"))


def require_pyxtal():
    try:
        from pyxtal.symmetry import Group  # noqa: F401
    except ImportError:
        sys.exit("error: PyXtal is not installed. Install it (and its pymatgen "
                 "dependency) with:\n    pip install pyxtal\n"
                 "then re-run this extractor.")
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
    """The operations of a rod group = those of its general position (maximum
    multiplicity, trivial site symmetry), as (rotation[9], translation[3]).

    Only pymatgen's SymmOp.affine_matrix (4x4) is read, which is stable across
    PyXtal/pymatgen versions.
    """
    general = max(group.Wyckoff_positions, key=lambda w: w.multiplicity)

    ops = []
    seen = set()
    for symm_op in general.ops:
        affine = symm_op.affine_matrix
        rotation = integerize_rotation(
            [[affine[r][c] for c in range(3)] for r in range(3)])
        translation = [rationalize_translation(affine[a][3]) for a in range(3)]
        key = (tuple(rotation), tuple(translation))
        if key not in seen:                 # distinct factor-group reps only
            seen.add(key)
            ops.append((rotation, translation))
    return ops


def main():
    out_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT
    Group = require_pyxtal()

    rows = []
    counts = []
    for number in range(1, NUM_ROD_GROUPS + 1):
        g = Group(number, dim=1)
        ops = operations_of(g)
        counts.append(len(ops))
        for rotation, translation in ops:
            rows.append([number, str(g.symbol), *rotation, *translation])

    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["number", "symbol", *["r%d" % i for i in range(9)],
                    *["t%d" % i for i in range(3)]])
        w.writerows(rows)

    print("rod groups: %d, operations: %d (min %d, max %d per group)"
          % (NUM_ROD_GROUPS, len(rows), min(counts), max(counts)))
    print("written: %s" % out_path)


if __name__ == "__main__":
    main()
