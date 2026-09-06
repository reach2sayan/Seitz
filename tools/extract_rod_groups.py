#!/usr/bin/env python3
"""PyXtal's rod-group operations -> tools/rod_group_data.csv.

The 75 rod groups are not in spglib; PyXtal's `Group(n, dim=1)` is the source.
Only the operations of the general position (and the HM symbol) are taken: the
C++ side derives the Wyckoff positions itself. Offline, run by hand; the build
transcribes the CSV with tools/transcribe_rod_groups.py.

    pip install pyxtal
    extract_rod_groups.py [out.csv]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pyxtal.symmetry import Group

NUM_ROD_GROUPS = 75
TRANSLATION_DENOMINATOR = 12
DEFAULT_OUT = Path(__file__).with_name("rod_group_data.csv")

def operations(number):
    """(number, symbol, r0..r8, t0..t2) rows: distinct factor-group representatives."""
    g = Group(number, dim=1)
    affine = np.array([op.affine_matrix for op in max(g.Wyckoff_positions, key=lambda w: w.multiplicity).ops])
    rotation, translation = affine[:, :3, :3], affine[:, :3, 3]
    assert np.allclose(rotation, np.rint(rotation), atol=1e-6)
    numerators = np.rint((translation % 1) * TRANSLATION_DENOMINATOR).astype(int) % TRANSLATION_DENOMINATOR
    rows = np.c_[np.rint(rotation).reshape(-1, 9).astype(int), numerators]
    return pd.DataFrame(rows, columns=[f"r{i}" for i in range(9)] + [f"t{i}" for i in range(3)]) \
        .drop_duplicates().assign(number=number, symbol=str(g.symbol))

def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_OUT
    t = pd.concat(map(operations, range(1, NUM_ROD_GROUPS + 1)), ignore_index=True)
    t = t[["number", "symbol", *t.columns[:12]]]
    t.to_csv(out, index=False)
    counts = t.groupby("number").size()
    print(f"rod groups: {NUM_ROD_GROUPS}, operations: {len(t)} (min {counts.min()}, max {counts.max()} per group)\nwritten: {out}")

if __name__ == "__main__":
    main()
