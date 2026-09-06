#!/usr/bin/env python3
"""Fill the metallic and van der Waals radius columns of tools/element_data.csv
from pymatgen's periodic table; the covalent column (Cordero et al. 2008) stays.
An element pymatgen has no value for is left empty (std::nullopt downstream).

    pip install pymatgen
    extract_element_radii.py [element_data.csv]
"""
import sys
from pathlib import Path

import pandas as pd
from pymatgen.core import Element

DEFAULT_CSV = Path(__file__).with_name("element_data.csv")
HEADER = """\
# Curated element reference data for Seitz crystal generation.
# covalent_radius: single-bond covalent radii (angstrom) from Cordero et
#   al., "Covalent radii revisited", Dalton Trans., 2008, 2832-2838.
# metallic_radius, vdw_radius: pymatgen's periodic table (angstrom),
#   extracted by tools/extract_element_radii.py; an empty field means
#   pymatgen tabulates no value and the accessor returns std::nullopt.
# Atomic volume is NOT tabulated: data/element_data.hpp derives it from
#   the covalent radius as the volume of a sphere.
"""


def main():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_CSV
    t = pd.read_csv(path, comment="#")
    elements = t.symbol.map(Element)
    assert (elements.map(lambda e: e.Z) == t.Z).all()
    t["metallic_radius"] = elements.map(lambda e: e.metallic_radius).astype(float).round(2)
    t["vdw_radius"] = elements.map(lambda e: e.van_der_waals_radius).astype(float).round(2)
    with open(path, "w", newline="") as f:
        f.write(HEADER)
        t.to_csv(f, index=False, float_format="%.2f")
    print(f"elements: {len(t)}; metallic radii: {t.metallic_radius.notna().sum()}; "
          f"vdW radii: {t.vdw_radius.notna().sum()}\nwritten: {path}")


if __name__ == "__main__":
    main()
