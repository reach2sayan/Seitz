#!/usr/bin/env python3
"""Stage the native NuGet package tree from two install prefixes.

`nuget pack contrib/nuget/Seitz.nuspec -BasePath <stage>` expects:

    <stage>/
        README.md, LICENSE
        build/native/Seitz.targets
        build/native/include/**                 seitz + boost + eigen3
        build/native/lib/x64/{Release,Debug}/seitz.lib

Headers come from the Release prefix; the two prefixes differ only in the
binaries. Each `seitz.lib` is every archive in that prefix's lib/ merged by
lib.exe -- the library itself plus the Boost archives it was built against,
which install alongside it because Dependencies.cmake leaves Boost's install
rules on. Merging is what lets contrib/nuget/Seitz.targets name one library
instead of tracking Boost's installed set.

Missing pieces stop the script: a package staged short would pack and publish
without complaint.

Usage:
    python tools/pack_nuget.py --release <prefix> --debug <prefix> --stage <dir>
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def merge(prefix, out):
    """Every .lib in `prefix/lib` as one archive at `out`."""
    libs = sorted((prefix / "lib").glob("*.lib"))
    if not any(lib.name == "seitz.lib" for lib in libs):
        print(f"missing: {prefix / 'lib' / 'seitz.lib'}", file=sys.stderr)
        return False
    out.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["lib.exe", f"/OUT:{out}", "/NOLOGO", *map(str, libs)], check=True)
    print(f"{out}: {out.stat().st_size / 1e6:.1f} MB from {len(libs)} archives")
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", type=Path, required=True, help="Release install prefix")
    ap.add_argument("--debug", type=Path, required=True, help="Debug install prefix")
    ap.add_argument("--stage", type=Path, required=True, help="output directory (recreated)")
    args = ap.parse_args()

    shutil.rmtree(args.stage, ignore_errors=True)
    native = args.stage / "build" / "native"
    native.mkdir(parents=True)

    shutil.copy2(ROOT / "contrib" / "nuget" / "Seitz.targets", native)
    shutil.copy2(ROOT / "README.md", args.stage)
    shutil.copy2(ROOT / "LICENSE", args.stage)
    # Headers only: lib/cmake and share/ describe a find_package() consumer,
    # which is not who this package is for.
    shutil.copytree(args.release / "include", native / "include")

    ok = all([merge(prefix, native / "lib" / "x64" / config / "seitz.lib")
              for config, prefix in (("Release", args.release), ("Debug", args.debug))])
    if not ok:
        return 1
    print(f"staged {args.stage}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
