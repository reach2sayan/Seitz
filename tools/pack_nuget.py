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


def headers(prefix, out):
    """The prefix's include/ as one flat tree, packable by `nuget pack`."""
    shutil.copytree(prefix / "include", out)

    # Boost's CMake install is versioned: include/boost-1_88/boost/... . One
    # include directory in Seitz.targets is worth more than a path that has to
    # be rewritten at every Boost bump, so the version dir is unwrapped here.
    versioned = [d for d in out.glob("boost-*") if d.is_dir()]
    if len(versioned) > 1:
        print(f"more than one versioned Boost include dir: {versioned}", file=sys.stderr)
        return False
    for d in versioned:
        for child in d.iterdir():
            shutil.move(str(child), out / child.name)
        d.rmdir()

    # `nuget pack` excludes dotfiles by default and can then die -- "String
    # cannot be empty. Parameter name: entryName" -- on the directory it just
    # emptied. Boost ships a boost/headers/.gitkeep that lands exactly there.
    # Drop what NuGet would drop, so it has nothing left to empty. (The pack
    # step also passes -NoDefaultExcludes, so this is what is actually removed,
    # rather than two rules disagreeing about it.)
    for path in out.rglob(".*"):
        if path.is_file():
            path.unlink()
    return True


def prune(stage):
    """Empty directories, deepest first. NuGet packs a directory as a zero-length
    entry name and dies on it, so none may reach `nuget pack` -- and one can
    appear at any depth once the dotfiles above are gone."""
    for path in sorted(stage.rglob("*"), key=lambda p: -len(p.parts)):
        if path.is_dir() and not any(path.iterdir()):
            path.rmdir()
    left = [p for p in stage.rglob("*") if p.is_dir() and not any(p.iterdir())]
    if left:
        print(f"empty directories survived pruning: {left}", file=sys.stderr)
        return False
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
    # LICENSE.txt, not LICENSE: NuGet decides whether a <file> target names a
    # file or a folder by looking for an extension, and an extensionless entry
    # at the package root is one of the shapes that reaches the entryName crash.
    # The nuspec's license is an SPDX expression either way; this copy is for
    # someone reading the package.
    shutil.copy2(ROOT / "LICENSE", args.stage / "LICENSE.txt")
    # Headers only: lib/cmake and share/ describe a find_package() consumer,
    # which is not who this package is for.
    ok = all([headers(args.release, native / "include"),
              *(merge(prefix, native / "lib" / "x64" / config / "seitz.lib")
                for config, prefix in (("Release", args.release), ("Debug", args.debug)))])
    ok = prune(args.stage) and ok
    if not ok:
        return 1
    print(f"staged {args.stage}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
