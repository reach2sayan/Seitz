#!/usr/bin/env python3
"""Diff the vendored mdspan against upstream, ignoring formatting.

Usage:
    python3 mdspan_vendor_diff.py                  # against the pinned commit
    python3 mdspan_vendor_diff.py --ref single-header   # against the branch tip
    python3 mdspan_vendor_diff.py --find-base      # which upstream commit is it?

Expected output today: eight tokens.  Six are formatting (a sorted #include, a
pair of swapped using-declarations); two are the LOCAL PATCH guard that keeps
at() out of a -fno-exceptions build.
"""
import argparse
import json
import re
import sys
import urllib.request
from difflib import SequenceMatcher
from pathlib import Path

from pygments.lexers import CppLexer
from pygments.token import Comment, Number, String, Text

REPO = "kokkos/mdspan"
BRANCH = "single-header"
# The commit include/seitz/third_party/mdspan.hpp was taken from.  Bump on re-vendor.
PINNED = "c3fc07b607db2a43cd65aabc12ffda9328833d95"
VENDORED = Path(__file__).resolve().parent.parent / "include/seitz/third_party/mdspan.hpp"

LEXER = CppLexer()

# Pygments lexes a whole `#if 0` block as one comment, which would hide a local
# patch inside one -- and the vendored header has two of them. `#if (0)` means
# the same and does not trip that rule. Both sides get it, so it is never itself
# a difference.
IF_ZERO = re.compile(r"^#if\s+0\s*$", re.M)


def normalize(src: str) -> list[str]:
    """One token per line, comments and whitespace dropped.

    A directive stays whole: it ends at its newline, so that one line break is
    meaning where every other one is formatting. Adjacent string literals are
    joined, because "abc " "def" and "abc def" are the same string to the
    compiler and clang-format moves the split around freely.

    Most of what follows is undoing Pygments' own line-sensitivity, since the
    whole point here is to be blind to formatting: it hands back a directive in
    raw chunks (re-lexed below, so spacing inside one cannot show up as a
    difference), a string literal in pieces, and a sign folded into the number
    after it -- but only when nothing separates them, so `x-1` and `x - 1`
    would differ.
    """
    out: list[str] = []
    directive: list[str] | None = None
    literal: str | None = None

    def flush() -> None:
        nonlocal directive, literal
        if directive is not None: out.append(" ".join(["#DIRECTIVE", *directive]))
        if literal is not None: out.append(literal)
        directive = literal = None

    for kind, value in LEXER.get_tokens(IF_ZERO.sub("#if (0)", src).replace("\\\n", "")):
        if kind in Comment.Preproc or kind in Comment.PreprocFile:
            # The leading `#` is dropped, which is also what keeps the re-lex
            # from finding a directive all over again.
            if value == "\n": flush()
            elif value != "#": directive = (directive or []) + [
                v for k, v in LEXER.get_tokens(value)
                if k not in Text and k not in Comment and v.strip()]
        elif kind in Comment:
            # A `//` comment ends the directive it sits on: Pygments leaves the
            # state without the newline that would otherwise close it.
            if "\n" in value: flush()
        elif kind in Text:
            pass                            # never ends a run of literals
        elif kind in String:
            literal = value if literal is None else (
                literal[:-1] + value[1:] if literal.endswith('"') and value.startswith('"')
                else literal + value)
        else:
            flush()
            out += [value[:1], value[1:]] if kind in Number and value[:1] in "+-" else [value]
    flush()
    return out


def diff(a: list[str], b: list[str], context: int = 3) -> list[str]:
    """Unified diff over token lists.  difflib.unified_diff would do, but its
    autojunk heuristic treats ';' and '::' -- thousands of them here -- as junk
    and stops finding the minimal edit; SequenceMatcher with autojunk off does."""
    ops = SequenceMatcher(None, a, b, autojunk=False).get_opcodes()
    out: list[str] = []
    for i, (tag, i1, i2, j1, j2) in enumerate(ops):
        if tag == "equal":
            head = i == 0
            tail = i == len(ops) - 1
            if i2 - i1 <= 2 * context and not head and not tail:
                out += [" " + l for l in a[i1:i2]]
            else:
                if not head:
                    out += [" " + l for l in a[i1 : i1 + context]]
                if not tail:
                    out.append(f"@@ {i2 - context + 1} @@")
                    out += [" " + l for l in a[max(i1, i2 - context) : i2]]
            continue
        out += ["-" + l for l in a[i1:i2]]
        out += ["+" + l for l in b[j1:j2]]
    return out

def fetch(ref: str, cache: Path) -> str:
    path = cache / f"{ref}.hpp"
    if not path.exists():
        url = f"https://raw.githubusercontent.com/{REPO}/{ref}/mdspan.hpp"
        cache.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(url, timeout=60) as r:
            path.write_bytes(r.read())
    return path.read_text()

def commits(limit: int) -> list[tuple[str, str]]:
    url = f"https://api.github.com/repos/{REPO}/commits?sha={BRANCH}&per_page={limit}"
    with urllib.request.urlopen(url, timeout=60) as r:
        return [(c["sha"], c["commit"]["committer"]["date"][:10]) for c in json.load(r)]

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", default=PINNED, help="upstream commit or branch (default: the pinned commit)")
    ap.add_argument("--file", type=Path, default=VENDORED, help="the vendored copy to check")
    ap.add_argument("--find-base", type=int, nargs="?", const=40, metavar="N",
                    help="rank the last N upstream commits by how close they are")
    ap.add_argument("--cache-dir", type=Path, default=Path("/tmp/mdspan-upstream"))
    args = ap.parse_args()

    ours = normalize(args.file.read_text())
    if args.find_base:
        ranked = []
        for sha, date in commits(args.find_base):
            n = sum(1 for l in diff(normalize(fetch(sha, args.cache_dir)), ours, context=0)
                    if l[:1] in "-+")
            ranked.append((n, date, sha))
            print(f"  {n:6d} tokens  {date}  {sha}", file=sys.stderr)
        n, date, sha = min(ranked)
        print(f"\nclosest: {sha} ({date}), {n} differing tokens")
        return 0

    theirs = normalize(fetch(args.ref, args.cache_dir))
    hunks = diff(theirs, ours)
    changed = sum(1 for l in hunks if l[:1] in "-+")
    print(f"--- upstream {args.ref[:10]}\n+++ {args.file}")
    print("\n".join(hunks))
    print(f"--- {len(theirs)} upstream tokens, {len(ours)} ours, {changed} differing")
    return 1 if changed else 0

if __name__ == "__main__":
    sys.exit(main())
