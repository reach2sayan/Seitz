#!/usr/bin/env python3
"""Transcribe spglib's msg_database.c data tables into constexpr C++ headers.

The magnetic space-group database (1651 UNI numbers) is pure generated data,
mirroring spg_database.c. We mechanically convert the six C arrays
  - magnetic_spacegroup_types[1652]              metadata
  - magnetic_spacegroup_hall_mapping[531][2]     {lo UNI, hi UNI} per Hall #
  - magnetic_spacegroup_uni_mapping[1652][2]     {#Hall settings, first Hall #}
  - magnetic_spacegroup_operation_index[1652][18][2]  {count, start} per setting
  - magnetic_symmetry_operations[]               encoded ops (timerev + spgdb)
  - alternative_transformations[1652][18][7]     encoded std transforms
into `inline constexpr std::array`s. The decoder/accessors are hand-written in
src/data/msg_database.cpp.

As with transcribe_spg_database.py, the output is split into two headers so the
encoded-operation packing (compile-time-only compaction) never reaches a public
header:
  - <metadata_out>: the catalog metadata + the two small mapping tables.
  - <ops_out>: the encoded ops, the operation index, and the alternative
    transformations; included only by src/data/msg_database.cpp.

Usage: transcribe_msg_database.py <msg_database.c> <metadata_out> <ops_out>
"""
import re
import sys


def extract_array_body(text, decl):
    """Return the text between the outermost braces of `decl ... = { ... }`."""
    start = text.index(decl)
    brace = text.index("{", start)
    depth, i, in_str = 0, brace, False
    while i < len(text):
        c = text[i]
        if in_str:
            if c == "\\":
                i += 2
                continue
            if c == '"':
                in_str = False
        elif c == '"':
            in_str = True
        elif c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return text[brace + 1:i]
        i += 1
    raise ValueError("unbalanced braces for " + decl)


def strip_comments(s):
    s = re.sub(r"/\*.*?\*/", "", s, flags=re.DOTALL)
    s = re.sub(r"//[^\n]*", "", s)
    return s


def parse_braced(s, i):
    """Recursive-descent parse of a C brace initializer of ints into nested
    python lists. s[i] must be '{'. Returns (list, next_index)."""
    assert s[i] == "{", s[i:i + 20]
    i += 1
    out = []
    while True:
        while i < len(s) and s[i] in " \t\r\n,":
            i += 1
        if s[i] == "}":
            return out, i + 1
        if s[i] == "{":
            sub, i = parse_braced(s, i)
            out.append(sub)
        else:
            j = i
            while j < len(s) and s[j] not in "{}, \t\r\n":
                j += 1
            out.append(int(s[i:j]))
            i = j


def parse_int_array(text, decl):
    """Parse a (possibly nested) int array body into nested python lists."""
    body = strip_comments(extract_array_body(text, decl))
    out, _ = parse_braced("{" + body + "}", 0)
    return out


def split_top_level(s):
    """Split on top-level commas, respecting strings and nested braces."""
    parts, depth, in_str, cur, i = [], 0, False, [], 0
    while i < len(s):
        c = s[i]
        if in_str:
            cur.append(c)
            if c == "\\":
                cur.append(s[i + 1])
                i += 2
                continue
            if c == '"':
                in_str = False
        elif c == '"':
            in_str = True
            cur.append(c)
        elif c in "{[":
            depth += 1
            cur.append(c)
        elif c in "}]":
            depth -= 1
            cur.append(c)
        elif c == "," and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(c)
        i += 1
    if "".join(cur).strip():
        parts.append("".join(cur))
    return parts


def parse_entries(body):
    """Yield each top-level {...} group's inner text (preserving strings)."""
    depth, in_str, cur, started, i = 0, False, [], False, 0
    while i < len(body):
        c = body[i]
        if in_str:
            cur.append(c)
            if c == "\\":
                cur.append(body[i + 1]); i += 2; continue
            if c == '"':
                in_str = False
        elif c == '"':
            in_str = True; cur.append(c)
        elif c == "{":
            depth += 1
            if depth == 1:
                started = True; cur = []
            else:
                cur.append(c)
        elif c == "}":
            depth -= 1
            if depth == 0 and started:
                yield "".join(cur); started = False
            else:
                cur.append(c)
        elif started:
            cur.append(c)
        i += 1


def parse_string_literal(tok):
    tok = tok.strip()
    return tok[tok.index('"') + 1:tok.rindex('"')]


def cpp_string(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def pad(rows, n, fill):
    rows = [list(r) for r in rows]
    assert len(rows) <= n, (len(rows), n)
    return rows + [list(fill) for _ in range(n - len(rows))]


def main():
    src, meta_out, ops_out = sys.argv[1], sys.argv[2], sys.argv[3]
    text = open(src).read()

    # --- metadata: magnetic_spacegroup_types[1652] ---
    types_body = extract_array_body(text, "magnetic_spacegroup_types[]")
    type_entries = list(parse_entries(strip_comments(types_body)))
    types = []
    for e in type_entries:
        f = split_top_level(e)
        assert len(f) == 6, (len(f), e)
        uni, litvin = int(f[0]), int(f[1])
        bns, og = parse_string_literal(f[2]), parse_string_literal(f[3])
        number, typ = int(f[4]), int(f[5])
        types.append((uni, litvin, bns, og, number, typ))
    assert len(types) == 1652, len(types)

    hall_mapping = parse_int_array(text, "magnetic_spacegroup_hall_mapping[][2]")
    uni_mapping = parse_int_array(text, "magnetic_spacegroup_uni_mapping[][2]")
    assert len(hall_mapping) == 531, len(hall_mapping)
    assert len(uni_mapping) == 1652, len(uni_mapping)

    op_index = parse_int_array(
        text, "magnetic_spacegroup_operation_index[][18][2]")
    op_index = [pad(e, 18, [0, 0]) for e in op_index]
    assert len(op_index) == 1652, len(op_index)

    ops = parse_int_array(text, "magnetic_symmetry_operations[]")
    assert all(isinstance(x, int) for x in ops)

    alt = parse_int_array(text, "alternative_transformations[][18][7]")
    alt = [pad([r + [0] * (7 - len(r)) for r in e], 18, [0] * 7) for e in alt]
    assert len(alt) == 1652, len(alt)

    gen_banner = (
        "// GENERATED by tools/transcribe_msg_database.py from spglib v2.7.0\n"
        "// msg_database.c. Do not edit by hand.\n")

    # --- metadata header (public): catalog metadata + small mapping tables ---
    with open(meta_out, "w") as fh:
        w = fh.write
        w("#pragma once\n\n")
        w(gen_banner)
        w("//\n")
        w("// Magnetic space-group metadata + the Hall<->UNI mapping tables (no\n")
        w("// encoded operations) so they can back the constexpr catalog in the\n")
        w("// public msg_database.hpp without leaking the operation encoding.\n\n")
        w("#include <array>\n\n")
        w("namespace spglib::data {\n\n")
        w("struct MagneticSpacegroupTypeRaw {\n")
        w("  int uni_number;\n")
        w("  int litvin_number;\n")
        w("  char const *bns_number;\n")
        w("  char const *og_number;\n")
        w("  int number;\n")
        w("  int type;\n")
        w("};\n\n")
        w("inline constexpr std::array<MagneticSpacegroupTypeRaw, %d>\n"
          "    kMagneticSpacegroupTypes = {{\n" % len(types))
        for (uni, litvin, bns, og, number, typ) in types:
            w("    {%d, %d, %s, %s, %d, %d},\n"
              % (uni, litvin, cpp_string(bns), cpp_string(og), number, typ))
        w("}};\n\n")

        w("// Per Hall number (0..530): {smallest UNI number, largest UNI number}\n")
        w("// of the magnetic space groups whose family space group is that Hall #.\n")
        w("inline constexpr std::array<std::array<int, 2>, %d>\n"
          "    kMagneticHallMapping = {{\n" % len(hall_mapping))
        for (lo, hi) in hall_mapping:
            w("    {{%d, %d}},\n" % (lo, hi))
        w("}};\n\n")

        w("// Per UNI number (0..1651): {number of Hall settings, smallest Hall #}.\n")
        w("inline constexpr std::array<std::array<int, 2>, %d>\n"
          "    kMagneticUniMapping = {{\n" % len(uni_mapping))
        for (n, first) in uni_mapping:
            w("    {{%d, %d}},\n" % (n, first))
        w("}};\n\n")
        w("} // namespace spglib::data\n")

    # --- operations header (private): encoded ops + index + alt transforms ---
    with open(ops_out, "w") as fh:
        w = fh.write
        w("#pragma once\n\n")
        w(gen_banner)
        w("//\n")
        w("// Encoded magnetic operations (compile-time-only compaction):\n")
        w("//   enc = timerev * 34012224 + spgdb_encoded_symmetry,\n")
        w("// where 34012224 = 3^9 * 12^3. Decoded in src/data/msg_database.cpp;\n")
        w("// must NOT be included by a public header.\n\n")
        w("#include <array>\n\n")
        w("namespace spglib::data {\n\n")

        w("inline constexpr std::array<int, %d> kMagneticOperations = {{\n"
          % len(ops))
        for i in range(0, len(ops), 8):
            w("    " + ", ".join(str(x) for x in ops[i:i + 8]) + ",\n")
        w("}};\n\n")

        w("// Per [UNI number][Hall-setting offset]: {operation count, start\n")
        w("// index into kMagneticOperations}. Rows padded to 18 with {0, 0}.\n")
        w("inline constexpr std::array<std::array<std::array<int, 2>, 18>, %d>\n"
          "    kMagneticOperationIndex = {{\n" % len(op_index))
        for e in op_index:
            w("    {{" + ", ".join("{{%d, %d}}" % (p[0], p[1]) for p in e)
              + "}},\n")
        w("}};\n\n")

        w("// Per [UNI number][Hall-setting offset]: up to 7 encoded alternative\n")
        w("// transformations (spgdb encoding), zero-terminated; rows padded to 7\n")
        w("// and settings padded to 18 with zeros.\n")
        w("inline constexpr std::array<std::array<std::array<int, 7>, 18>, %d>\n"
          "    kAlternativeTransformations = {{\n" % len(alt))
        for e in alt:
            w("    {{" + ", ".join("{{%s}}" % ", ".join(str(x) for x in r)
                                   for r in e) + "}},\n")
        w("}};\n\n")
        w("} // namespace spglib::data\n")

    print("types: %d, hall_mapping: %d, uni_mapping: %d, op_index: %d, "
          "ops: %d, alt: %d" % (len(types), len(hall_mapping),
                                len(uni_mapping), len(op_index), len(ops),
                                len(alt)))


if __name__ == "__main__":
    main()
