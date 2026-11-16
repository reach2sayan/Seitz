#pragma once

// The hand-built crystals the suites share: one definition each, so a test file
// names the structure it exercises instead of carrying its own copy of it.

#include "helpers.hpp"

#include <seitz/core/cell.hpp>

#include <cmath>

namespace seitz::test {

// ---- lattices (columns are the basis vectors) ----

[[nodiscard]] inline Matrix3d columns(Vector3d const &a, Vector3d const &b,
                                      Vector3d const &c) {
  Matrix3d l;
  l << a, b, c;
  return l;
}

[[nodiscard]] inline Matrix3d cubic(double a) {
  return a * Matrix3d::Identity();
}

[[nodiscard]] inline Matrix3d orthorhombic(double a, double b, double c) {
  return Matrix3d{Vector3d{a, b, c}.asDiagonal()};
}

[[nodiscard]] inline Matrix3d tetragonal(double a, double c) {
  return orthorhombic(a, a, c);
}

// a = b, gamma = 120; also the layer setting when c is a vacuum gap.
[[nodiscard]] inline Matrix3d hexagonal(double a, double c) {
  return columns(Vector3d(a, 0.0, 0.0),
                 Vector3d(-a / 2.0, a * std::sqrt(3.0) / 2.0, 0.0),
                 Vector3d(0.0, 0.0, c));
}

// ---- structures, with the space group each one exercises ----

[[nodiscard]] inline Cell primitive_cubic(double a) { // Pm-3m
  return make_cell(a, {{0.0, 0.0, 0.0}}, {0});
}

[[nodiscard]] inline Cell bcc(double a) { // Im-3m
  return make_cell(a, {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, {0, 0});
}

[[nodiscard]] inline Cell fcc(double a) { // Fm-3m
  return make_cell(
      a, {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}, {0.5, 0.0, 0.5}, {0.0, 0.5, 0.5}},
      {0, 0, 0, 0});
}

// Conventional NaCl: 4 Na + 4 Cl (Fm-3m), non-primitive.
[[nodiscard]] inline Cell rock_salt(double a) {
  return make_cell(a,
                   {{0.0, 0.0, 0.0},
                    {0.5, 0.5, 0.0},
                    {0.5, 0.0, 0.5},
                    {0.0, 0.5, 0.5},
                    {0.5, 0.5, 0.5},
                    {0.0, 0.0, 0.5},
                    {0.0, 0.5, 0.0},
                    {0.5, 0.0, 0.0}},
                   {0, 0, 0, 0, 1, 1, 1, 1});
}

// Conventional NaCl with real species (Na 11, Cl 17), the atom order the
// standardization tests were written against.
[[nodiscard]] inline Cell nacl(double a = 5.64) {
  return make_cell(a,
                   {{0, 0, 0},
                    {0, 0.5, 0.5},
                    {0.5, 0, 0.5},
                    {0.5, 0.5, 0},
                    {0.5, 0.5, 0.5},
                    {0.5, 0, 0},
                    {0, 0.5, 0},
                    {0, 0, 0.5}},
                   {11, 11, 11, 11, 17, 17, 17, 17});
}

[[nodiscard]] inline Cell diamond(double a) { // Fd-3m (non-symmorphic)
  return make_cell(a,
                   {{0.0, 0.0, 0.0},
                    {0.0, 0.5, 0.5},
                    {0.5, 0.0, 0.5},
                    {0.5, 0.5, 0.0},
                    {0.25, 0.25, 0.25},
                    {0.25, 0.75, 0.75},
                    {0.75, 0.25, 0.75},
                    {0.75, 0.75, 0.25}},
                   {0, 0, 0, 0, 0, 0, 0, 0});
}

[[nodiscard]] inline Cell zincblende(double a) { // F-43m
  Cell const d = diamond(a);
  return Cell(d.lattice(), d.positions(), {0, 0, 0, 0, 1, 1, 1, 1});
}

[[nodiscard]] inline Cell bcc_tetragonal() { // I4/mmm
  return make_cell(tetragonal(4.0, 6.0), {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}},
                   {0, 0});
}

[[nodiscard]] inline Cell rutile() { // P4_2/mnm (non-symmorphic)
  return make_cell(tetragonal(4.6, 3.0),
                   {{0.0, 0.0, 0.0},
                    {0.5, 0.5, 0.5},
                    {0.3, 0.3, 0.0},
                    {0.7, 0.7, 0.0},
                    {0.8, 0.2, 0.5},
                    {0.2, 0.8, 0.5}},
                   {0, 0, 1, 1, 1, 1});
}

[[nodiscard]] inline Cell orthorhombic_p() { // Pmmm
  return make_cell(orthorhombic(4.0, 5.0, 6.0), {{0.0, 0.0, 0.0}}, {0});
}

[[nodiscard]] inline Cell c_centered_ortho() { // Cmmm
  return make_cell(orthorhombic(4.0, 6.0, 7.0),
                   {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}}, {0, 0});
}

// Unique axis b, beta != 90: exercises the monoclinic matcher.
[[nodiscard]] inline Cell monoclinic_p() {
  return make_cell(columns(Vector3d(5.0, 0.0, 0.0), Vector3d(0.0, 6.0, 0.0),
                           Vector3d(1.5, 0.0, 7.0)),
                   {{0.0, 0.0, 0.0}}, {0});
}

[[nodiscard]] inline Cell c_centered_monocli() { // C2/m, unique b
  double const beta = 1.91; // radians, ~109.5 deg
  return make_cell(columns(Vector3d(5.0, 0.0, 0.0), Vector3d(0.0, 6.0, 0.0),
                           Vector3d(7.0 * std::cos(beta), 0.0,
                                    7.0 * std::sin(beta))),
                   {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}}, {0, 0});
}

[[nodiscard]] inline Cell hcp() { // P6_3/mmc
  return make_cell(hexagonal(3.2, 5.2),
                   {{1.0 / 3.0, 2.0 / 3.0, 0.25}, {2.0 / 3.0, 1.0 / 3.0, 0.75}},
                   {0, 0});
}

[[nodiscard]] inline Cell wurtzite() { // P6_3mc (polar)
  double const u = 0.375;
  return make_cell(hexagonal(3.2, 5.2),
                   {{1.0 / 3.0, 2.0 / 3.0, 0.0},
                    {2.0 / 3.0, 1.0 / 3.0, 0.5},
                    {1.0 / 3.0, 2.0 / 3.0, u},
                    {2.0 / 3.0, 1.0 / 3.0, 0.5 + u}},
                   {0, 0, 1, 1});
}

// R-centred hexagonal cell (trigonal/rhombohedral).
[[nodiscard]] inline Cell rhombohedral_r() {
  return make_cell(hexagonal(3.0, 12.0),
                   {{0.0, 0.0, 0.0},
                    {1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0},
                    {2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}},
                   {0, 0, 0});
}

// P1: two unlike atoms in general positions.
[[nodiscard]] inline Cell triclinic() {
  return make_cell(columns(Vector3d(4.0, 0.0, 0.0), Vector3d(0.7, 4.3, 0.0),
                           Vector3d(0.5, 0.9, 5.1)),
                   {{0.1, 0.2, 0.3}, {0.6, 0.55, 0.27}}, {0, 1});
}

} // namespace seitz::test
