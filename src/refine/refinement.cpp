#include "refine/refinement.hpp"

#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/lattice.hpp>
#include <cppcrystal/core/point_group.hpp>
#include <cppcrystal/data/spg_database.hpp>
#include "math/fractional.hpp"
#include "math/integer_matrix.hpp"
#include "symmetry/pointgroup.hpp"

#include <cmath>
#include <ranges>

namespace cppcrystal::refine {

using data::operations_from_database;

namespace {

// ---- idealized conventional lattice setters --------------------------------
// Each fills `m` (columns = basis vectors) from the metric's lengths/angles.

[[nodiscard]] double len(Matrix3d const &g, int i) {
  return std::sqrt(g(i, i));
}

[[nodiscard]] Matrix3d set_tricli(Matrix3d const &g) {
  double const a = len(g, 0), b = len(g, 1), c = len(g, 2);
  double const alpha = std::acos(g(1, 2) / b / c);
  double const beta = std::acos(g(0, 2) / a / c);
  double const gamma = std::acos(g(0, 1) / a / b);
  double const cg = std::cos(gamma), cb = std::cos(beta), ca = std::cos(alpha);
  double const sg = std::sin(gamma);

  Matrix3d m = Matrix3d::Zero();
  m(0, 0) = a;
  m(0, 1) = b * cg;
  m(0, 2) = c * cb;
  m(1, 1) = b * sg;
  m(1, 2) = c * (ca - cb * cg) / sg;
  m(2, 2) =
      c * std::sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg) / sg;
  return m;
}

[[nodiscard]] Matrix3d set_monocli(Matrix3d const &g, std::string_view choice) {
  // choice is e.g. "b", "-b", "b1", "c", "a"; the axis letter selects the
  // unique axis (the leading '-' is a sign convention, skipped).
  std::size_t const pos = (!choice.empty() && choice[0] == '-') ? 1 : 0;
  char const axis = pos < choice.size() ? choice[pos] : 'b';
  double const a = len(g, 0), b = len(g, 1), c = len(g, 2);

  Matrix3d m = Matrix3d::Zero();
  switch (axis) {
  case 'a': {
    double const angle = std::acos(g(1, 2) / b / c);
    m(0, 2) = c;
    m(1, 0) = a;
    m(0, 1) = b * std::cos(angle);
    m(2, 1) = b * std::sin(angle);
    break;
  }
  case 'b': {
    double const angle = std::acos(g(0, 2) / a / c);
    m(0, 0) = a;
    m(1, 1) = b;
    m(0, 2) = c * std::cos(angle);
    m(2, 2) = c * std::sin(angle);
    break;
  }
  case 'c': {
    double const angle = std::acos(g(0, 1) / a / b);
    m(0, 1) = b;
    m(1, 2) = c;
    m(0, 0) = a * std::cos(angle);
    m(2, 0) = a * std::sin(angle);
    break;
  }
  default:
    break;
  }
  return m;
}

[[nodiscard]] Matrix3d set_ortho(Matrix3d const &g) {
  Matrix3d m = Matrix3d::Zero();
  m(0, 0) = len(g, 0);
  m(1, 1) = len(g, 1);
  m(2, 2) = len(g, 2);
  return m;
}

[[nodiscard]] Matrix3d set_tetra(Matrix3d const &g) {
  double const ab = (len(g, 0) + len(g, 1)) / 2;
  Matrix3d m = Matrix3d::Zero();
  m(0, 0) = ab;
  m(1, 1) = ab;
  m(2, 2) = len(g, 2);
  return m;
}

[[nodiscard]] Matrix3d set_trigo(Matrix3d const &g) {
  double const ab = len(g, 0) + len(g, 1);
  Matrix3d m = Matrix3d::Zero();
  m(0, 0) = ab / 2;
  m(0, 1) = -ab / 4;
  m(1, 1) = ab / 4 * std::sqrt(3.0);
  m(2, 2) = len(g, 2);
  return m;
}

[[nodiscard]] Matrix3d set_rhomb(Matrix3d const &g) {
  double const a = len(g, 0), b = len(g, 1), c = len(g, 2);
  double const angle =
      std::acos((g(0, 1) / a / b + g(0, 2) / a / c + g(1, 2) / b / c) / 3);
  double const ahex = 2 * (a + b + c) / 3 * std::sin(angle / 2);
  double const chex = (a + b + c) / 3 * std::sqrt(3 * (1 + 2 * std::cos(angle)));

  Matrix3d m = Matrix3d::Zero();
  m(0, 0) = ahex / 2;
  m(1, 0) = ahex / (2 * std::sqrt(3.0));
  m(2, 0) = chex / 3;
  m(0, 1) = -ahex / 2;
  m(1, 1) = ahex / (2 * std::sqrt(3.0));
  m(2, 1) = chex / 3;
  m(0, 2) = 0;
  m(1, 2) = -ahex / std::sqrt(3.0);
  m(2, 2) = chex / 3;
  return m;
}

[[nodiscard]] Matrix3d set_cubic(Matrix3d const &g) {
  double const a = (len(g, 0) + len(g, 1) + len(g, 2)) / 3;
  Matrix3d m = Matrix3d::Zero();
  m(0, 0) = a;
  m(1, 1) = a;
  m(2, 2) = a;
  return m;
}

} // namespace

Lattice conventional_lattice(spacegroup::Spacegroup const &sg) {
  Matrix3d const g = Lattice{sg.bravais_lattice}.metric();
  Holohedry const holohedry =
      symmetry::pointgroup_by_number(sg.type.pointgroup_number).holohedry;

  return Lattice{[&]() -> Matrix3d {
    switch (holohedry) {
    case Holohedry::triclinic:
      return set_tricli(g);
    case Holohedry::monoclinic:
      return set_monocli(g, sg.type.choice);
    case Holohedry::orthorhombic:
      return set_ortho(g);
    case Holohedry::tetragonal:
      return set_tetra(g);
    case Holohedry::trigonal:
      return (!sg.type.choice.empty() && sg.type.choice[0] == 'R')
                 ? set_rhomb(g)
                 : set_trigo(g);
    case Holohedry::hexagonal:
      return set_trigo(g);
    case Holohedry::cubic:
      return set_cubic(g);
    default:
      return Matrix3d::Zero();
    }
  }()};
}

spacegroup::Spacegroup find_similar_bravais_lattice(spacegroup::Spacegroup sg,
                                                    double symprec) {
  Operations const conv_sym = operations_from_database(sg.type.hall_number);
  Matrix3d const std_lattice = conventional_lattice(sg).matrix();

  // Find the proper-rotation setting closest (Frobenius) to the idealized one.
  double min_length = sg.bravais_lattice.norm();
  std::optional<std::size_t> rot_i = std::nullopt;
  Matrix3d rot_lat = sg.bravais_lattice;
  for (const auto& [i, sym] : conv_sym | std::views::enumerate) {
    Matrix3i const &rot = sym.rotation;
    if (rot.determinant() < 0) {
      continue;
    }
    Matrix3d const candidate = sg.bravais_lattice * rot.cast<double>();
    double const length = (candidate - std_lattice).norm();
    if (length < min_length - symprec) {
      rot_lat = candidate;
      rot_i = i;
      min_length = length;
    }
  }

  if (!rot_i.has_value()) {
    return sg;
  }

  // Adjust the origin shift to the shortest equivalent for the chosen rotation
  // (3D path, lattice_rank = 3).
  Matrix3i const &chosen = conv_sym[rot_i.value()].rotation;
  double shortest_length = 2.0;
  Vector3d shortest_p = sg.origin_shift;
  for (auto const &op :
       conv_sym | std::views::filter([&chosen](auto const &op) {
         return op.rotation == chosen;
       })) {
    Matrix3d const inv = op.rotation.cast<double>().inverse();
    Vector3d p = inv * sg.origin_shift - inv * op.translation;
    p = math::nearest_offset(p); // p[j] -= nint(p[j]) for all 3 axes
    double const length = p.norm();
    if (length < shortest_length - symprec) {
      shortest_length = length;
      shortest_p = math::wrap_to_unit_cell(p);
    }
  }

  sg.origin_shift = shortest_p;
  sg.bravais_lattice = rot_lat;
  return sg;
}

} // namespace cppcrystal::refine
