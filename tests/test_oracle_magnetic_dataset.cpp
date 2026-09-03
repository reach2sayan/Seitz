// Oracle test for the full magnetic dataset (spglib.c get_magnetic_dataset):
// cppcrystal::get_magnetic_dataset must reproduce spg_get_magnetic_dataset — the
// magnetic space-group identity, the standardized cell (lattice + positions +
// site tensors), equivalent atoms, and the transformation to the standardized
// setting.

#include "oracle.hpp"

#include <cppcrystal/core/magnetic_cell.hpp>
#include <cppcrystal/magnetic_dataset.hpp>

#include <cppcrystal/math/integer_matrix.hpp>

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <vector>

namespace {

using cppcrystal::Cell;
using cppcrystal::Lattice;
using cppcrystal::CollinearTensors;
using cppcrystal::MagneticCell;
using cppcrystal::Matrix3d;
using cppcrystal::Matrix3i;
using cppcrystal::NoncollinearTensors;
using cppcrystal::noncollinear_tensors;
using cppcrystal::Positions;
using cppcrystal::SiteTensor;
using cppcrystal::TensorKind;
using cppcrystal::SiteTensors;
using cppcrystal::Vector3d;

Cell make_cell(double a, std::vector<std::array<double, 3>> const &pos,
               std::vector<int> const &types) {
  Positions p(static_cast<Eigen::Index>(pos.size()), 3);
  for (std::size_t i = 0; i < pos.size(); ++i) {
    p.row(static_cast<Eigen::Index>(i)) =
        Eigen::RowVector3d(pos[i][0], pos[i][1], pos[i][2]);
  }
  return Cell(Lattice{a * Matrix3d::Identity()}, p, types);
}

// Flat tensor array (spglib layout) from a MagneticCell.
std::vector<double> flat_tensors(MagneticCell const &mcell) {
  std::vector<double> out;
  Eigen::Index const n = mcell.size();
  if (mcell.rank() == SiteTensor::collinear) {
    for (Eigen::Index i = 0; i < n; ++i)
      out.push_back(mcell.scalar(i));
  } else {
    for (Eigen::Index i = 0; i < n; ++i)
      for (int s = 0; s < 3; ++s)
        out.push_back(mcell.vector(i)[s]);
  }
  return out;
}

bool fractional_overlap(Vector3d const &a, Vector3d const &b, double symprec) {
  Vector3d d = a - b;
  d = d.unaryExpr([](double x) { return x - std::round(x); });
  return d.cwiseAbs().maxCoeff() <= symprec;
}

bool same_lattice(Matrix3d const &a, Matrix3d const &b) {
  Matrix3d const rel = a.inverse() * b;
  Matrix3i const rounded = cppcrystal::math::round_to_int(rel);
  return (rel - rounded.cast<double>()).cwiseAbs().maxCoeff() < 1e-6 &&
         std::abs(rounded.determinant()) == 1;
}

void check(MagneticCell const &input, TensorKind kind, double symprec) {
  bool const is_axial = kind == TensorKind::axial;
  MagneticCell const mcell(input.cell(), input.tensors(), kind);
  auto const got = cppcrystal::get_magnetic_dataset(mcell, {symprec});
  REQUIRE(got);

  // Reference dataset.
  cppcrystal::oracle::CCell c(mcell.cell());
  std::vector<double> const tensors = flat_tensors(mcell);
  int const rank = mcell.rank() == SiteTensor::collinear ? 0 : 1;
  SpglibMagneticDataset *ref = spg_get_magnetic_dataset(
      c.lattice, c.pos(), c.types.data(), tensors.data(), rank, c.num_atom(),
      is_axial ? 1 : 0, symprec);
  REQUIRE(ref != nullptr);

  // Magnetic space-group identity.
  CHECK(got->uni_number == ref->uni_number);
  CHECK(static_cast<int>(got->msg_type) == ref->msg_type);
  CHECK(got->hall_number == ref->hall_number);
  CHECK(static_cast<int>(got->std_types.size()) == ref->n_std_atoms);

  // Transformation to the standardized setting.
  Matrix3d ref_tmat;
  cppcrystal::oracle::from_c_lattice(ref_tmat, ref->transformation_matrix);
  CHECK((got->transformation_matrix - ref_tmat).cwiseAbs().maxCoeff() < 1e-5);
  Vector3d const ref_shift(ref->origin_shift[0], ref->origin_shift[1],
                           ref->origin_shift[2]);
  CHECK((got->origin_shift - ref_shift).cwiseAbs().maxCoeff() < 1e-5);
  Matrix3d ref_std_rot;
  cppcrystal::oracle::from_c_lattice(ref_std_rot, ref->std_rotation_matrix);
  CHECK((got->std_rotation_matrix - ref_std_rot).cwiseAbs().maxCoeff() < 1e-5);

  // Standardized lattice.
  Matrix3d ref_std_lat;
  cppcrystal::oracle::from_c_lattice(ref_std_lat, ref->std_lattice);
  CHECK((got->std_lattice - ref_std_lat).cwiseAbs().maxCoeff() < 1e-5);

  // Equivalent atoms (per input atom).
  REQUIRE(static_cast<int>(got->equivalent_atoms.size()) == ref->n_atoms);
  for (int i = 0; i < ref->n_atoms; ++i)
    CHECK(got->equivalent_atoms[static_cast<std::size_t>(i)] ==
          ref->equivalent_atoms[i]);

  // Standardized cell atoms, matched as a set (type + position + tensor).
  int const n_std = ref->n_std_atoms;
  bool const collinear = mcell.rank() == SiteTensor::collinear;
  for (int r = 0; r < n_std; ++r) {
    Vector3d const ref_pos(ref->std_positions[r][0], ref->std_positions[r][1],
                           ref->std_positions[r][2]);
    bool found = false;
    for (int g = 0; g < n_std && !found; ++g) {
      auto const ug = static_cast<std::size_t>(g);
      if (got->std_types[ug] != ref->std_types[r])
        continue;
      if (!fractional_overlap(got->std_positions.row(g).transpose(), ref_pos,
                              symprec))
        continue;
      bool tensor_ok = true;
      if (collinear) {
        tensor_ok = std::abs(std::get<CollinearTensors>(got->std_tensors)[ug] -
                             ref->std_tensors[r]) < 1e-4;
      } else {
        auto const v = std::get<NoncollinearTensors>(got->std_tensors).row(ug);
        for (int s = 0; s < 3; ++s)
          tensor_ok = tensor_ok &&
                      std::abs(v[s] - ref->std_tensors[3 * r + s]) < 1e-4;
      }
      found = tensor_ok;
    }
    INFO("standardized atom " << r << " not matched");
    CHECK(found);
  }

  // Magnetic operations, matched as a set.
  REQUIRE(static_cast<int>(got->operations.size()) == ref->n_operations);
  for (int s = 0; s < ref->n_operations; ++s) {
    Matrix3i r;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        r(i, j) = ref->rotations[s][i][j];
    Vector3d const t(ref->translations[s][0], ref->translations[s][1],
                     ref->translations[s][2]);
    bool const tr = ref->time_reversals[s] != 0;
    bool found = false;
    for (auto const &op : got->operations) {
      if (op.rotation == r && op.time_reversal == tr &&
          fractional_overlap(op.translation, t, symprec)) {
        found = true;
        break;
      }
    }
    INFO("operation " << s << " not matched");
    CHECK(found);
  }

  Matrix3d ref_prim;
  cppcrystal::oracle::from_c_lattice(ref_prim, ref->primitive_lattice);
  CHECK(same_lattice(got->primitive_lattice, ref_prim));

  spg_free_magnetic_dataset(ref);
}

} // namespace

TEST_CASE("magnetic dataset: collinear ferromagnet", "[oracle][magnetic]") {
  Cell const cell = make_cell(4.0, {{0.0, 0.0, 0.0}}, {0});
  MagneticCell const mcell(cell, SiteTensors{CollinearTensors{1.0}});
  check(mcell, TensorKind::axial, 1e-5);
}

TEST_CASE("magnetic dataset: collinear antiferromagnet", "[oracle][magnetic]") {
  Cell const cell =
      make_cell(4.0, {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, {0, 0});
  MagneticCell const mcell(cell, SiteTensors{CollinearTensors{1.0, -1.0}});
  check(mcell, TensorKind::axial, 1e-5);
}

TEST_CASE("magnetic dataset: non-collinear antiferromagnet",
          "[oracle][magnetic]") {
  Cell const cell =
      make_cell(4.0, {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, {0, 0});
  MagneticCell const mcell(
      cell, SiteTensors{noncollinear_tensors(
                {Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, -1.0)})});
  check(mcell, TensorKind::axial, 1e-5);
}

TEST_CASE("magnetic dataset: type-II grey group", "[oracle][magnetic]") {
  Cell const cell =
      make_cell(4.0, {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, {0, 0});
  MagneticCell const mcell(cell, SiteTensors{CollinearTensors{1.0, 1.0}});
  check(mcell, TensorKind::axial, 1e-5);
}
