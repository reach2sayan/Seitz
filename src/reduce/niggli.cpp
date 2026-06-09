#include <spglib/reduce/niggli.hpp>

#include <array>
#include <cmath>

// Krivy-Gruber Niggli reduction, a faithful port of niggli.c (space-group path,
// aperiodic_axis = -1). The cell is represented by the six metric parameters
//   A = a.a, B = b.b, C = c.c, xi = 2 b.c, eta = 2 a.c, zeta = 2 a.b
// and each step right-multiplies the lattice by an integer matrix `tmat`,
// transforming the basis columns. The control flow (which steps restart the
// pass) mirrors the reference exactly.
namespace spglib::reduce {

namespace {

constexpr int kMaxAttempts = 1000; // spglib default (SPGLIB_NUM_ATTEMPTS)

struct State {
  Matrix3d lattice;
  Matrix3d tmat = Matrix3d::Identity();
  double A = 0, B = 0, C = 0, xi = 0, eta = 0, zeta = 0, eps = 0;
  int l = 0, m = 0, n = 0;
};

void set_angle_types(State &p) {
  p.l = p.xi < -p.eps ? -1 : (p.xi > p.eps ? 1 : 0);
  p.m = p.eta < -p.eps ? -1 : (p.eta > p.eps ? 1 : 0);
  p.n = p.zeta < -p.eps ? -1 : (p.zeta > p.eps ? 1 : 0);
}

void set_parameters(State &p) {
  Matrix3d const g = p.lattice.transpose() * p.lattice;
  p.A = g(0, 0);
  p.B = g(1, 1);
  p.C = g(2, 2);
  p.xi = 2 * g(1, 2);
  p.eta = 2 * g(0, 2);
  p.zeta = 2 * g(0, 1);
  set_angle_types(p);
}

void reset(State &p) {
  p.lattice = p.lattice * p.tmat;
  set_parameters(p);
}

bool step1(State &p) {
  if (p.A > p.B + p.eps || (!(std::fabs(p.A - p.B) > p.eps) &&
                            std::fabs(p.xi) > std::fabs(p.eta) + p.eps)) {
    p.tmat << 0, -1, 0, -1, 0, 0, 0, 0, -1;
    return true;
  }
  return false;
}

bool step2(State &p) {
  if (p.B > p.C + p.eps || (!(std::fabs(p.B - p.C) > p.eps) &&
                            std::fabs(p.eta) > std::fabs(p.zeta) + p.eps)) {
    p.tmat << -1, 0, 0, 0, 0, -1, 0, -1, 0;
    return true;
  }
  return false;
}

bool step3(State &p) {
  if (p.l * p.m * p.n == 1) {
    int const i = p.l == -1 ? -1 : 1;
    int const j = p.m == -1 ? -1 : 1;
    int const k = p.n == -1 ? -1 : 1;
    p.tmat << i, 0, 0, 0, j, 0, 0, 0, k;
    return true;
  }
  return false;
}

bool step4(State &p) {
  if (p.l == -1 && p.m == -1 && p.n == -1) {
    return false;
  }
  if (p.l * p.m * p.n == 0 || p.l * p.m * p.n == -1) {
    int i = 1, j = 1, k = 1, r = -1; // r: which of i/j/k is the "free" axis
    // clang-format off
    if (p.l == 1) { i = -1; }
    if (p.l == 0) { r = 0;  }
    if (p.m == 1) { j = -1; }
    if (p.m == 0) { r = 1;  }
    if (p.n == 1) { k = -1; }
    if (p.n == 0) { r = 2;  }
    if (i * j * k == -1) {
      if (r == 0) { i = -1; }
      if (r == 1) { j = -1; }
      if (r == 2) { k = -1; }
    }
    // clang-format on
    p.tmat << i, 0, 0, 0, j, 0, 0, 0, k;
    return true;
  }
  return false;
}

bool step5(State &p) {
  if (std::fabs(p.xi) > p.B + p.eps ||
      (!(std::fabs(p.B - p.xi) > p.eps) && 2 * p.eta < p.zeta - p.eps) ||
      (!(std::fabs(p.B + p.xi) > p.eps) && p.zeta < -p.eps)) {
    p.tmat << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    if (p.xi > 0) {
      p.tmat(1, 2) = -1;
    }
    if (p.xi < 0) {
      p.tmat(1, 2) = 1;
    }
    return true;
  }
  return false;
}

bool step6(State &p) {
  if (std::fabs(p.eta) > p.A + p.eps ||
      (!(std::fabs(p.A - p.eta) > p.eps) && 2 * p.xi < p.zeta - p.eps) ||
      (!(std::fabs(p.A + p.eta) > p.eps) && p.zeta < -p.eps)) {
    p.tmat << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    if (p.eta > 0) {
      p.tmat(0, 2) = -1;
    }
    if (p.eta < 0) {
      p.tmat(0, 2) = 1;
    }
    return true;
  }
  return false;
}

bool step7(State &p) {
  if (std::fabs(p.zeta) > p.A + p.eps ||
      (!(std::fabs(p.A - p.zeta) > p.eps) && 2 * p.xi < p.eta - p.eps) ||
      (!(std::fabs(p.A + p.zeta) > p.eps) && p.eta < -p.eps)) {
    p.tmat << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    if (p.zeta > 0) {
      p.tmat(0, 1) = -1;
    }
    if (p.zeta < 0) {
      p.tmat(0, 1) = 1;
    }
    return true;
  }
  return false;
}

bool step8(State &p) {
  if (p.xi + p.eta + p.zeta + p.A + p.B < -p.eps ||
      (!(std::fabs(p.xi + p.eta + p.zeta + p.A + p.B) > p.eps) &&
       2 * (p.A + p.eta) + p.zeta > p.eps)) {
    p.tmat << 1, 0, 1, 0, 1, 1, 0, 0, 1;
    return true;
  }
  return false;
}

} // namespace

Result<Matrix3d> niggli_reduce(Matrix3d const &lattice, double eps) {
  using StepFn = bool (*)(State &);
  std::array<StepFn, 8> const steps{step1, step2, step3, step4,
                                    step5, step6, step7, step8};

  State p;
  p.lattice = lattice;
  p.eps = eps;
  set_parameters(p);

  bool converged = false;
  for (int attempt = 0; attempt < kMaxAttempts && !converged; ++attempt) {
    int j = 0;
    for (; j < 8; ++j) {
      if (steps[static_cast<std::size_t>(j)](p)) {
        reset(p);
        // Steps 2, 5, 6, 7, 8 restart the pass from step 1.
        if (j == 1 || j == 4 || j == 5 || j == 6 || j == 7)
          break;
      }
    }
    if (j == 8) {
      converged = true;
    }
  }

  if (!converged) {
    return leaf::new_error(e_niggli_failed{});
  }
  return p.lattice;
}

} // namespace spglib::reduce
