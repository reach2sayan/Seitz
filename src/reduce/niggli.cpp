#include <cppcrystal/core/lattice.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>

// Krivy-Gruber Niggli reduction (space-group path). The cell is represented by
// the six metric parameters
//   A = a.a, B = b.b, C = c.c, xi = 2 b.c, eta = 2 a.c, zeta = 2 a.b
// and each step right-multiplies the lattice by an integer matrix `tmat`,
// transforming the basis columns. Which steps restart the pass is part of the
// algorithm (see the restart set below) and must not be changed.
namespace cppcrystal {

namespace {

constexpr int kMaxAttempts = 1000;

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

// Step 4 sign choice: flip every axis whose angle type is positive; if that
// leaves an odd number of flips, also flip the last axis whose type is zero.
[[nodiscard]] std::array<int, 3> step4_signs(std::array<int, 3> const &types) {
  std::array<int, 3> signs{};
  std::ranges::transform(types, signs.begin(),
                         [](int t) { return t == 1 ? -1 : 1; });
  if (signs[0] * signs[1] * signs[2] == -1) {
    if (auto const free = std::ranges::find_last(types, 0); !free.empty()) {
      signs[static_cast<std::size_t>(free.begin() - types.begin())] = -1;
    }
  }
  return signs;
}

bool step4(State &p) {
  if (p.l == -1 && p.m == -1 && p.n == -1) {
    return false;
  }
  if (p.l * p.m * p.n == 0 || p.l * p.m * p.n == -1) {
    auto const [i, j, k] = step4_signs({p.l, p.m, p.n});
    p.tmat << i, 0, 0, 0, j, 0, 0, 0, k;
    return true;
  }
  return false;
}

// Steps 5-7 share one shape: a doubled off-diagonal metric term too large
// against a diagonal one (or equal to it, with an ordering tie-break on the
// other two terms) shears the lattice by one unit in cell (row, col) against
// the sign of the term.
using Term = double State::*;
bool shear_step(State &p, Term term, Term diag, Term first, Term second,
                Index row, Index col) {
  double const t = p.*term;
  double const d = p.*diag;
  if (std::fabs(t) > d + p.eps ||
      (!(std::fabs(d - t) > p.eps) && 2 * p.*first < p.*second - p.eps) ||
      (!(std::fabs(d + t) > p.eps) && p.*second < -p.eps)) {
    p.tmat = Matrix3d::Identity();
    p.tmat(row, col) = t > 0 ? -1.0 : (t < 0 ? 1.0 : 0.0);
    return true;
  }
  return false;
}

bool step5(State &p) {
  return shear_step(p, &State::xi, &State::B, &State::eta, &State::zeta, 1, 2);
}

bool step6(State &p) {
  return shear_step(p, &State::eta, &State::A, &State::xi, &State::zeta, 0, 2);
}

bool step7(State &p) {
  return shear_step(p, &State::zeta, &State::A, &State::xi, &State::eta, 0, 1);
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

Result<Lattice> Lattice::niggli(double eps) const {
  Matrix3d const &lattice = basis_;
  struct Step {
    bool (*apply)(State &);
    bool restarts; // firing this step restarts the pass from step 1
  };
  // Restart set {step2, step5, step6, step7, step8}: those steps restart the
  // pass, while step1/3/4 fall through to the next step.
  constexpr std::array<Step, 8> steps{{
      {step1, false},
      {step2, true},
      {step3, false},
      {step4, false},
      {step5, true},
      {step6, true},
      {step7, true},
      {step8, true},
  }};

  State p;
  p.lattice = lattice;
  p.eps = eps;
  set_parameters(p);

  bool converged = false;
  for (int attempt = 0; attempt < kMaxAttempts && !converged; ++attempt) {
    bool restarted = false;
    for (auto const &[apply, restarts] : steps) {
      if (std::invoke(apply, p)) {
        reset(p);
        if (restarts) {
          restarted = true;
          break;
        }
      }
    }
    // A full pass with no restarting step means the cell is Niggli-reduced.
    converged = !restarted;
  }

  if (!converged) {
    return leaf::new_error(e_niggli_failed{});
  }
  return Lattice{p.lattice};
}

} // namespace cppcrystal
