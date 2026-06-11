#include <spglib/group/wyckoff_position.hpp>

#include <spglib/math/fractional.hpp>

#include <algorithm>
#include <vector>

namespace spglib::group {

namespace {
// Tolerance for comparing exact (rational) database coordinates mod 1.
constexpr double kGroupTol = 1e-5;

[[nodiscard]] bool same_point(Vector3d const &a, Vector3d const &b) {
  return math::frac_displacement(a, b).cwiseAbs().maxCoeff() < kGroupTol;
}
} // namespace

char WyckoffPosition::letter() const noexcept {
  return letter_ < 26 ? static_cast<char>('a' + letter_)
                      : static_cast<char>('A' + (letter_ - 26));
}

Positions WyckoffPosition::get_all_positions(Vector3d const &xyz) const {
  // Project the input onto the position's canonical coordinate so a point only
  // approximately on the position still yields the exact orbit.
  Vector3d const canonical = repr_rotation_.cast<double>() * xyz +
                             repr_translation_;

  std::vector<Vector3d> orbit;
  orbit.reserve(orbit_ops_.size());
  for (auto const &op : orbit_ops_) {
    Vector3d const image = math::wrap_to_unit_cell(op.apply(canonical));
    bool const seen = std::ranges::any_of(
        orbit, [&](Vector3d const &p) { return same_point(p, image); });
    if (!seen) {
      orbit.push_back(image);
    }
  }

  Positions out(static_cast<Index>(orbit.size()), 3);
  for (Index i = 0; i < static_cast<Index>(orbit.size()); ++i) {
    out.row(i) = orbit[static_cast<std::size_t>(i)].transpose();
  }
  return out;
}

} // namespace spglib::group
