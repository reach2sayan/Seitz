#pragma once

#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/periodicity.hpp>
#include <cppcrystal/core/types.hpp>
#include <cppcrystal/generate/assignments.hpp>
#include <cppcrystal/group/point_group.hpp>
#include <cppcrystal/group/rod_group.hpp>
#include <cppcrystal/group/space_group.hpp>

#include <cstddef>
#include <random>
#include <ranges>
#include <string_view>
#include <vector>

namespace cppcrystal::generate {

// One generated structure: a cell that describes itself, plus the Wyckoff
// assignment it realizes. The cell's periodicity says which family produced it
// — all-periodic for a 3D crystal, one aperiodic axis for a layer, one periodic
// axis for a rod, none for a cluster — so a cluster is a Cell like any other:
// its "lattice" is the metric its point group is isometric in and its positions
// are fractional in that metric. There is no second, Cartesian shape.
struct Generated {
  Cell cell;
  Assignment<group::Wyckoff> assignment;
};

// The fractional box a seed coordinate is drawn from, per axis. The seed is
// projected onto the Wyckoff locus before its orbit is built, so this only has
// to say where in the cell the structure should live: the full repeat along a
// periodic axis, a band centred on 0 along an aperiodic one, where an
// axis-flipping operation must map the structure onto itself (image at -z, not
// at 1 - z) rather than away from it.
struct SeedBox {
  Vector3d low{Vector3d::Zero()};
  Vector3d high{Vector3d::Ones()};

  [[nodiscard]] Vector3d sample(std::mt19937_64 &rng) const;
};

// How one group family realizes its geometry: what to call it in an error, how
// its cell is periodic, where in that cell a seed coordinate belongs, and the
// random metric one attempt is made in. Everything else about generation — the
// assignment enumeration, the general-position restriction, the shuffled
// assignment/attempt search, the distance acceptance — is family-independent
// and lives in Generator.
//
// A layer group is a SpaceGroup whose Hall key names the layer family; the
// traits read that from the group rather than taking a second entry point.
template <class G> struct GroupTraits;

template <> struct GroupTraits<group::SpaceGroup> {
  [[nodiscard]] static std::string_view kind(group::SpaceGroup const &g) noexcept;
  [[nodiscard]] static CellPeriodicity
  periodicity(group::SpaceGroup const &g) noexcept;
  [[nodiscard]] static SeedBox seed_box(group::SpaceGroup const &g) noexcept;
  [[nodiscard]] static Matrix3d lattice(group::SpaceGroup const &g,
                                        Composition const &comp,
                                        GenerateOptions const &options,
                                        int attempt, std::mt19937_64 &rng);
};

template <> struct GroupTraits<group::PointGroup> {
  [[nodiscard]] static std::string_view kind(group::PointGroup const &) noexcept;
  [[nodiscard]] static CellPeriodicity
  periodicity(group::PointGroup const &) noexcept;
  [[nodiscard]] static SeedBox seed_box(group::PointGroup const &) noexcept;
  [[nodiscard]] static Matrix3d lattice(group::PointGroup const &g,
                                        Composition const &comp,
                                        GenerateOptions const &options,
                                        int attempt, std::mt19937_64 &rng);
};

template <> struct GroupTraits<group::RodGroup> {
  [[nodiscard]] static std::string_view kind(group::RodGroup const &) noexcept;
  [[nodiscard]] static CellPeriodicity
  periodicity(group::RodGroup const &g) noexcept;
  [[nodiscard]] static SeedBox seed_box(group::RodGroup const &g) noexcept;
  [[nodiscard]] static Matrix3d lattice(group::RodGroup const &g,
                                        Composition const &comp,
                                        GenerateOptions const &options,
                                        int attempt, std::mt19937_64 &rng);
};

// Assignments considered per search: enough for any real composition, a bound
// for pathological ones.
inline constexpr std::size_t kMaxAssignments = 1000;

// The random-structure generator, over any group family that offers Wyckoff
// positions. `G` picks the geometry through GroupTraits<G>; the search itself
// is the same for a 3D crystal, a layer, a rod and a cluster.
//
// Non-owning: the group must outlive the generator (SpaceGroup::of hands back a
// shared flyweight, the others are owned by the caller).
template <class G> class Generator {
public:
  explicit Generator(G const &group, GenerateOptions options = {}) noexcept
      : group_(&group), options_(options) {}

  // Whether `comp` can be placed on this group at all.
  [[nodiscard]] bool compatible(Composition const &comp) const {
    return assignable(group_->wyckoffs(), comp);
  }

  // The first `max` valid Wyckoff assignments of `comp` (see
  // enumerate_assignments for the rules); fewer when there are fewer, exactly
  // `max` when the enumeration was cut short.
  [[nodiscard]] std::vector<Assignment<group::Wyckoff>>
  assignments(Composition const &comp,
              std::size_t max = kMaxAssignments) const {
    return std::ranges::to<std::vector<Assignment<group::Wyckoff>>>(
        enumerate_assignments(group_->wyckoffs(), comp) |
        std::views::take(max));
  }

  // A random structure with this group's symmetry and `comp`'s composition,
  // deterministic in options.seed. Errors when the composition has no Wyckoff
  // assignment here, when Placement::general_only cannot be met, or when no
  // distance-valid structure turns up within the attempt budget.
  [[nodiscard]] Result<Generated> operator()(Composition const &comp) const;

private:
  G const *group_;
  GenerateOptions options_;
};

extern template class Generator<group::SpaceGroup>;
extern template class Generator<group::PointGroup>;
extern template class Generator<group::RodGroup>;

} // namespace cppcrystal::generate
