#pragma once

#include <seitz/core/cell.hpp>
#include <seitz/core/error.hpp>
#include <seitz/core/periodicity.hpp>
#include <seitz/core/types.hpp>
#include <seitz/generate/assignments.hpp>
#include <seitz/group/point_group.hpp>
#include <seitz/group/rod_group.hpp>
#include <seitz/group/space_group.hpp>

#include <concepts>
#include <cstddef>
#include <random>
#include <ranges>
#include <string_view>
#include <vector>

#pragma GCC visibility push(default)

namespace seitz::generate {

// One generated structure:
struct Generated {
  Cell cell;
  Assignment<group::Wyckoff> assignment;
};

// The fractional box a seed coordinate is drawn from, per axis. The seed is
// projected onto the Wyckoff locus first, so this only fixes where in the cell
// the structure lives: the full repeat along a periodic axis, a band centred
// on 0 along an aperiodic one, where an axis-flipping op must map the
// structure onto itself (image at -z, not 1 - z).
struct SeedBox {
  Vector3d low{Vector3d::Zero()};
  Vector3d high{Vector3d::Ones()};
  [[nodiscard]] Vector3d sample(std::mt19937_64 &rng) const;
};

// How one group family realizes its geometry: its name in an error, its cell
// periodicity, the seed box, and the random metric of one attempt. Everything
// else -- assignment enumeration, the general-position restriction, the
// shuffled assignment/attempt search, distance acceptance -- is
// family-independent and lives in Generator.
//
// A layer group is a SpaceGroup whose Hall key names the layer family; the
// traits read that off the group, not a second entry point.
template <class G> struct GroupTraits;

template <> struct GroupTraits<group::SpaceGroup> {
  [[nodiscard]] static std::string_view
  kind(group::SpaceGroup const &g) noexcept;
  [[nodiscard]] static CellPeriodicity
  periodicity(group::SpaceGroup const &g) noexcept;
  [[nodiscard]] static SeedBox seed_box(group::SpaceGroup const &g) noexcept;
  [[nodiscard]] static Matrix3d lattice(group::SpaceGroup const &g,
                                        Composition const &comp,
                                        GenerateOptions const &options,
                                        int attempt, std::mt19937_64 &rng);
};

template <> struct GroupTraits<group::PointGroup> {
  [[nodiscard]] static std::string_view
  kind(group::PointGroup const &) noexcept;
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

// A group the generator can place a composition on: it has Wyckoff positions
// and a GroupTraits<G> specialization.
template <class G>
concept GeneratableGroup = requires(G const &g, Composition const &comp,
                                    GenerateOptions const &options,
                                    std::mt19937_64 &rng) {
  { g.wyckoffs() } -> std::ranges::input_range;
  { GroupTraits<G>::kind(g) } -> std::convertible_to<std::string_view>;
  { GroupTraits<G>::periodicity(g) } -> std::same_as<CellPeriodicity>;
  { GroupTraits<G>::seed_box(g) } -> std::same_as<SeedBox>;
  { GroupTraits<G>::lattice(g, comp, options, 0, rng) } -> std::same_as<Matrix3d>;
};

// Assignments per search: ample for any real composition, a bound on
// pathological ones.
inline constexpr std::size_t kMaxAssignments = 1000;

// The random-structure generator
template <GeneratableGroup G> class Generator {
public:
  explicit Generator(G const &group, GenerateOptions options = {}) noexcept
      : group_{&group}, options_{options} {}
  // Whether `comp` can be placed on this group at all, the option's fixed
  // sites included.
  [[nodiscard]] bool compatible(Composition const &comp) const;

  // The first `max` valid Wyckoff assignments of `comp`, each beginning with
  // the fixed sites; fewer when there are fewer, exactly `max` when the
  // enumeration was cut short.
  [[nodiscard]] std::vector<Assignment<group::Wyckoff>>
  assignments(Composition const &comp,
              std::size_t max = kMaxAssignments) const;

  // A random structure with this group's symmetry and `comp`'s composition,
  // deterministic in options.seed. The returned assignment carries every
  // orbit's generating coordinate.
  // Errors : (1) the composition has no Wyckoff assignment here (the fixed
  //              sites counted),
  //          (2) a fixed site names no position of this group,
  //          (3) Placement::general_only cannot be met,
  //          (4) the caller's lattice is not invariant under the group,
  //          (5) or no distance-valid structure turns up within the attempt
  //              budget.
  [[nodiscard]] Result<Generated> operator()(Composition const &comp) const;

private:
  // options_.sites resolved to placements on this group; errors on a letter
  // the group does not have.
  [[nodiscard]] Result<Assignment<group::Wyckoff>> fixed_placements() const;

  G const *group_;
  GenerateOptions options_;
};

extern template class Generator<group::SpaceGroup>;
extern template class Generator<group::PointGroup>;
extern template class Generator<group::RodGroup>;

} // namespace seitz::generate

#pragma GCC visibility pop
