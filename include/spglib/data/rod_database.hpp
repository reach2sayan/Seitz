#pragma once

#include <spglib/core/symmetry_operation.hpp>

#include <string_view>

namespace spglib::data {

// Decoder for the rod-group (1D-periodic subperiodic group) operation tables in
// the generated header data/rod_group_tables.hpp (produced by
// tools/transcribe_rod_groups.py from PyXtal's Group(n, dim=1)).
//
// Rod groups are not in spglib, so unlike the layer groups (reached through
// spglib's negative-Hall datasets) they have no SpacegroupType metadata and no
// site-symmetry database — only their symmetry operations are tabulated, and the
// rod Wyckoff positions are derived in-house from those operations (group::
// RodGroup), the same way group::PointGroup derives the 0D point-group Wyckoffs.

// Number of rod groups (75).
[[nodiscard]] int num_rod_groups() noexcept;

// The single periodic axis (0=a, 1=b, 2=c). PyXtal's convention is c (= 2): the
// orbit expansion folds only this component; the other two are aperiodic.
[[nodiscard]] int rod_periodic_axis() noexcept;

// Hermann-Mauguin symbol of rod group `rod_number` (1..75); empty out of range.
[[nodiscard]] std::string_view rod_symbol(int rod_number) noexcept;

// The symmetry operations of rod group `rod_number` (1..75): one period's worth
// of factor-group representatives (the pure unit translation along the periodic
// axis is implicit, as in the space-group operation tables). Empty out of range.
// Port-style analogue of data::operations_from_database for the rod tables.
[[nodiscard]] SymmetryOperations rod_operations_from_database(int rod_number);

} // namespace spglib::data
