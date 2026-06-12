#pragma once

// Umbrella header for the CppCrystal port of spglib.
// Phase 0 exposes only the foundations; module headers are added here as the
// port progresses (cell, symmetry, spacegroup, magnetic, kpoints, ...).

#include <spglib/analysis/magnetic_symmetry_analyzer.hpp>
#include <spglib/analysis/operation_set.hpp>
#include <spglib/analysis/symmetry_analyzer.hpp>
#include <spglib/core/cell.hpp>
#include <spglib/core/centering.hpp>
#include <spglib/core/error.hpp>
#include <spglib/core/magnetic_cell.hpp>
#include <spglib/core/magnetic_symmetry_operation.hpp>
#include <spglib/core/overlap.hpp>
#include <spglib/core/point_group.hpp>
#include <spglib/core/symmetry_operation.hpp>
#include <spglib/core/tolerance.hpp>
#include <spglib/core/types.hpp>
#include <spglib/core/version.hpp>
#include <spglib/data/msg_database.hpp>
#include <spglib/data/spg_database.hpp>
#include <spglib/dataset.hpp>
#include <spglib/generate/crystal_builder.hpp>
#include <spglib/generate/random_lattice.hpp>
#include <spglib/generate/wyckoff_combinations.hpp>
#include <spglib/group/space_group.hpp>
#include <spglib/group/wyckoff_position.hpp>
#include <spglib/standardize.hpp>
#include <spglib/kpoint/brillouin_zone.hpp>
#include <spglib/kpoint/grid.hpp>
#include <spglib/kpoint/reciprocal_mesh.hpp>
#include <spglib/kpoint/reciprocal_mesh_builder.hpp>
#include <spglib/magnetic/magnetic_spacegroup.hpp>
#include <spglib/magnetic_dataset.hpp>
#include <spglib/math/fractional.hpp>
#include <spglib/math/integer_matrix.hpp>
#include <spglib/reduce/delaunay.hpp>
#include <spglib/reduce/niggli.hpp>
#include <spglib/spacegroup/hall_symbol.hpp>
#include <spglib/spacegroup/spacegroup.hpp>
#include <spglib/warmup.hpp>
#include <spglib/spin/spin.hpp>
#include <spglib/symmetry/find_symmetry.hpp>
#include <spglib/symmetry/pointgroup.hpp>
#include <spglib/symmetry/primitive.hpp>
