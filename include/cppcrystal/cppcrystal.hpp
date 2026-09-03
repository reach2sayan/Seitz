#pragma once

// Umbrella header: the full public API. The determination pipeline
// (symmetry/, spacegroup/, refine/, spin/, magnetic/, and the numeric and
// container helpers behind them) is private to src/ and deliberately absent.

#include <cppcrystal/analysis/analyzer.hpp>
#include <cppcrystal/analysis/dataset.hpp>
#include <cppcrystal/analysis/magnetic_symmetry_analyzer.hpp>
#include <cppcrystal/analysis/symmetry_analyzer.hpp>
#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/fractional.hpp>
#include <cppcrystal/core/lattice.hpp>
#include <cppcrystal/core/magnetic_cell.hpp>
#include <cppcrystal/core/magnetic_symmetry_operation.hpp>
#include <cppcrystal/core/periodicity.hpp>
#include <cppcrystal/core/point_group.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/core/types.hpp>
#include <cppcrystal/core/version.hpp>
#include <cppcrystal/data/msg_database.hpp>
#include <cppcrystal/data/spg_database.hpp>
#include <cppcrystal/generate/crystal_builder.hpp>
#include <cppcrystal/generate/random_lattice.hpp>
#include <cppcrystal/generate/wyckoff_combinations.hpp>
#include <cppcrystal/group/space_group.hpp>
#include <cppcrystal/group/wyckoff.hpp>
#include <cppcrystal/kpoint/brillouin_zone.hpp>
#include <cppcrystal/kpoint/grid.hpp>
#include <cppcrystal/kpoint/reciprocal_mesh.hpp>
#include <cppcrystal/kpoint/reciprocal_mesh_builder.hpp>
#include <cppcrystal/spacegroup_match.hpp>
#include <cppcrystal/warmup.hpp>
