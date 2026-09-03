#pragma once

// Umbrella header: the whole public API, and nothing else. The determination
// pipeline (symmetry/, spacegroup/, refine/, spin/, magnetic/), the generated
// operation tables and the numeric and container helpers behind them all live
// under src/ and are unreachable from here by construction — a translation unit
// that compiles against this header compiles against the installed library.

// Vocabulary: the values everything else is phrased in.
#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/fractional.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/lattice.hpp>
#include <cppcrystal/core/magnetic_cell.hpp>
#include <cppcrystal/core/magnetic_symmetry_operation.hpp>
#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/periodicity.hpp>
#include <cppcrystal/core/point_group.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/core/types.hpp>
#include <cppcrystal/core/version.hpp>

// Determination: a structure in, its symmetry out.
#include <cppcrystal/analysis/analyzer.hpp>
#include <cppcrystal/analysis/dataset.hpp>
#include <cppcrystal/analysis/magnetic_symmetry_analyzer.hpp>
#include <cppcrystal/analysis/symmetry_analyzer.hpp>
#include <cppcrystal/spacegroup_match.hpp>

// The catalogs: the databases as queryable objects, structure-free.
#include <cppcrystal/data/element_data.hpp>
#include <cppcrystal/data/msg_database.hpp>
#include <cppcrystal/data/sitesym_database.hpp>
#include <cppcrystal/data/spg_database.hpp>
#include <cppcrystal/group/point_group.hpp>
#include <cppcrystal/group/rod_group.hpp>
#include <cppcrystal/group/space_group.hpp>
#include <cppcrystal/group/subgroup_graph.hpp>
#include <cppcrystal/group/wyckoff.hpp>

// Construction: random structures with a prescribed symmetry.
#include <cppcrystal/generate/assignments.hpp>
#include <cppcrystal/generate/distance_check.hpp>
#include <cppcrystal/generate/generator.hpp>

// Reciprocal space.
#include <cppcrystal/kpoint/mesh.hpp>

#include <cppcrystal/warmup.hpp>
