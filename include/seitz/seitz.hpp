#pragma once

// Umbrella header: the whole public API and nothing else. The determination
// pipeline (symmetry/, spacegroup/, refine/, spin/, magnetic/), the generated
// tables and the helpers behind them live under src/ and are unreachable from
// here, so a TU compiling against this header compiles against the installed
// library.

// Vocabulary: the values everything else is phrased in.
#include <seitz/core/cell.hpp>
#include <seitz/core/error.hpp>
#include <seitz/core/fractional.hpp>
#include <seitz/core/keys.hpp>
#include <seitz/core/lattice.hpp>
#include <seitz/core/magnetic_cell.hpp>
#include <seitz/core/magnetic_symmetry_operation.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/periodicity.hpp>
#include <seitz/core/point_group.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/core/tolerance.hpp>
#include <seitz/core/types.hpp>
#include <seitz/core/version.hpp>

// Determination: a structure in, its symmetry out.
#include <seitz/analysis/analyzer.hpp>
#include <seitz/analysis/dataset.hpp>
#include <seitz/analysis/magnetic_symmetry_analyzer.hpp>
#include <seitz/analysis/symmetry_analyzer.hpp>
#include <seitz/spacegroup_match.hpp>

// The catalogs: the databases as queryable objects, structure-free.
#include <seitz/data/element_data.hpp>
#include <seitz/data/msg_database.hpp>
#include <seitz/data/spg_database.hpp>
#include <seitz/group/point_group.hpp>
#include <seitz/group/rod_group.hpp>
#include <seitz/group/space_group.hpp>
#include <seitz/group/subgroup_graph.hpp>
#include <seitz/group/wyckoff.hpp>

// Construction: random structures with a prescribed symmetry.
#include <seitz/generate/assignments.hpp>
#include <seitz/generate/distance_check.hpp>
#include <seitz/generate/generator.hpp>

// Reciprocal space.
// Alloy configurational thermodynamics: symmetry-distinct cluster enumeration
// and the Cluster Variation Method.
#include <seitz/alloy/cluster.hpp>
#include <seitz/alloy/clusters_pool.hpp>
#include <seitz/alloy/cvm.hpp>
#include <seitz/alloy/parent_lattice.hpp>
#include <seitz/alloy/site_basis.hpp>

#include <seitz/kpoint/mesh.hpp>

#include <seitz/warmup.hpp>
