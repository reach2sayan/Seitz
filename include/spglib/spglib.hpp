#pragma once

// Umbrella header for the CppCrystal port of spglib.
// Phase 0 exposes only the foundations; module headers are added here as the
// port progresses (cell, symmetry, spacegroup, magnetic, kpoints, ...).

#include <spglib/core/cell.hpp>
#include <spglib/core/error.hpp>
#include <spglib/core/overlap.hpp>
#include <spglib/core/point_group.hpp>
#include <spglib/core/symmetry_operation.hpp>
#include <spglib/core/tolerance.hpp>
#include <spglib/core/types.hpp>
#include <spglib/core/version.hpp>
#include <spglib/data/spg_database.hpp>
#include <spglib/math/fractional.hpp>
#include <spglib/math/integer_matrix.hpp>
#include <spglib/reduce/delaunay.hpp>
#include <spglib/reduce/niggli.hpp>
#include <spglib/spacegroup/hall_symbol.hpp>
#include <spglib/symmetry/find_symmetry.hpp>
#include <spglib/symmetry/pointgroup.hpp>
#include <spglib/symmetry/primitive.hpp>
