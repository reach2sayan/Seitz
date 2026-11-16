#include <seitz/analysis/symmetry_analyzer.hpp>
#include <seitz/data/spg_database.hpp>

#include <cstdio>

// Rock salt in its conventional cell: Fm-3m, number 225.
int main() {
  using namespace seitz;
  Positions positions(8, 3);
  positions << 0, 0, 0, 0, .5, .5, .5, 0, .5, .5, .5, 0, //
      .5, .5, .5, .5, 0, 0, 0, .5, 0, 0, 0, .5;
  Cell cell{Lattice{Matrix3d{5.64 * Matrix3d::Identity()}}, positions,
            Types{11, 11, 11, 11, 17, 17, 17, 17}};
  auto const sa = analysis::SymmetryAnalyzer::from_cell(std::move(cell));
  auto const hall = sa.hall();
  if (!hall) {
    return 1;
  }
  int const number = data::spacegroup_type(*hall).number;
  std::printf("rock salt: space group %d\n", number);
  return number == 225 ? 0 : 1;
}
