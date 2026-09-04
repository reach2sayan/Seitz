// The public API in one file: everything here compiles against <cppcrystal/…>
// alone, with none of src/ on the include path, so this program is also the
// guard that the installed interface is self-contained.

#include <cppcrystal/cppcrystal.hpp>

#include <print>

using namespace cppcrystal;

int main() {
  std::println("CppCrystal {} — modern C++23 crystallography library "
               "(validated against reference spglib v{}.{}.{})",
               version_string(), kReferenceSpglibVersion.major,
               kReferenceSpglibVersion.minor, kReferenceSpglibVersion.patch);

  // 1. Determination: a structure in, its symmetry out. The analyzer memoizes,
  //    so it is the object you keep, not a call you repeat.
  Positions positions(2, 3);
  positions << 0.0, 0.0, 0.0, 0.5, 0.5, 0.5;
  Cell const bcc{Lattice{Matrix3d::Identity() * 3.0}, positions, Types{26, 26}};

  auto const analyzer = analysis::SymmetryAnalyzer::from_cell(bcc);
  return leaf::try_handle_all(
      [&]() -> Result<int> {
        BOOST_LEAF_AUTO(hall, analyzer.hall());
        BOOST_LEAF_AUTO(ops, analyzer.operations());
        auto const &type = data::spacegroup_type(hall);
        // international_short is a string_view; %s would have read past its
        // end for any table entry that was not a literal.
        std::println("  bcc Fe -> {} (No. {}), {} operations",
                     type.international_short, type.number, ops.size());

        // 2. The catalog face: a group as a standalone object, no structure.
        group::SpaceGroup const &sg = group::SpaceGroup::of(hall);
        std::println("  {} has {} Wyckoff positions, general position {}x",
                     type.international_short, sg.wyckoffs().size(),
                     sg.wyckoffs().back().multiplicity());

        // 3. Construction: a random structure with a prescribed symmetry.
        BOOST_LEAF_AUTO(generated, generate::Generator(sg, {.seed = 20260903U})(
                                       generate::Composition{{26, 2}}));
        std::println("  generated a {}-atom cell of volume {:.2f} A^3",
                     generated.cell.size(), generated.cell.lattice().volume());
        return 0;
      },
      [](leaf::error_info const &) {
        std::println("  determination failed");
        return 1;
      });
}
