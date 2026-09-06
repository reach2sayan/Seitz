#include <seitz/core/keys.hpp>
#include <seitz/core/periodicity.hpp>
#include <seitz/data/spg_database.hpp>
#include <seitz/group/space_group.hpp>
#include <seitz/group/subgroup_graph.hpp>
#include <seitz/group/wyckoff.hpp>

#include "casters.hpp" // borrowed_list, copied_list, to_str
#include "errors.hpp"  // unwrap

#include <pybind11/eigen.h>
#include <pybind11/native_enum.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <vector>

namespace seitz::python {

namespace {

using group::GroupBase;
using group::SpaceGroup;
using group::SubgroupGraph;
using group::SubgroupEdge;
using group::SubgroupKind;
using group::Wyckoff;

// The GroupBase face, defined once and mixed into each family's class_. The
// families are unrelated concrete types that happen to answer the same
// questions -- there is no runtime hierarchy in the C++ and there is none here.
template <class G> void bind_group_base(py::class_<G> &cls) {
  // Lambdas, not member pointers: every GroupBase accessor is noexcept, part of
  // the type since C++17, and pybind11's method_adaptor has no noexcept
  // overload -- a member pointer falls through to the pass-through overload and
  // binds GroupBase, an unregistered type, as self. It fails at call time.

  cls.def_property_readonly("number",
                            [](G const &self) { return self.number(); })
      .def_property_readonly(
          "symbol", [](G const &self) { return to_str(self.symbol()); })
      .def_property_readonly(
          "order", [](G const &self) { return self.order(); },
          "Number of operations in the group.")
      .def_property_readonly(
          "operations",
          [](G const &self) { return copied_list(self.operations()); })
      // Borrowed, not copied: a Wyckoff's ADDRESS is its identity --
      // generate::Placed stores a Wyckoff const * into exactly this storage --
      // so each element is parented on the group and keeps it alive.
      .def_property_readonly(
          "wyckoffs",
          [](py::object const &self) {
            return borrowed_list(self, self.cast<G const &>().wyckoffs());
          },
          py::doc("The Wyckoff positions, ascending by letter; the last is the "
                  "general position."))
      .def(
          "wyckoff",
          [](py::object const &self, char letter) {
            Wyckoff const *w =
                unwrap([&] { return self.cast<G const &>().wyckoff(letter); });
            return py::cast(*w, py::return_value_policy::reference_internal,
                            self);
          },
          py::arg("letter"),
          py::doc(
              "The position with that letter. Raises if this group has none."))
      .def("__len__", [](G const &self) { return self.wyckoffs().size(); });
}

} // namespace

void bind_group(py::module_ &m) {
  // ---- Wyckoff -----------------------------------------------------------
  //
  // No constructor: Wyckoff's is private, and instances only ever come from a
  // group, which is also what owns them.
  py::class_<Wyckoff>(
      m, "Wyckoff",
      "One Wyckoff position: an orbit class parameterised by an "
      "affine locus in the fractional cell.")
      .def_property_readonly("multiplicity", &Wyckoff::multiplicity)
      .def_property_readonly("degrees_of_freedom", &Wyckoff::degrees_of_freedom,
                             "Free coordinates of the locus, 0..3.")
      .def_property_readonly(
          "letter",
          [](Wyckoff const &self) { return std::string(1, self.letter()); },
          "Wyckoff letter; 'a' is the most special.")
      .def_property_readonly(
          "site_symmetry",
          [](Wyckoff const &self) { return to_str(self.site_symmetry()); },
          "Tabulated symbol; empty for a derived (point/rod) position.")
      .def_property_readonly(
          "operations",
          [](Wyckoff const &self) { return copied_list(self.operations()); },
          "The site-symmetry group: the operations fixing a generic point.")
      .def(
          "sample",
          [](Wyckoff const &self, std::vector<double> const &params) {
            return self.sample(std::span<double const>{params});
          },
          py::arg("params"),
          py::doc("A point on the locus from `degrees_of_freedom` free "
                  "parameters. Extra parameters are ignored."))
      .def("canonical", &Wyckoff::canonical, py::arg("xyz"),
           py::doc("Project a point onto the locus, so an approximate "
                   "coordinate still yields the exact orbit."))
      .def("orbit", &Wyckoff::orbit, py::arg("xyz"),
           py::arg_v("periodicity", all_periodic(), "all_periodic()"),
           py::doc("The full orbit of `xyz` as an (N, 3) array, one row per "
                   "image, projected onto the locus first."))
      .def("__repr__", [](Wyckoff const &self) {
        return "Wyckoff('" + std::string(1, self.letter()) +
               "', multiplicity=" + std::to_string(self.multiplicity()) + ")";
      });

  // ---- SpaceGroup --------------------------------------------------------
  //
  // Layer groups come through this class too, the family carried by the key --
  // there is no LayerGroup type, in Python or in C++.
  //
  // `explicit SpaceGroup(HallNumber)` is deliberately not bound: it builds a
  // private unshared instance, i.e. a duplicate of what the flyweight shares.
  py::class_<SpaceGroup> space_group(
      m, "SpaceGroup",
      "A space or layer group as a standalone, structure-free object. One "
      "shared immutable instance per Hall setting.");
  bind_group_base(space_group);
  space_group
      .def_static("of", &SpaceGroup::of, py::arg("hall"),
                  // A Boost.Flyweight: one immutable object per setting, alive
                  // for the life of the process. `reference` with no parent is
                  // the right policy -- there is no owner to keep alive.
                  py::return_value_policy::reference,
                  py::doc("The group of a Hall setting. Total: a HallNumber "
                          "cannot name a setting that does not exist."))
      .def_static(
          "from_number",
          [](GroupFamily family, int number) {
            return unwrap(
                [&] { return SpaceGroup::from_number(family, number); });
          },
          py::arg("family"), py::arg("number"),
          py::return_value_policy::reference,
          py::doc("The group of an international number, in its default Hall "
                  "setting. 1..230 for space groups, 1..80 for layer groups."))
      .def_property_readonly("hall", &SpaceGroup::hall)
      .def_property_readonly("type", &SpaceGroup::type,
                             py::return_value_policy::reference)
      .def_property_readonly("centering", &SpaceGroup::centering)
      .def("__repr__", [](SpaceGroup const &self) {
        return "SpaceGroup(" + std::to_string(self.number()) + ", '" +
               std::string(self.symbol()) + "')";
      });

  // ---- the subgroup graph ------------------------------------------------
  //
  // All-static and all-constexpr on the C++ side, so it is bound as free
  // functions in a submodule rather than a class nobody would instantiate.
  py::native_enum<SubgroupKind>(
      m, "SubgroupKind", "enum.IntEnum",
      "How a maximal subgroup sits in its supergroup: the same lattice and a "
      "smaller point group (translationengleiche), or the same point group and "
      "a smaller lattice (klassengleiche).")
      .value("translationengleiche", SubgroupKind::translationengleiche)
      .value("klassengleiche", SubgroupKind::klassengleiche)
      .finalize();

  py::class_<SubgroupEdge>(
      m, "SubgroupEdge",
      "One maximal-subgroup relation: `sub` is a maximal subgroup of `super` "
      "(space-group numbers; `sub` in the Hall setting `hall`) of the given "
      "kind and index. The subgroup's conventional cell in the supergroup's "
      "frame: (a_H b_H c_H) = (a_G b_G c_G) @ basis, x_G = basis @ x_H + "
      "origin.")
      .def_readonly("id", &SubgroupEdge::id)
      .def_readonly("super", &SubgroupEdge::super)
      .def_readonly("sub", &SubgroupEdge::sub)
      .def_readonly("hall", &SubgroupEdge::hall)
      .def_readonly("kind", &SubgroupEdge::kind)
      .def_readonly("index", &SubgroupEdge::index)
      .def_readonly("basis", &SubgroupEdge::basis)
      .def_readonly("origin", &SubgroupEdge::origin)
      .def("in_subgroup_frame", &SubgroupEdge::in_subgroup_frame, py::arg("op"),
           py::doc("An operation of the supergroup's frame expressed in the "
                   "subgroup's, or None if it does not map the subgroup's "
                   "lattice onto itself."))
      .def("__repr__", [](SubgroupEdge const &self) {
        return "SubgroupEdge(" + std::to_string(self.super) + " -> " +
               std::to_string(self.sub) + ", kind=" +
               (self.kind == SubgroupKind::translationengleiche ? "t" : "k") +
               ", index=" + std::to_string(self.index) + ")";
      });

  py::module_ subgroups = m.def_submodule(
      "subgroups",
      "The maximal-subgroup graph of the 230 space groups, translationengleiche "
      "and klassengleiche.");
  auto const as_edges = [](std::ranges::input_range auto ids) {
    return std::vector<SubgroupEdge>{
        std::from_range, ids | std::views::transform(SubgroupGraph::edge)};
  };
  subgroups.def(
      "maximal_subgroups",
      [as_edges](int number, std::optional<SubgroupKind> kind) {
        return as_edges(SubgroupGraph::maximal_subgroups(number, kind));
      },
      py::arg("number"), py::arg("kind") = py::none(),
      py::doc("The maximal subgroups of `number`, of one kind or both. Empty "
              "out of range."));
  subgroups.def(
      "minimal_supergroups",
      [as_edges](int number, std::optional<SubgroupKind> kind) {
        return as_edges(SubgroupGraph::minimal_supergroups(number, kind));
      },
      py::arg("number"), py::arg("kind") = py::none(),
      py::doc("The groups of which `number` is a maximal subgroup."));
  subgroups.def("edge", &SubgroupGraph::edge, py::arg("id"),
                py::doc("The relation with that id."));
  subgroups.def("is_subgroup", &SubgroupGraph::is_subgroup, py::arg("sub"),
                py::arg("super"),
                py::doc("Whether `sub` is reachable from `super` by a chain of "
                        "maximal-subgroup steps of any kind. True when they are "
                        "equal."));
  subgroups.def(
      "path",
      [as_edges](int super, int sub,
                 std::optional<SubgroupKind> kind) -> py::object {
        auto const chain = SubgroupGraph::path(super, sub, kind);
        if (!chain) {
          return py::none();
        }
        return py::cast(as_edges(*chain));
      },
      py::arg("super"), py::arg("sub"), py::arg("kind") = py::none(),
      py::doc("A shortest chain of relations from `super` down to `sub` along "
              "steps of the given kind (both when None); empty when they are "
              "equal, None when there is none."));
  subgroups.attr("K_NUM_SPACE_GROUPS") = group::kNumSpaceGroups;
}

} // namespace seitz::python
