/// \file test_mesh.cpp
/// \brief Mesh generation, validation diagnostics, boundary edges, sub-meshes.
#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/Selector.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"
#include "sparlab/mesh/SubMesh.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <set>

using namespace sparlab;
using Catch::Approx;

TEST_CASE("structured generator produces the documented numbering", "[mesh]") {
  StructuredMeshSpec spec;
  spec.nx = 3;
  spec.ny = 2;
  spec.lx = 3.0;
  spec.ly = 2.0;
  const Mesh mesh = make_structured_quad_mesh(spec);

  REQUIRE(mesh.num_nodes() == 12);
  REQUIRE(mesh.num_elements() == 6);
  REQUIRE(mesh.nodes_per_elem() == 4);

  // node(i, j) = j * (nx + 1) + i
  REQUIRE(mesh.node(0).isApprox(Vector2(0.0, 0.0)));
  REQUIRE(mesh.node(3).isApprox(Vector2(3.0, 0.0)));
  REQUIRE(mesh.node(4).isApprox(Vector2(0.0, 1.0)));
  REQUIRE(mesh.node(11).isApprox(Vector2(3.0, 2.0)));

  // element(0) is the lower-left cell with counter-clockwise nodes 0,1,5,4
  const Index* e0 = mesh.element_nodes(0);
  REQUIRE(e0[0] == 0);
  REQUIRE(e0[1] == 1);
  REQUIRE(e0[2] == 5);
  REQUIRE(e0[3] == 4);

  for (Index e = 0; e < mesh.num_elements(); ++e) {
    REQUIRE(mesh.element_area(e) == Approx(1.0));
  }
  REQUIRE(mesh.element_centroid(0).isApprox(Vector2(0.5, 0.5)));

  const Eigen::Vector4d bb = mesh.bounding_box();
  REQUIRE(bb(0) == Approx(0.0));
  REQUIRE(bb(1) == Approx(0.0));
  REQUIRE(bb(2) == Approx(3.0));
  REQUIRE(bb(3) == Approx(2.0));

  REQUIRE(mesh.structured_info().has_value());
  REQUIRE(mesh.structured_info()->uniform);
  REQUIRE(structured_node_index(*mesh.structured_info(), 2, 1) == 6);
  REQUIRE(structured_element_index(*mesh.structured_info(), 2, 1) == 5);
}

TEST_CASE("structured generator rejects degenerate specifications",
          "[mesh][diagnostics]") {
  StructuredMeshSpec spec;
  spec.nx = 0;
  REQUIRE_THROWS_AS(make_structured_quad_mesh(spec), ConfigError);

  spec.nx = 2;
  spec.ny = -1;
  REQUIRE_THROWS_AS(make_structured_quad_mesh(spec), ConfigError);

  spec.ny = 2;
  spec.lx = 0.0;
  REQUIRE_THROWS_AS(make_structured_quad_mesh(spec), ConfigError);

  spec.lx = 1.0;
  spec.ly = -1.0;
  REQUIRE_THROWS_AS(make_structured_quad_mesh(spec), ConfigError);
}

TEST_CASE("mesh validation catches invalid connectivity and geometry",
          "[mesh][diagnostics]") {
  Eigen::Matrix2Xd coords(2, 4);
  coords << 0, 1, 1, 0,
            0, 0, 1, 1;

  SECTION("out-of-range node index") {
    const Mesh mesh(coords, {0, 1, 2, 9}, ElementType::Quad4);
    REQUIRE_THROWS_AS(mesh.validate(), MeshError);
  }

  SECTION("clockwise ordering gives a negative area") {
    const Mesh mesh(coords, {0, 3, 2, 1}, ElementType::Quad4);
    REQUIRE_THROWS_AS(mesh.validate(), MeshError);
  }

  SECTION("repeated node inside one element") {
    const Mesh mesh(coords, {0, 1, 1, 3}, ElementType::Quad4);
    REQUIRE_THROWS_AS(mesh.validate(), MeshError);
  }

  SECTION("connectivity length is not a multiple of the stride") {
    REQUIRE_THROWS_AS(Mesh(coords, {0, 1, 2}, ElementType::Quad4), MeshError);
  }

  SECTION("orphan nodes are reported") {
    Eigen::Matrix2Xd five(2, 5);
    five << 0, 1, 1, 0, 5,
            0, 0, 1, 1, 5;
    const Mesh mesh(five, {0, 1, 2, 3}, ElementType::Quad4);
    REQUIRE_THROWS_AS(mesh.validate(), MeshError);
  }
}

TEST_CASE("boundary edges are the edges owned by one element", "[mesh]") {
  StructuredMeshSpec spec;
  spec.nx = 3;
  spec.ny = 2;
  const Mesh mesh = make_structured_quad_mesh(spec);
  const std::vector<Mesh::BoundaryEdge> edges = mesh.boundary_edges();

  // A 3 x 2 grid has 2*(3+2) = 10 boundary edges.
  REQUIRE(edges.size() == 10);
  for (const Mesh::BoundaryEdge& edge : edges) {
    REQUIRE(edge.element >= 0);
    REQUIRE(edge.element < mesh.num_elements());
    REQUIRE(edge.local_edge >= 0);
    REQUIRE(edge.local_edge < 4);
  }
}

TEST_CASE("perturbed mesh keeps the boundary and stays valid", "[mesh]") {
  StructuredMeshSpec spec;
  spec.nx = 5;
  spec.ny = 4;
  spec.lx = 2.0;
  spec.ly = 1.0;
  const Mesh uniform = make_structured_quad_mesh(spec);
  const Mesh perturbed = make_perturbed_quad_mesh(spec, 0.3, 42u);

  REQUIRE(perturbed.num_nodes() == uniform.num_nodes());
  REQUIRE(perturbed.num_elements() == uniform.num_elements());
  REQUIRE_FALSE(perturbed.structured_info()->uniform);
  REQUIRE_NOTHROW(perturbed.validate());

  // Total area is preserved exactly: only interior nodes moved.
  Scalar area_uniform = 0.0;
  Scalar area_perturbed = 0.0;
  for (Index e = 0; e < uniform.num_elements(); ++e) {
    area_uniform += uniform.element_area(e);
    area_perturbed += perturbed.element_area(e);
  }
  REQUIRE(area_perturbed == Approx(area_uniform));

  // The bounding box is unchanged.
  REQUIRE(perturbed.bounding_box().isApprox(uniform.bounding_box()));

  // Deterministic for a fixed seed.
  const Mesh again = make_perturbed_quad_mesh(spec, 0.3, 42u);
  REQUIRE(again.coordinates().isApprox(perturbed.coordinates()));

  REQUIRE_THROWS_AS(make_perturbed_quad_mesh(spec, 0.5), ConfigError);
  REQUIRE_THROWS_AS(make_perturbed_quad_mesh(spec, -0.1), ConfigError);
}

TEST_CASE("selectors pick the expected nodes and elements", "[selector]") {
  StructuredMeshSpec spec;
  spec.nx = 4;
  spec.ny = 4;
  spec.lx = 1.0;
  spec.ly = 1.0;
  const Mesh mesh = make_structured_quad_mesh(spec);

  SECTION("degenerate box selects a line of nodes") {
    SelectorGroup group;
    group.name = "left";
    Selector box;
    box.kind = SelectorKind::Box;
    box.xmax = 0.0;
    group.members.push_back(box);
    REQUIRE(group.select_nodes(mesh).size() == 5);
  }

  SECTION("circle selects elements by centroid") {
    SelectorGroup group;
    group.name = "hole";
    Selector circle;
    circle.kind = SelectorKind::Circle;
    circle.center = Vector2(0.5, 0.5);
    circle.radius = 0.2;
    group.members.push_back(circle);
    const std::vector<Index> elements = group.select_elements(mesh);
    REQUIRE(elements.size() == 4);  // the four cells touching the centre
    for (Index e : elements) {
      REQUIRE((mesh.element_centroid(e) - circle.center).norm() <= 0.2 + 1.0e-9);
    }
  }

  SECTION("annulus excludes the inner disc") {
    SelectorGroup group;
    Selector annulus;
    annulus.kind = SelectorKind::Annulus;
    annulus.center = Vector2(0.5, 0.5);
    annulus.inner_radius = 0.2;
    annulus.radius = 0.4;
    group.members.push_back(annulus);
    const std::vector<Index> elements = group.select_elements(mesh);
    REQUIRE_FALSE(elements.empty());
    for (Index e : elements) {
      const Scalar r = (mesh.element_centroid(e) - annulus.center).norm();
      REQUIRE(r >= 0.2 - 1.0e-9);
      REQUIRE(r <= 0.4 + 1.0e-9);
    }
  }

  SECTION("union of two boxes and inversion") {
    SelectorGroup group;
    Selector left;
    left.kind = SelectorKind::Box;
    left.xmax = 0.0;
    Selector right;
    right.kind = SelectorKind::Box;
    right.xmin = 1.0;
    group.members = {left, right};
    REQUIRE(group.select_nodes(mesh).size() == 10);

    group.invert = true;
    REQUIRE(group.select_nodes(mesh).size() == mesh.num_nodes() - 10);
  }

  SECTION("nearest node picks exactly one") {
    SelectorGroup group;
    Selector nearest;
    nearest.kind = SelectorKind::NearestNode;
    nearest.point = Vector2(0.99, 0.01);
    group.members.push_back(nearest);
    const std::vector<Index> nodes = group.select_nodes(mesh);
    REQUIRE(nodes.size() == 1);
    REQUIRE(mesh.node(nodes.front()).isApprox(Vector2(1.0, 0.0)));
  }

  SECTION("explicit ids are range-checked") {
    SelectorGroup group;
    group.name = "bad";
    Selector ids;
    ids.kind = SelectorKind::NodeIds;
    ids.ids = {0, 1, 9999};
    group.members.push_back(ids);
    REQUIRE_THROWS_AS(group.select_nodes(mesh), ConfigError);
  }

  SECTION("node-only selectors are rejected for element regions") {
    SelectorGroup group;
    Selector nearest;
    nearest.kind = SelectorKind::NearestNode;
    nearest.point = Vector2::Zero();
    group.members.push_back(nearest);
    REQUIRE_THROWS_AS(group.select_elements(mesh), ConfigError);
  }
}

TEST_CASE("sub-mesh extraction renumbers nodes compactly", "[mesh][submesh]") {
  StructuredMeshSpec spec;
  spec.nx = 4;
  spec.ny = 3;
  const Mesh mesh = make_structured_quad_mesh(spec);

  const std::vector<Index> subset = {0, 1, 4, 5};
  const SubMeshResult sub = extract_element_subset(mesh, subset);

  REQUIRE(sub.mesh.num_elements() == 4);
  REQUIRE(sub.element_map == subset);
  REQUIRE(static_cast<std::size_t>(sub.mesh.num_nodes()) == sub.node_map.size());
  REQUIRE_NOTHROW(sub.mesh.validate());

  // Coordinates survive the renumbering.
  for (Index n = 0; n < sub.mesh.num_nodes(); ++n) {
    REQUIRE(sub.mesh.node(n).isApprox(mesh.node(sub.node_map[static_cast<std::size_t>(n)])));
  }
  // old_to_new is consistent with node_map.
  for (std::size_t n = 0; n < sub.node_map.size(); ++n) {
    REQUIRE(sub.old_to_new_node[static_cast<std::size_t>(sub.node_map[n])] ==
            static_cast<Index>(n));
  }
  REQUIRE_THROWS_AS(extract_element_subset(mesh, {}), MeshError);
  REQUIRE_THROWS_AS(extract_element_subset(mesh, {999}), MeshError);
}

TEST_CASE("edge connectivity ignores corner-only contact", "[mesh][submesh]") {
  StructuredMeshSpec spec;
  spec.nx = 2;
  spec.ny = 2;
  const Mesh mesh = make_structured_quad_mesh(spec);

  // A checkerboard pair touching only at one corner must count as two groups.
  Vector density = Vector::Zero(mesh.num_elements());
  density(0) = 1.0;  // lower-left
  density(3) = 1.0;  // upper-right
  const Vector volumes = Vector::Ones(mesh.num_elements());

  const TopologyInterpretation report =
      interpret_density_as_solid(mesh, density, volumes, 0.5,
                                 /*largest_component_only=*/true);
  REQUIRE(report.elements_above_threshold == 2);
  REQUIRE(report.components_above_threshold == 2);
  REQUIRE(report.elements_retained == 1);
  REQUIRE(report.volume_discarded_as_islands == Approx(1.0));

  // Keeping everything reports two groups but retains both cells.
  const TopologyInterpretation keep_all =
      interpret_density_as_solid(mesh, density, volumes, 0.5, false);
  REQUIRE(keep_all.elements_retained == 2);
  REQUIRE(keep_all.volume_discarded_as_islands == Approx(0.0));
}

TEST_CASE("density interpretation rejects impossible thresholds",
          "[mesh][submesh][diagnostics]") {
  StructuredMeshSpec spec;
  spec.nx = 2;
  spec.ny = 2;
  const Mesh mesh = make_structured_quad_mesh(spec);
  const Vector volumes = Vector::Ones(mesh.num_elements());

  REQUIRE_THROWS_AS(
      interpret_density_as_solid(mesh, Vector::Zero(mesh.num_elements()), volumes, 0.5),
      MeshError);
  REQUIRE_THROWS_AS(
      interpret_density_as_solid(mesh, Vector::Ones(mesh.num_elements()), volumes, 0.0),
      ConfigError);
  REQUIRE_THROWS_AS(
      interpret_density_as_solid(mesh, Vector::Ones(mesh.num_elements()), volumes, 1.0),
      ConfigError);
  REQUIRE_THROWS_AS(interpret_density_as_solid(mesh, Vector::Ones(2), volumes, 0.5),
                    MeshError);
}

TEST_CASE("element components are found for a fully connected mesh",
          "[mesh][submesh]") {
  StructuredMeshSpec spec;
  spec.nx = 5;
  spec.ny = 4;
  const Mesh mesh = make_structured_quad_mesh(spec);
  const auto components = element_components_by_edge(mesh);
  REQUIRE(components.size() == 1);
  REQUIRE(components.front().size() == static_cast<std::size_t>(mesh.num_elements()));
}
