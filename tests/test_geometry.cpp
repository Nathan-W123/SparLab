/// \file test_geometry.cpp
/// \brief Before/after geometry export: boundary surfaces, STL layout and the
///        watertightness check that ties the surface back to the cells.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/io/Config.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/io/Json.hpp"
#include "sparlab/io/ResultWriter.hpp"
#include "sparlab/io/StlWriter.hpp"
#include "sparlab/mesh/SubMesh.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cstdio>
#include <fstream>
#include <sys/stat.h>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;

namespace {

long file_size(const std::string& path) {
  struct stat st {};
  if (stat(path.c_str(), &st) != 0) return -1;
  return static_cast<long>(st.st_size);
}

/// Every triangle's outward normal must point away from the centre of a
/// convex solid.
void require_outward(const TriangleSurface& surface, const Vector3& centre) {
  for (const Triangle& t : surface.triangles()) {
    const Vector3 centroid = (t.a + t.b + t.c) / 3.0;
    REQUIRE(t.normal().dot(centroid - centre) > 0.0);
  }
}

}  // namespace

TEST_CASE("a 2-D mesh is extruded into a closed surface of the right volume",
          "[geometry][stl][verification]") {
  StructuredMeshSpec spec;
  spec.nx = 6;
  spec.ny = 3;
  spec.lx = 0.6;
  spec.ly = 0.3;
  const Mesh mesh = make_perturbed_quad_mesh(spec, 0.25, 5u);
  const Scalar thickness = 0.02;
  const TriangleSurface surface = boundary_surface(mesh, thickness);

  // 2 caps x 2 triangles per cell + 2 triangles per boundary edge.
  REQUIRE(surface.faces.size() == static_cast<std::size_t>(4 * mesh.num_elements() +
                                                           2 * mesh.boundary_faces().size()));
  REQUIRE(surface.points.size() == static_cast<std::size_t>(2 * mesh.num_nodes()));
  const SurfaceStats stats = surface_stats(surface);
  REQUIRE(stats.closed);
  REQUIRE(stats.unmatched_edges == 0);
  REQUIRE(stats.enclosed_volume == Approx(0.6 * 0.3 * thickness).epsilon(1.0e-12));
  REQUIRE(stats.area == Approx(2.0 * 0.6 * 0.3 + 2.0 * (0.6 + 0.3) * thickness)
                            .epsilon(1.0e-12));
  REQUIRE(stats.bounds.lower.isApprox(Vector3::Zero()));
  REQUIRE(stats.bounds.upper.isApprox(Vector3(0.6, 0.3, thickness)));
  require_outward(surface, Vector3(0.3, 0.15, 0.5 * thickness));
  REQUIRE_THROWS_AS(boundary_surface(mesh, 0.0), ConfigError);
}

TEST_CASE("a hex mesh's boundary faces form a closed outward surface",
          "[geometry][stl][verification]") {
  StructuredMeshSpec spec;
  spec.nx = 4;
  spec.ny = 3;
  spec.nz = 2;
  spec.lx = 0.4;
  spec.ly = 0.3;
  spec.lz = 0.2;
  const Mesh mesh = make_perturbed_hex_mesh(spec, 0.2, 9u);
  const TriangleSurface surface = boundary_surface(mesh);
  REQUIRE(surface.faces.size() == 2 * mesh.boundary_faces().size());
  const SurfaceStats stats = surface_stats(surface);
  REQUIRE(stats.closed);
  // The outer faces of a perturbed grid are planar (boundary nodes do not
  // move), so the enclosed volume is exact.
  REQUIRE(stats.enclosed_volume == Approx(0.4 * 0.3 * 0.2).epsilon(1.0e-12));
  REQUIRE(stats.area == Approx(2.0 * (0.4 * 0.3 + 0.3 * 0.2 + 0.2 * 0.4)).epsilon(1.0e-12));
  require_outward(surface, Vector3(0.2, 0.15, 0.1));

  // A thresholded, non-convex structure (an L of cells) is still watertight
  // and encloses exactly its cells.
  Vector density = Vector::Zero(mesh.num_elements());
  const StructuredGridInfo& info = *mesh.structured_info();
  for (Index k = 0; k < info.nz; ++k) {
    for (Index i = 0; i < info.nx; ++i) density(structured_element_index(info, i, 0, k)) = 1.0;
    for (Index j = 0; j < info.ny; ++j) density(structured_element_index(info, 0, j, k)) = 1.0;
  }
  Vector volumes(mesh.num_elements());
  for (Index e = 0; e < mesh.num_elements(); ++e) volumes(e) = mesh.element_measure(e);
  const TopologyInterpretation interp =
      interpret_density_as_solid(mesh, density, volumes, 0.5);
  REQUIRE(interp.components_above_threshold == 1);
  const TriangleSurface after = boundary_surface(interp.sub.mesh);
  const SurfaceStats after_stats = surface_stats(after);
  REQUIRE(after_stats.closed);
  // The cut runs through distorted cells whose faces are bilinear patches;
  // two flat triangles approximate each, so the enclosed volume agrees with
  // the cell volume only to the size of that geometric error (measured here
  // at 1.5e-4), not to round-off. An open or inverted surface would be off by
  // orders of magnitude more.
  REQUIRE(after_stats.enclosed_volume == Approx(interp.volume_retained).epsilon(1.0e-3));
  REQUIRE(after_stats.enclosed_volume != Approx(interp.volume_retained).epsilon(1.0e-9));

  // On the uniform grid the same cut is exact.
  const Mesh uniform = make_structured_hex_mesh(spec);
  Vector uniform_volumes(uniform.num_elements());
  for (Index e = 0; e < uniform.num_elements(); ++e) {
    uniform_volumes(e) = uniform.element_measure(e);
  }
  const TopologyInterpretation exact =
      interpret_density_as_solid(uniform, density, uniform_volumes, 0.5);
  REQUIRE(surface_stats(boundary_surface(exact.sub.mesh)).enclosed_volume ==
          Approx(exact.volume_retained).epsilon(1.0e-12));

  // Dropping one triangle opens the surface, and the check says so.
  TriangleSurface open = after;
  open.faces.pop_back();
  const SurfaceStats open_stats = surface_stats(open);
  REQUIRE_FALSE(open_stats.closed);
  REQUIRE(open_stats.unmatched_edges == 3);
}

TEST_CASE("binary STL files have the standard layout and round-trip",
          "[geometry][stl][io]") {
  ensure_directory("results/_test_tmp");
  const std::string path = "results/_test_tmp/box.stl";
  StructuredMeshSpec spec;
  spec.nx = spec.ny = spec.nz = 2;
  const Mesh mesh = make_structured_hex_mesh(spec);
  const TriangleSurface surface = boundary_surface(mesh);
  write_stl(path, surface, "unit test");

  // 80-byte header, 4-byte count, 50 bytes per facet.
  REQUIRE(file_size(path) == 84 + 50 * static_cast<long>(surface.faces.size()));
  std::ifstream in(path, std::ios::binary);
  char head[80];
  in.read(head, 80);
  REQUIRE(std::string(head, 5) != "solid");
  REQUIRE(std::string(head, 20) == "SparLab binary STL: ");

  const std::vector<Triangle> back = read_stl(path);
  const std::vector<Triangle> expected = surface.triangles();
  REQUIRE(back.size() == expected.size());
  // Single precision round trip of coordinates that are multiples of 0.5.
  Scalar volume = 0.0;
  for (std::size_t i = 0; i < back.size(); ++i) {
    REQUIRE(back[i].a.isApprox(expected[i].a));
    REQUIRE(back[i].b.isApprox(expected[i].b));
    REQUIRE(back[i].c.isApprox(expected[i].c));
    volume += back[i].a.dot(back[i].b.cross(back[i].c)) / 6.0;
  }
  REQUIRE(volume == Approx(1.0));
  REQUIRE(surface_stats(surface).enclosed_volume == Approx(1.0));
  std::remove(path.c_str());
  REQUIRE_THROWS_AS(write_stl("/definitely/not/a/dir/x.stl", surface), IoError);
  REQUIRE_THROWS_AS(read_stl("/definitely/not/here.stl"), IoError);
}

TEST_CASE("the result writer exports before/after geometry with its checks",
          "[geometry][io][writers]") {
  const json::Value doc = json::parse(R"({
    "name": "geometry_case",
    "mesh": { "nx": 10, "ny": 5, "lx": 0.5, "ly": 0.25 },
    "material": { "youngs_modulus": 70e9, "poisson_ratio": 0.3, "density": 2700 },
    "model": { "thickness": 0.01 },
    "boundary_conditions": [ { "fix": ["x","y"], "region": { "box": { "xmax": 0.0 } } } ],
    "load_cases": [ { "name": "tip",
      "point_loads": [ { "force": [0,-100], "region": { "box": { "xmin": 0.5, "ymax": 0.0 } } } ] } ],
    "topology": { "enabled": true, "volume_fraction": 0.5,
                  "optimizer": { "max_iterations": 3 } }
  })", "geometry");
  const Configuration config = parse_configuration(doc, "geometry", true);
  FemModel model = build_model(config);
  const DesignDomain domain = build_design_domain(config, model);

  const std::string dir = "results/_test_tmp/geometry";
  ResultWriter writer(dir, config);
  const json::Value before = writer.write_geometry(model.mesh(), domain.initial_design(),
                                                   model.thickness(), "structure_before",
                                                   "design domain");
  REQUIRE(before.find("num_triangles")->number_value() ==
          Approx(4.0 * model.mesh().num_elements() + 2.0 * 30.0));
  REQUIRE(before.find("enclosed_volume_m3")->number_value() ==
          Approx(0.5 * 0.25 * 0.01).epsilon(1.0e-12));
  REQUIRE(before.find("volume_relative_mismatch")->number_value() < 1.0e-12);
  REQUIRE(before.find("extruded_thickness_m")->number_value() == Approx(0.01));
  REQUIRE(file_size(path_join(dir, "structure_before.stl")) > 84);
  REQUIRE(file_size(path_join(dir, "structure_before.vtk")) > 0);

  // A thresholded half-height strip as the "after" structure.
  // Three of the five element rows (cells are 0.05 m tall).
  Vector density = Vector::Zero(model.mesh().num_elements());
  for (Index e = 0; e < model.mesh().num_elements(); ++e) {
    if (model.mesh().element_centroid(e).y() < 0.13) density(e) = 0.9;
  }
  const TopologyInterpretation interp = interpret_density_as_solid(
      model.mesh(), density, domain.element_volumes(), 0.5);
  Vector retained(interp.sub.mesh.num_elements());
  for (std::size_t i = 0; i < interp.sub.element_map.size(); ++i) {
    retained(static_cast<Eigen::Index>(i)) = density(interp.sub.element_map[i]);
  }
  const json::Value after = writer.write_geometry(interp.sub.mesh, retained,
                                                  model.thickness(), "structure_after",
                                                  "interpreted structure");
  REQUIRE(after.find("enclosed_volume_m3")->number_value() ==
          Approx(interp.volume_retained).epsilon(1.0e-12));
  REQUIRE(after.find("cell_volume_m3")->number_value() ==
          Approx(0.5 * 0.15 * 0.01).epsilon(1.0e-12));
  REQUIRE(after.find("volume_relative_mismatch")->number_value() < 1.0e-12);
  REQUIRE(after.find("closed_surface")->bool_value());
  REQUIRE(before.find("closed_surface")->bool_value());

  // Length mismatches are refused.
  REQUIRE_THROWS_AS(writer.write_geometry(model.mesh(), Vector::Ones(3), 0.01, "x", "y"),
                    IoError);
}
