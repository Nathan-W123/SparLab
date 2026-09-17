/// \file TestSupport.hpp
/// \brief Shared fixtures for the SparLab test suite.
#pragma once

#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/FemModel.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"

#include <cmath>

namespace sparlab {
namespace testing {

/// Reference aluminium-like material used by most tests.
inline IsotropicMaterial default_material(Scalar poisson = 0.3) {
  return IsotropicMaterial(70.0e9, poisson, 2700.0, "test_aluminium");
}

/// Geometry and loading of the cantilever used by the beam-theory tests.
struct CantileverCase {
  Scalar length = 1.0;
  Scalar height = 0.1;
  Scalar thickness = 0.01;
  Scalar youngs = 70.0e9;
  Scalar poisson = 0.0;
  Scalar density = 2700.0;
  Scalar tip_load = -1000.0;

  /// Euler-Bernoulli tip deflection [m].
  Scalar euler_bernoulli_tip() const {
    const Scalar i = thickness * height * height * height / 12.0;
    return tip_load * length * length * length / (3.0 * youngs * i);
  }

  /// Timoshenko tip deflection (rectangular shear coefficient 5/6) [m].
  Scalar timoshenko_tip() const {
    const Scalar i = thickness * height * height * height / 12.0;
    const Scalar a = thickness * height;
    const Scalar g = youngs / (2.0 * (1.0 + poisson));
    return tip_load * length * length * length / (3.0 * youngs * i) +
           tip_load * length / ((5.0 / 6.0) * g * a);
  }
};

/// Build the standard clamped-root cantilever with a tip resultant spread over
/// the tip edge nodes.
inline FemModel make_cantilever(const CantileverCase& c, Index nx, Index ny) {
  StructuredMeshSpec spec;
  spec.nx = nx;
  spec.ny = ny;
  spec.lx = c.length;
  spec.ly = c.height;

  FemModel model(make_structured_quad_mesh(spec),
                 IsotropicMaterial(c.youngs, c.poisson, c.density, "cantilever"),
                 c.thickness, StressState::PlaneStress, IntegrationOptions());

  DisplacementConstraint root;
  root.region.name = "root";
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  root.region.members.push_back(box);
  root.fix_x = true;
  root.fix_y = true;
  model.constraints().push_back(root);

  LoadCaseSpec load;
  load.name = "tip";
  PointLoadSpec tip;
  tip.region.name = "tip_edge";
  Selector tip_box;
  tip_box.kind = SelectorKind::Box;
  tip_box.xmin = c.length;
  tip.region.members.push_back(tip_box);
  tip.force = Vector2(0.0, c.tip_load);
  tip.distribute_total = true;
  load.point_loads.push_back(tip);
  model.load_case_specs().push_back(load);

  model.finalize();
  return model;
}

/// Mean vertical displacement over the tip edge [m].
inline Scalar tip_deflection(const FemModel& model, const Vector& u) {
  const StructuredGridInfo& info = *model.mesh().structured_info();
  Scalar sum = 0.0;
  for (Index j = 0; j <= info.ny; ++j) {
    const Index node = structured_node_index(info, info.nx, j);
    sum += u(node * kDofsPerNode + 1);
  }
  return sum / static_cast<Scalar>(info.ny + 1);
}

/// A tiny square plate clamped on its left edge, used for solver comparisons.
inline FemModel make_small_plate(Index nx = 4, Index ny = 3, Scalar poisson = 0.3) {
  StructuredMeshSpec spec;
  spec.nx = nx;
  spec.ny = ny;
  spec.lx = 0.4;
  spec.ly = 0.3;

  FemModel model(make_structured_quad_mesh(spec), default_material(poisson), 0.005,
                 StressState::PlaneStress, IntegrationOptions());

  DisplacementConstraint bc;
  bc.region.name = "left";
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  bc.region.members.push_back(box);
  bc.fix_x = true;
  bc.fix_y = true;
  model.constraints().push_back(bc);

  LoadCaseSpec load;
  load.name = "corner";
  PointLoadSpec p;
  p.region.name = "corner_node";
  Selector nearest;
  nearest.kind = SelectorKind::NearestNode;
  nearest.point = Vector2(spec.lx, 0.0);
  p.region.members.push_back(nearest);
  p.force = Vector2(500.0, -2000.0);
  load.point_loads.push_back(p);
  model.load_case_specs().push_back(load);

  model.finalize();
  return model;
}

/// The three infinitesimal rigid-body modes of a node set, about its centroid.
/// Columns: x-translation, y-translation, rotation.
inline Matrix rigid_body_modes(const Mesh& mesh) {
  const Index nn = mesh.num_nodes();
  Vector2 centroid = Vector2::Zero();
  for (Index n = 0; n < nn; ++n) centroid += mesh.node(n);
  centroid /= static_cast<Scalar>(nn);

  Matrix modes = Matrix::Zero(nn * kDofsPerNode, 3);
  for (Index n = 0; n < nn; ++n) {
    const Vector2 r = mesh.node(n) - centroid;
    modes(n * kDofsPerNode + 0, 0) = 1.0;
    modes(n * kDofsPerNode + 1, 1) = 1.0;
    modes(n * kDofsPerNode + 0, 2) = -r.y();
    modes(n * kDofsPerNode + 1, 2) = r.x();
  }
  return modes;
}

/// Nodal coordinates of a single unit square Q4, counter-clockwise.
inline Eigen::Matrix<Scalar, 2, Eigen::Dynamic> unit_square_coords(Scalar a = 1.0,
                                                                   Scalar b = 1.0) {
  Eigen::Matrix<Scalar, 2, Eigen::Dynamic> coords(2, 4);
  coords << 0.0, a, a, 0.0,
            0.0, 0.0, b, b;
  return coords;
}

/// A distorted (but convex) quadrilateral for isoparametric-mapping tests.
inline Eigen::Matrix<Scalar, 2, Eigen::Dynamic> distorted_quad_coords() {
  Eigen::Matrix<Scalar, 2, Eigen::Dynamic> coords(2, 4);
  coords << 0.10, 2.05, 1.80, 0.35,
            0.05, 0.20, 1.35, 1.10;
  return coords;
}

}  // namespace testing
}  // namespace sparlab
