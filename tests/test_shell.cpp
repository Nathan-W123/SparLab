/// \file test_shell.cpp
/// \brief The MITC4 shell element: rigid-body motions, the zero-energy
///        modes, mass and rotary inertia, loads, and the geometric
///        stiffness and its derivative, on flat, rotated and warped cells;
///        and shell models: the membrane and bending patch tests, a curved
///        panel moved rigidly, directors at folds, assembled loads and a
///        plate without shear locking.
#include "TestSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/elements/Shell4.hpp"
#include "sparlab/fem/Loads.hpp"
#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/io/CalculixWriter.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/io/Config.hpp"
#include "sparlab/io/MeshReader.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <sstream>

using namespace sparlab;
using namespace sparlab::testing;
using Catch::Approx;
using Catch::Matchers::ContainsSubstring;

namespace {

/// A rectangle a x b in the x-y plane.
Matrix flat_cell(Scalar a, Scalar b) {
  Matrix x(3, 4);
  x << 0.0, a, a, 0.0,
       0.0, 0.0, b, b,
       0.0, 0.0, 0.0, 0.0;
  return x;
}

/// A distorted, warped cell turned out of every coordinate plane.
Matrix warped_cell() {
  Matrix x(3, 4);
  x << 0.00, 1.05, 0.95, -0.10,
       0.00, 0.10, 0.85, 0.75,
       0.00, 0.06, -0.04, 0.08;
  const Matrix3 rot = (Eigen::AngleAxisd(0.4, Vector3::UnitX()) *
                       Eigen::AngleAxisd(-0.7, Vector3::UnitY()) *
                       Eigen::AngleAxisd(0.3, Vector3::UnitZ()))
                          .toRotationMatrix();
  return rot * x;
}

/// The geometry with directors: each node's normal tilted a little (as
/// averaging over neighbours does), so they are not the element's own.
Matrix with_tilted_directors(const Matrix& x) {
  Matrix g(6, 4);
  g.topRows(3) = x;
  for (int k = 0; k < 4; ++k) {
    const Scalar r = k == 0 || k == 3 ? -1.0 : 1.0;
    const Scalar s = k < 2 ? -1.0 : 1.0;
    Vector3 n = Shell4Element::normal(x, r, s);
    n += 0.05 * Vector3(std::cos(1.0 + k), std::sin(2.0 * k), 0.3);
    g.block<3, 1>(3, k) = n.normalized();
  }
  return g;
}

/// The six rigid-body motions of the element's nodes (translations, then
/// rotations about the global axes through the origin).
Matrix rigid_modes(const Matrix& geometry) {
  Matrix r = Matrix::Zero(24, 6);
  for (int k = 0; k < 4; ++k) {
    const Vector3 x = geometry.block<3, 1>(0, k);
    for (int m = 0; m < 3; ++m) {
      r(6 * k + m, m) = 1.0;
      Vector3 w = Vector3::Zero();
      w(m) = 1.0;
      r.block<3, 1>(6 * k, 3 + m) = w.cross(x);
      r.block<3, 1>(6 * k + 3, 3 + m) = w;
    }
  }
  return r;
}

}  // namespace

TEST_CASE("the shell stiffness leaves the six rigid-body motions free and no other",
          "[shell][element]") {
  const Matrix3 d = default_material().plane_stress_matrix();
  const IntegrationOptions opts;
  for (const Matrix& geometry : {flat_cell(1.0, 0.6), warped_cell(),
                                 with_tilted_directors(warped_cell())}) {
    for (const Scalar t : {0.1, 0.001}) {
      INFO("rows " << geometry.rows() << ", t = " << t);
      const Shell4Element element;
      const Matrix ke = element.stiffness(geometry, d, t, opts);
      REQUIRE((ke - ke.transpose()).cwiseAbs().maxCoeff() <= 1e-15 * ke.cwiseAbs().maxCoeff());
      const Matrix r = rigid_modes(geometry);
      REQUIRE((ke * r).norm() <= 1e-12 * ke.norm() * r.norm());
      // Exactly six zero-energy modes: MITC4 has no spurious ones, and the
      // drilling penalty takes the rotations about the normal.
      Eigen::SelfAdjointEigenSolver<Matrix> eig(ke);
      const Vector lambda = eig.eigenvalues();
      const Scalar top = lambda.maxCoeff();
      int zero = 0;
      for (Eigen::Index i = 0; i < lambda.size(); ++i) {
        REQUIRE(lambda(i) >= -1e-12 * top);
        if (lambda(i) < 1e-12 * top) ++zero;
      }
      REQUIRE(zero == 6);
    }
  }
  // Without the drilling penalty the rotation about the normal is free at
  // every node of a flat cell: four more zero-energy modes.
  const Shell4Element free(0.0);
  Eigen::SelfAdjointEigenSolver<Matrix> eig(free.stiffness(flat_cell(1.0, 0.6), d, 0.1, opts));
  const Vector lambda = eig.eigenvalues();
  int zero = 0;
  for (Eigen::Index i = 0; i < lambda.size(); ++i) zero += lambda(i) < 1e-12 * lambda.maxCoeff();
  REQUIRE(zero == 10);
}

TEST_CASE("the shell mass holds rho t A per direction and the rotary inertia rho t^3 A / 12",
          "[shell][element]") {
  const Scalar rho = 2700.0;
  const IntegrationOptions opts;
  const Shell4Element element;
  const Scalar a = 1.2;
  const Scalar b = 0.7;
  const Scalar t = 0.05;
  const Matrix me = element.consistent_mass(flat_cell(a, b), rho, t, opts);
  REQUIRE((me - me.transpose()).cwiseAbs().maxCoeff() <= 1e-15 * me.cwiseAbs().maxCoeff());
  const Matrix r = rigid_modes(flat_cell(a, b));
  for (int m = 0; m < 3; ++m) {
    REQUIRE(r.col(m).dot(me * r.col(m)) == Approx(rho * t * a * b).epsilon(1e-13));
  }
  // A uniform rotation about x or y turns the thickness through the plane.
  for (int m = 0; m < 2; ++m) {
    Vector theta = Vector::Zero(24);
    for (int k = 0; k < 4; ++k) theta(6 * k + 3 + m) = 1.0;
    REQUIRE(theta.dot(me * theta) == Approx(rho * t * t * t / 12.0 * a * b).epsilon(1e-13));
  }
  // On a warped cell the thickness runs along the interpolated director,
  // which is not quite normal to the surface: the translational mass is rho
  // times the volume of that solid, X = sum N_k (x_k + zeta t/2 V_k) -
  // measured here by central differences of X on a 6 x 6 x 6 Gauss grid.
  const Matrix w = warped_cell();
  const Matrix mw = element.consistent_mass(w, rho, t, opts);
  const Matrix rw = rigid_modes(w);
  Matrix v(3, 4);
  for (int k = 0; k < 4; ++k) {
    v.col(k) = Shell4Element::normal(w, k == 0 || k == 3 ? -1.0 : 1.0, k < 2 ? -1.0 : 1.0);
  }
  const auto position = [&](Scalar r, Scalar s, Scalar z) {
    const Eigen::Vector4d n(0.25 * (1 - r) * (1 - s), 0.25 * (1 + r) * (1 - s),
                            0.25 * (1 + r) * (1 + s), 0.25 * (1 - r) * (1 + s));
    return Vector3((w + (0.5 * z * t) * v) * n);
  };
  const std::array<Scalar, 6> gx{-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
                                 0.2386191860831969, 0.6612093864662645, 0.9324695142031521};
  const std::array<Scalar, 6> gw{0.1713244923791704, 0.3607615730481386, 0.4679139345726910,
                                 0.4679139345726910, 0.3607615730481386, 0.1713244923791704};
  const Scalar step = 1e-6;
  Scalar volume = 0.0;
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      for (int k = 0; k < 6; ++k) {
        const Scalar r = gx[static_cast<std::size_t>(i)];
        const Scalar s = gx[static_cast<std::size_t>(j)];
        const Scalar z = gx[static_cast<std::size_t>(k)];
        Matrix3 jac;
        jac.col(0) = (position(r + step, s, z) - position(r - step, s, z)) / (2 * step);
        jac.col(1) = (position(r, s + step, z) - position(r, s - step, z)) / (2 * step);
        jac.col(2) = (position(r, s, z + step) - position(r, s, z - step)) / (2 * step);
        volume += gw[static_cast<std::size_t>(i)] * gw[static_cast<std::size_t>(j)] *
                  gw[static_cast<std::size_t>(k)] * jac.determinant();
      }
    }
  }
  for (int m = 0; m < 3; ++m) {
    REQUIRE(rw.col(m).dot(mw * rw.col(m)) == Approx(rho * volume).epsilon(1e-8));
  }
  // It differs from rho t A by the tilt of the directors (here 0.7 %).
  REQUIRE(std::abs(volume / (t * Shell4Element::area(w)) - 1.0) < 0.01);
}

TEST_CASE("shell loads: a pressure against the normal, an edge traction over its area",
          "[shell][element]") {
  const IntegrationOptions opts;
  const Shell4Element element;
  const Matrix w = warped_cell();
  const Scalar p = 3.0e5;
  const Vector f = element.pressure_load(w, p, opts);
  Vector3 total = Vector3::Zero();
  for (int k = 0; k < 4; ++k) {
    total += f.segment<3>(6 * k);
    REQUIRE(f.segment<3>(6 * k + 3).norm() == 0.0);
  }
  // The resultant of a uniform pressure on a warped surface is -p times its
  // vector area, 1/2 (x3 - x1) x (x4 - x2).
  const Vector3 area_vector = 0.5 * Vector3(w.col(2) - w.col(0)).cross(Vector3(w.col(3) - w.col(1)));
  REQUIRE((total + p * area_vector).norm() <= 1e-12 * p * area_vector.norm());

  const Vector3 traction(1.0e6, -2.0e6, 0.5e6);
  const Scalar t = 0.02;
  const Vector fe = element.boundary_traction(w, 1, traction, t, opts);
  const Scalar length = (Vector3(w.col(2)) - Vector3(w.col(1))).norm();
  REQUIRE((fe.segment<3>(6) - 0.5 * length * t * traction).norm() <= 1e-9 * traction.norm());
  REQUIRE((fe.segment<3>(12) - 0.5 * length * t * traction).norm() <= 1e-9 * traction.norm());
  REQUIRE(fe.segment<3>(0).norm() == 0.0);
  REQUIRE(fe.segment<3>(18).norm() == 0.0);
}

TEST_CASE("the shell geometric stiffness: the in-plane stress on the slope of w, and its "
          "derivative",
          "[shell][element][buckling]") {
  const Matrix3 d = default_material(0.3).plane_stress_matrix();
  const IntegrationOptions opts;
  const Shell4Element element;
  const Scalar a = 1.0;
  const Scalar b = 0.5;
  const Scalar t = 0.01;
  // A uniform strain eps along x (free sideways): sigma_11 = E eps.
  const Scalar eps = -1.0e-4;
  Vector ue = Vector::Zero(24);
  const Matrix x = flat_cell(a, b);
  for (int k = 0; k < 4; ++k) {
    ue(6 * k) = eps * x(0, k);
    ue(6 * k + 1) = -0.3 * eps * x(1, k);
  }
  const Matrix kg = element.geometric_stiffness(x, d, ue, 1.0, t, opts);
  REQUIRE((kg - kg.transpose()).cwiseAbs().maxCoeff() <= 1e-14 * kg.cwiseAbs().maxCoeff());
  // w = x / a (a uniform slope along x): phi^T K_G phi = sigma_11 t A / a^2.
  Vector phi = Vector::Zero(24);
  for (int k = 0; k < 4; ++k) phi(6 * k + 2) = x(0, k) / a;
  const Scalar sigma = 70.0e9 * eps;
  REQUIRE(phi.dot(kg * phi) == Approx(sigma * t * a * b / (a * a)).epsilon(1e-12));

  // phi^T K_G(u) phi = g(phi)^T u for any u and phi (K_G is linear in u).
  std::mt19937 rng(7u);
  std::uniform_real_distribution<Scalar> dist(-1.0, 1.0);
  const Matrix w = with_tilted_directors(warped_cell());
  for (int trial = 0; trial < 3; ++trial) {
    Vector u(24);
    Vector m(24);
    for (int i = 0; i < 24; ++i) {
      u(i) = 1e-4 * dist(rng);
      m(i) = dist(rng);
    }
    const Scalar direct = m.dot(element.geometric_stiffness(w, d, u, 2.0, t, opts) * m);
    const Vector g = element.geometric_stiffness_derivative(w, d, m, 2.0, t, opts);
    REQUIRE(g.dot(u) == Approx(direct).epsilon(1e-11));
  }
}

namespace {

Matrix3 turned() {
  return (Eigen::AngleAxisd(0.7, Vector3(1.0, 2.0, -0.5).normalized()) *
          Eigen::AngleAxisd(-0.4, Vector3::UnitY()))
      .toRotationMatrix();
}

/// The distorted plate a x b of n x n cells, turned by `rot`.
Mesh turned_plate(const Matrix3& rot, Scalar a, Scalar b, Index n) {
  ShellMeshSpec spec;
  spec.n1 = n;
  spec.n2 = n;
  spec.lx = a;
  spec.ly = b;
  spec.perturbation = 0.2;
  spec.seed = 3u;
  const Mesh flat = make_structured_shell_mesh(spec);
  Mesh mesh(rot * flat.coordinates(), flat.connectivity(), ElementType::Shell4);
  mesh.set_node_normals(rot * flat.node_normals());
  return mesh;
}

/// Prescribe every DOF of `nodes` from the nodal field `u` (6 per node).
std::vector<DisplacementConstraint> prescribe(const std::vector<Index>& nodes, const Vector& u) {
  std::vector<DisplacementConstraint> out;
  for (Index n : nodes) {
    DisplacementConstraint c;
    c.region.name = "node" + std::to_string(n);
    c.region.members.emplace_back();
    c.region.members.back().kind = SelectorKind::NodeIds;
    c.region.members.back().ids.push_back(n);
    for (int k = 0; k < 6; ++k) c.set(k, true, u(6 * n + k));
    out.push_back(c);
  }
  return out;
}

Vector solve_prescribed(FemModel& model) {
  LoadCaseSpec lc;
  lc.name = "prescribed";
  lc.prescribed_displacement_only = true;
  model.load_case_specs() = {lc};
  model.finalize();
  Assembler assembler(model);
  StaticAnalysisOptions options;
  options.linear.type = LinearSolverType::SimplicialLdlt;
  StaticAnalysis analysis(model, assembler, options);
  return analysis.solve_all().front().displacement;
}

/// The in-plane tensor with components `c` = (11, 22, 12) on the basis
/// (a, b), as a 3 x 3 tensor.
Matrix3 tensor(const Vector3& c, const Vector3& a, const Vector3& b) {
  return c(0) * a * a.transpose() + c(1) * b * b.transpose() +
         c(2) * (a * b.transpose() + b * a.transpose());
}

}  // namespace

TEST_CASE("shell patch tests: constant membrane and bending states on a distorted mesh in "
          "a turned plane are exact",
          "[shell][model][patch]") {
  const Scalar a = 1.0;
  const Scalar b = 0.8;
  const Scalar t = 0.02;
  const Scalar e_mod = 70.0e9;
  const Scalar nu = 0.3;
  const Matrix3 rot = turned();
  const Vector3 e1 = rot.col(0);
  const Vector3 e2 = rot.col(1);
  const Vector3 n = rot.col(2);
  Mesh mesh = turned_plate(rot, a, b, 4);
  // u = A x + c in the plane, w quadratic: constant strains and curvatures.
  Eigen::Matrix2d amat;
  amat << 2.0e-4, -3.0e-4, 1.0e-4, -1.5e-4;
  const Scalar kx = 2.0e-3;
  const Scalar ky = -1.0e-3;
  const Scalar kxy = 1.5e-3;
  const Scalar omega = 0.5 * (amat(1, 0) - amat(0, 1));
  const Index nn = mesh.num_nodes();
  Vector exact(6 * nn);
  for (Index k = 0; k < nn; ++k) {
    const Vector3 local = rot.transpose() * mesh.node(k);
    const Scalar x = local(0);
    const Scalar y = local(1);
    const Eigen::Vector2d u = amat * Eigen::Vector2d(x, y) + Eigen::Vector2d(1.0e-4, -2.0e-4);
    const Scalar w = 0.5 * kx * x * x + 0.5 * ky * y * y + kxy * x * y + 3.0e-4 - 2.0e-4 * x;
    const Scalar wx = kx * x + kxy * y - 2.0e-4;
    const Scalar wy = ky * y + kxy * x;
    exact.segment<3>(6 * k) = rot * Vector3(u(0), u(1), w);
    exact.segment<3>(6 * k + 3) = rot * Vector3(wy, -wx, omega);
  }
  std::vector<Index> boundary;
  std::vector<Index> interior;
  for (Index k = 0; k < nn; ++k) {
    const Vector3 local = rot.transpose() * mesh.node(k);
    const bool edge = std::abs(local(0)) < 1e-12 || std::abs(local(0) - a) < 1e-12 ||
                      std::abs(local(1)) < 1e-12 || std::abs(local(1) - b) < 1e-12;
    (edge ? boundary : interior).push_back(k);
  }
  REQUIRE(interior.size() == 9);
  FemModel model(std::move(mesh), IsotropicMaterial(e_mod, nu, 2700.0), t, StressState::Shell,
                 IntegrationOptions());
  model.constraints() = prescribe(boundary, exact);
  const Vector u = solve_prescribed(model);
  const Scalar scale = exact.cwiseAbs().maxCoeff();
  for (Index k : interior) {
    INFO("node " << k);
    REQUIRE((u.segment<6>(6 * k) - exact.segment<6>(6 * k)).cwiseAbs().maxCoeff() <=
            1e-10 * scale);
  }
  // The resultants: N = t D eps, M = -D_b (kappa), Q = 0 in every element.
  const Scalar d_b = e_mod * t * t * t / (12.0 * (1.0 - nu * nu));
  const Scalar c = e_mod * t / (1.0 - nu * nu);
  const Vector3 eps(amat(0, 0), amat(1, 1), 0.5 * (amat(0, 1) + amat(1, 0)));
  const Matrix3 n_exact = tensor(
      Vector3(c * (eps(0) + nu * eps(1)), c * (eps(1) + nu * eps(0)), c * (1.0 - nu) * eps(2)),
      e1, e2);
  const Matrix3 m_exact =
      tensor(Vector3(-d_b * (kx + nu * ky), -d_b * (ky + nu * kx), -d_b * (1.0 - nu) * kxy), e1,
             e2);
  Assembler assembler(model);
  const ShellField field = recover_shell_resultants(model, assembler, u);
  for (Index e = 0; e < model.mesh().num_elements(); ++e) {
    const ShellResultants& r = field.element[static_cast<std::size_t>(e)];
    INFO("element " << e);
    REQUIRE(std::abs(r.e3.dot(n)) == Approx(1.0).epsilon(1e-12));
    const Matrix3 n_h = tensor(r.membrane, r.e1, r.e2);
    const Matrix3 m_h = tensor(r.moment, r.e1, r.e2);
    REQUIRE((n_h - n_exact).cwiseAbs().maxCoeff() <= 1e-9 * n_exact.cwiseAbs().maxCoeff());
    REQUIRE((m_h - m_exact).cwiseAbs().maxCoeff() <= 1e-9 * m_exact.cwiseAbs().maxCoeff());
    REQUIRE(r.shear.norm() <= 1e-9 * m_exact.cwiseAbs().maxCoeff() / t);
  }
}

TEST_CASE("a curved shell panel moved rigidly follows the motion and carries no stress",
          "[shell][model]") {
  ShellMeshSpec spec;
  spec.shape = ShellShape::Cylinder;
  spec.axis = 1;
  spec.radius = 2.0;
  spec.length = 1.5;
  spec.angle_start = -20.0;
  spec.angle_end = 55.0;
  spec.n1 = 5;
  spec.n2 = 4;
  spec.origin = Vector3(0.3, -0.2, 0.1);
  Mesh mesh = make_structured_shell_mesh(spec);
  const Vector3 omega(0.3e-3, -0.5e-3, 0.8e-3);
  const Vector3 shift(1.0e-3, 2.0e-3, -1.0e-3);
  const Index nn = mesh.num_nodes();
  Vector exact(6 * nn);
  std::vector<Index> boundary;
  std::vector<Index> interior;
  for (Index k = 0; k < nn; ++k) {
    exact.segment<3>(6 * k) = shift + omega.cross(mesh.node(k));
    exact.segment<3>(6 * k + 3) = omega;
    const Index i = k % (spec.n1 + 1);
    const Index j = k / (spec.n1 + 1);
    (i == 0 || i == spec.n1 || j == 0 || j == spec.n2 ? boundary : interior).push_back(k);
  }
  FemModel model(std::move(mesh), default_material(), 0.05, StressState::Shell,
                 IntegrationOptions());
  model.constraints() = prescribe(boundary, exact);
  const Vector u = solve_prescribed(model);
  for (Index k : interior) {
    INFO("node " << k);
    REQUIRE((u.segment<6>(6 * k) - exact.segment<6>(6 * k)).cwiseAbs().maxCoeff() <=
            1e-12);
  }
  Assembler assembler(model);
  const ShellField field = recover_shell_resultants(model, assembler, u);
  // Against the membrane force of a strain as large as the motion, E t |omega|.
  const Scalar unit = 70.0e9 * 0.05 * omega.norm();
  for (const ShellResultants& r : field.element) {
    REQUIRE(r.membrane.norm() <= 1e-9 * unit);
    REQUIRE(r.moment.norm() <= 1e-9 * unit * 0.05);
    REQUIRE(r.shear.norm() <= 1e-9 * unit);
  }
  REQUIRE(field.element_strain_energy.maxCoeff() <= 1e-20 * unit * unit);
}

TEST_CASE("shell directors: exact normals, averages over a smooth surface, own normals at "
          "a fold",
          "[shell][model]") {
  // A cylinder panel without its normals: at an interior node the two cells'
  // normals average to the radial direction; at a free edge the one cell
  // keeps its own, half a cell's angle off.
  ShellMeshSpec spec;
  spec.shape = ShellShape::Cylinder;
  spec.radius = 1.0;
  spec.length = 1.0;
  spec.angle_start = 0.0;
  spec.angle_end = 60.0;
  spec.n1 = 4;
  spec.n2 = 2;
  const Mesh exact = make_structured_shell_mesh(spec);
  const Scalar half_cell = 7.5 * 3.14159265358979323846 / 180.0;
  {
    FemModel model(Mesh(exact.coordinates(), exact.connectivity(), ElementType::Shell4),
                   default_material(), 0.01, StressState::Shell, IntegrationOptions());
    for (Index e = 0; e < model.mesh().num_elements(); ++e) {
      const Matrix g = model.element_geometry(e);
      REQUIRE(g.rows() == 6);
      const Index* nodes = model.mesh().element_nodes(e);
      for (int k = 0; k < 4; ++k) {
        const Vector3 v = g.block<3, 1>(3, k);
        const Vector3 radial = exact.node_normals().col(nodes[k]);
        const Index i = nodes[k] % (spec.n1 + 1);
        if (i == 0 || i == spec.n1) {
          REQUIRE(std::acos(std::min(1.0, v.dot(radial))) == Approx(half_cell).epsilon(1e-12));
        } else {
          REQUIRE((v - radial).norm() <= 1e-14);
        }
      }
    }
  }
  {
    FemModel model(Mesh(exact), default_material(), 0.01, StressState::Shell,
                   IntegrationOptions());
    for (Index e = 0; e < model.mesh().num_elements(); ++e) {
      const Matrix g = model.element_geometry(e);
      const Index* nodes = model.mesh().element_nodes(e);
      for (int k = 0; k < 4; ++k) {
        REQUIRE((Vector3(g.block<3, 1>(3, k)) - Vector3(exact.node_normals().col(nodes[k])))
                    .norm() <= 1e-15);
      }
    }
  }
  // Cells 15 degrees apart are one surface; 25 degrees apart, a fold, unless
  // the fold angle is raised.
  spec.n1 = 2;
  spec.angle_end = 50.0;
  const Mesh coarse = make_structured_shell_mesh(spec);
  const Mesh bare(coarse.coordinates(), coarse.connectivity(), ElementType::Shell4);
  FemModel folded(bare, default_material(), 0.01, StressState::Shell, IntegrationOptions());
  const Index shared = 1;  // the middle node of the first row, in cells 0 and 1
  const auto director = [&](const FemModel& model, Index e, Index node) {
    const Index* nodes = model.mesh().element_nodes(e);
    for (int k = 0; k < 4; ++k) {
      if (nodes[k] == node) return Vector3(model.element_geometry(e).block<3, 1>(3, k));
    }
    return Vector3(Vector3::Zero());
  };
  REQUIRE((director(folded, 0, shared) - Shell4Element::normal(bare.element_coordinates(0), 1.0, -1.0))
              .norm() <= 1e-14);
  REQUIRE(director(folded, 0, shared).dot(director(folded, 1, shared)) ==
          Approx(std::cos(25.0 * 3.14159265358979323846 / 180.0)).epsilon(1e-12));
  FemModel smooth(bare, default_material(), 0.01, StressState::Shell, IntegrationOptions());
  ShellOptions options;
  options.fold_angle_deg = 30.0;
  smooth.set_shell_options(options);
  REQUIRE((director(smooth, 0, shared) - director(smooth, 1, shared)).norm() <= 1e-15);
  REQUIRE((director(smooth, 0, shared) - Vector3(coarse.node_normals().col(shared))).norm() <=
          1e-14);
  options.fold_angle_deg = 95.0;
  REQUIRE_THROWS_AS(smooth.set_shell_options(options), ConfigError);

  // Normals the mesh does not follow - 80 degrees off a flat plate's at its
  // centre node - would tilt the fibres nearly flat: refused; 70 degrees off
  // is taken as it stands.
  ShellMeshSpec plate;
  plate.n1 = 2;
  plate.n2 = 2;
  const Mesh flat = make_structured_shell_mesh(plate);
  const auto tilted = [&](Scalar degrees) {
    Mesh mesh(flat.coordinates(), flat.connectivity(), ElementType::Shell4);
    Matrix normals = flat.node_normals();
    const Scalar angle = degrees * 3.14159265358979323846 / 180.0;
    normals.col(4) = Vector3(std::sin(angle), 0.0, std::cos(angle));
    mesh.set_node_normals(normals);
    return FemModel(mesh, default_material(), 0.01, StressState::Shell, IntegrationOptions());
  };
  REQUIRE_THROWS_WITH(tilted(80.0), ContainsSubstring("at most 75 allowed"));
  const FemModel steep = tilted(70.0);
  REQUIRE(Vector3(steep.element_geometry(0).block<3, 1>(3, 2)).z() ==
          Approx(std::cos(70.0 * 3.14159265358979323846 / 180.0)).epsilon(1e-14));
}

TEST_CASE("shell model loads: a pressure over a panel, gravity through the volume",
          "[shell][model][loads]") {
  ShellMeshSpec spec;
  spec.shape = ShellShape::Sphere;
  spec.radius = 3.0;
  spec.angle_start = 10.0;
  spec.angle_end = 80.0;
  spec.polar_start = 40.0;
  spec.polar_end = 110.0;
  spec.n1 = 5;
  spec.n2 = 4;
  Mesh mesh = make_structured_shell_mesh(spec);
  Vector3 area_vector = Vector3::Zero();
  Scalar area = 0.0;
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    const Matrix x = mesh.element_coordinates(e);
    area_vector += 0.5 * Vector3(x.col(2) - x.col(0)).cross(Vector3(x.col(3) - x.col(1)));
    area += mesh.element_measure(e);
  }
  const Scalar t = 0.04;
  FemModel model(std::move(mesh), IsotropicMaterial(70.0e9, 0.3, 2700.0), t, StressState::Shell,
                 IntegrationOptions());
  DisplacementConstraint held;
  held.region.name = "equator";
  Selector s;
  s.kind = SelectorKind::Box;
  s.zmax = -1.0;
  held.region.members.push_back(s);
  for (int k = 0; k < 6; ++k) held.set(k, true);
  model.constraints() = {held};
  LoadCaseSpec lc;
  lc.name = "loads";
  PressureLoadSpec p;
  p.region.name = "all";
  p.region.members.push_back(Selector());
  p.pressure = 2.0e5;
  lc.pressures.push_back(p);
  model.load_case_specs() = {lc};
  LoadCaseSpec weight;
  weight.name = "weight";
  weight.gravity = Vector3(1.0, -2.0, -9.81);
  model.load_case_specs().push_back(weight);
  model.finalize();
  const auto resultant = [&](const Vector& f) {
    Vector3 r = Vector3::Zero();
    for (Index k = 0; k < model.mesh().num_nodes(); ++k) r += f.segment<3>(6 * k);
    return r;
  };
  // A pressure presses against the outward normals: its resultant is
  // -p times the vector area of the (bilinear) surface.
  const Vector3 fp = resultant(model.load_vectors()[0]);
  REQUIRE((fp + p.pressure * area_vector).norm() <= 1e-12 * p.pressure * area);
  // Self-weight is rho g times the volume the shell's mass matrix integrates:
  // that of the solid its directors sweep, which on facets of a curved
  // surface falls short of t A by the chords' tilt against the normals - a
  // discretisation error of order h^2.
  Assembler assembler(model);
  const auto swept_mass = [](const FemModel& m, const Assembler& a) {
    const SparseMatrix mass = a.assemble_mass(MassType::Consistent);
    Vector ones = Vector::Zero(m.dofs().num_dofs());
    for (Index k = 0; k < m.mesh().num_nodes(); ++k) ones(6 * k) = 1.0;
    return ones.dot(mass * ones);
  };
  const Scalar mass = swept_mass(model, assembler);
  const Vector3 fg = resultant(model.load_vectors()[1]);
  REQUIRE((fg - mass * weight.gravity).norm() <= 1e-12 * mass * weight.gravity.norm());
  REQUIRE(assembler.total_mass() == Approx(2700.0 * t * area).epsilon(1e-12));
  const Scalar coarse_gap = 1.0 - mass / assembler.total_mass();
  REQUIRE(coarse_gap > 0.0);
  REQUIRE(coarse_gap < 0.025);
  spec.n1 *= 2;
  spec.n2 *= 2;
  FemModel fine(make_structured_shell_mesh(spec), IsotropicMaterial(70.0e9, 0.3, 2700.0), t,
                StressState::Shell, IntegrationOptions());
  const Assembler fine_assembler(fine);
  const Scalar fine_gap = 1.0 - swept_mass(fine, fine_assembler) / fine_assembler.total_mass();
  INFO("t A - volume: " << coarse_gap << " then " << fine_gap << " of t A");
  REQUIRE(fine_gap > 0.0);
  REQUIRE(fine_gap < coarse_gap / 3.5);
}

TEST_CASE("a simply supported plate under pressure converges to Mindlin's and Kirchhoff's "
          "deflection without shear locking",
          "[shell][model][plate]") {
  const Scalar a = 1.0;
  const Scalar e_mod = 1.0e9;
  const Scalar nu = 0.3;
  const Scalar q = 1.0e3;
  // Navier's series for the centre deflection, w D / (q a^4) = 0.00406235;
  // Mindlin's adds the Marcus moment over the shear stiffness,
  // 0.0736713 q a^2 / (k G t) (hard simple support).
  const Scalar navier = 0.00406235;
  for (const Scalar t : {0.1, 0.001}) {
    ShellMeshSpec spec;
    spec.n1 = 16;
    spec.n2 = 16;
    spec.lx = a;
    spec.ly = a;
    FemModel model(make_structured_shell_mesh(spec), IsotropicMaterial(e_mod, nu, 0.0), t,
                   StressState::Shell, IntegrationOptions());
    DisplacementConstraint x_edges;
    x_edges.region.name = "x_edges";
    for (const Scalar x : {0.0, a}) {
      Selector s;
      s.kind = SelectorKind::Box;
      s.xmin = x;
      s.xmax = x;
      x_edges.region.members.push_back(s);
    }
    DisplacementConstraint y_edges;
    y_edges.region.name = "y_edges";
    for (const Scalar y : {0.0, a}) {
      Selector s;
      s.kind = SelectorKind::Box;
      s.ymin = y;
      s.ymax = y;
      y_edges.region.members.push_back(s);
    }
    for (int k : {0, 1, 2}) {
      x_edges.set(k, true);
      y_edges.set(k, true);
    }
    x_edges.set(3, true);  // the rotation along the edge: a hard support
    y_edges.set(4, true);
    model.constraints() = {x_edges, y_edges};
    LoadCaseSpec lc;
    lc.name = "pressure";
    PressureLoadSpec p;
    p.region.name = "all";
    p.region.members.push_back(Selector());
    p.pressure = q;
    lc.pressures.push_back(p);
    model.load_case_specs() = {lc};
    model.finalize();
    Assembler assembler(model);
    StaticAnalysis analysis(model, assembler);
    const StaticSolution sol = analysis.solve_all().front();
    const Index centre = 8 * 17 + 8;
    REQUIRE((model.mesh().node(centre) - Vector3(0.5 * a, 0.5 * a, 0.0)).norm() <= 1e-15);
    const Scalar d_b = e_mod * t * t * t / (12.0 * (1.0 - nu * nu));
    const Scalar g = e_mod / (2.0 * (1.0 + nu));
    const Scalar kirchhoff = navier * q * a * a * a * a / d_b;
    const Scalar mindlin = kirchhoff + 0.0736713 * q * a * a / (5.0 / 6.0 * g * t);
    const Scalar w = -sol.displacement(6 * centre + 2);
    INFO("t / a = " << t << ", w / w_Mindlin = " << w / mindlin);
    REQUIRE(w / mindlin == Approx(1.0).margin(2.0e-3));
    REQUIRE(w < mindlin);
  }
}

TEST_CASE("shell decks: a generated cylinder with sections, rotations and moments, and the "
          "analyses a shell refuses",
          "[shell][io][config]") {
  const std::string deck = R"({
    "name": "panel",
    "mesh": { "type": "structured_shell", "shape": "cylinder", "axis": "y", "radius": 2.0,
              "length": 1.0, "angles": [0.0, 60.0], "n_around": 6, "n_along": 4 },
    "material": { "youngs_modulus": 70e9, "poisson_ratio": 0.3, "density": 2700.0 },
    "model": { "thickness": 0.01,
               "shell": { "drilling_stiffness": 1e-4, "fold_angle": 15.0,
                          "sections": [ { "name": "thick_band", "thickness": 0.02,
                                          "region": { "box": { "ymax": 0.5 } } } ] } },
    "boundary_conditions": [
      { "name": "root", "fix": ["x", "y", "z", "rx", "ry", "rz"],
        "region": { "box": { "ymax": 0.0 } } },
      { "name": "turned", "fix": ["rx"], "rotation": [0.001, 0.0, 0.0],
        "region": { "box": { "ymin": 1.0 } } }
    ],
    "load_cases": [ { "name": "tip", "point_loads": [
      { "force": [0.0, 0.0, 10.0], "moment": [1.0, 2.0, 3.0], "distribution": "per_node",
        "region": { "box": { "ymin": 1.0, "zmin": 1.99 } } } ] } ]
  })";
  const Configuration config = parse_configuration(json::parse(deck, "panel.json"), "panel", true);
  REQUIRE(config.is_shell());
  REQUIRE(config.dim() == 3);
  REQUIRE(config.shell_mesh.shape == ShellShape::Cylinder);
  REQUIRE(config.shell_mesh.axis == 1);
  REQUIRE(config.describe_mesh() == "structured_shell cylinder 6 x 4");
  REQUIRE(config.constraints[1].fix_rx);
  REQUIRE(config.constraints[1].value_rx == 0.001);
  REQUIRE(config.load_cases[0].point_loads[0].moment == Vector3(1.0, 2.0, 3.0));
  const FemModel model = build_model(config);
  REQUIRE(model.mesh().element_type() == ElementType::Shell4);
  REQUIRE(model.dofs_per_node() == 6);
  REQUIRE(model.shell_options().drilling_factor == 1e-4);
  REQUIRE(model.shell_options().fold_angle_deg == 15.0);
  Index thick = 0;
  for (Index e = 0; e < model.mesh().num_elements(); ++e) {
    const bool in_band = model.mesh().element_centroid(e).y() <= 0.5;
    REQUIRE(model.thickness_of(e) == (in_band ? 0.02 : 0.01));
    thick += in_band ? 1 : 0;
  }
  REQUIRE(thick == 12);
  // The moment lands on the rotations of the loaded node, the force on its
  // translations.
  const Vector& f = model.load_vectors()[0];
  Index loaded = -1;
  for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
    if (f.segment<6>(6 * n).norm() > 0.0) loaded = n;
  }
  REQUIRE(loaded >= 0);
  REQUIRE(Vector(f.segment<6>(6 * loaded)).isApprox(
      (Vector(6) << 0.0, 0.0, 10.0, 1.0, 2.0, 3.0).finished()));

  // What a shell does not take is refused with the reason.
  const auto refused = [&](const std::string& extra) {
    std::string text = deck;
    text.insert(text.rfind('}'), "," + extra);
    try {
      parse_configuration(json::parse(text, "panel.json"), "panel", true);
    } catch (const ConfigError& e) {
      return std::string(e.what());
    }
    return std::string();
  };
  REQUIRE_THAT(refused(R"("nonlinear": { "enabled": true })"),
               ContainsSubstring("small rotations"));
  REQUIRE_THAT(refused(R"("topology": { "enabled": true })"),
               ContainsSubstring("thickness"));
  std::string with_temperature = deck;
  with_temperature.insert(with_temperature.find("\"point_loads\""),
                          R"("temperature": { "uniform": 20.0 }, )");
  REQUIRE_THROWS_WITH(parse_configuration(json::parse(with_temperature, "t.json"), "t", true),
                      ContainsSubstring("thermal strain"));
  std::string no_thickness = deck;
  no_thickness.replace(no_thickness.find("\"thickness\": 0.01,"),
                       std::string("\"thickness\": 0.01,").size(), "");
  REQUIRE_THROWS_WITH(parse_configuration(json::parse(no_thickness, "t.json"), "t", true),
                      ContainsSubstring("'model.thickness' is required"));
  // Rotations and moments belong to shells.
  const std::string solid = R"({
    "mesh": { "type": "structured_hex", "nx": 2, "ny": 1, "nz": 1, "lx": 1, "ly": 1, "lz": 1 },
    "material": { "youngs_modulus": 1e9, "poisson_ratio": 0.3 },
    "boundary_conditions": [ { "fix": ["x", "rx"], "region": { "box": { "xmax": 0 } } } ],
    "load_cases": [ { "point_loads": [ { "force": [1, 0, 0], "region": { "box": { "xmin": 1 } } } ] } ]
  })";
  REQUIRE_THROWS_WITH(parse_configuration(json::parse(solid, "s.json"), "s", true),
                      ContainsSubstring("carry translations only"));
}

TEST_CASE("shell meshes from files: S4R cells in space, their node order kept; other shells "
          "refused",
          "[shell][io][unstructured]") {
  // Two S4R cells folded along x = 1 - one in z = 0, one rising in z -
  // the second listed the other way round (its normal points down).
  const std::string inp =
      "*NODE\n"
      "1, 0, 0, 0\n2, 1, 0, 0\n3, 1, 1, 0\n4, 0, 1, 0\n5, 1, 0, 1\n6, 1, 1, 1\n"
      "*ELEMENT, TYPE=S4R, ELSET=SKIN\n"
      "11, 1, 2, 3, 4\n12, 2, 5, 6, 3\n"
      "*SHELL SECTION, ELSET=SKIN, MATERIAL=AL\n0.002\n";
  std::istringstream in(inp);
  MeshReadReport report;
  const Mesh mesh = read_abaqus_inp(in, "skin.inp", ".", MeshReadOptions(), &report);
  REQUIRE(mesh.element_type() == ElementType::Shell4);
  REQUIRE(mesh.dim() == 3);
  REQUIRE(report.dimension == 3);
  REQUIRE(report.reoriented == 0);
  REQUIRE(mesh.num_elements() == 2);
  REQUIRE(mesh.element_sets().at("SKIN") == std::vector<Index>({0, 1}));
  REQUIRE(mesh.node(4) == Vector3(1.0, 0.0, 1.0));
  REQUIRE(Shell4Element::normal(mesh.element_coordinates(0), 0.0, 0.0).isApprox(Vector3::UnitZ()));
  REQUIRE(Shell4Element::normal(mesh.element_coordinates(1), 0.0, 0.0).isApprox(-Vector3::UnitX()));
  REQUIRE(mesh.element_measure(1) == Approx(1.0));
  // At the fold both cells keep their own normals as directors.
  const FemModel model(mesh, default_material(), 0.002, StressState::Shell, IntegrationOptions());
  REQUIRE(Vector3(model.element_geometry(0).block<3, 1>(3, 1)).isApprox(Vector3::UnitZ()));
  REQUIRE(Vector3(model.element_geometry(1).block<3, 1>(3, 0)).isApprox(-Vector3::UnitX()));

  const auto read = [](const std::string& text) {
    std::istringstream stream(text);
    return read_abaqus_inp(stream, "bad.inp");
  };
  REQUIRE_THROWS_WITH(read("*NODE\n1, 0, 0, 0\n2, 1, 0, 0\n3, 0, 1, 0\n"
                           "*ELEMENT, TYPE=S3\n1, 1, 2, 3\n"),
                      ContainsSubstring("MITC4"));
  REQUIRE_THROWS_WITH(read("*NODE\n1, 0, 0, 0\n2, 1, 0, 0\n3, 1, 1, 0\n4, 0, 1, 0\n"
                           "5, 0.5, 0, 0\n6, 1, 0.5, 0\n7, 0.5, 1, 0\n8, 0, 0.5, 0\n"
                           "*ELEMENT, TYPE=S8R\n1, 1, 2, 3, 4, 5, 6, 7, 8\n"),
                      ContainsSubstring("MITC4"));
  // A plane quadrilateral file read as a shell.
  MeshReadOptions as_shell;
  as_shell.shell = true;
  std::istringstream plane("*NODE\n1, 0, 0\n2, 2, 0\n3, 2, 1\n4, 0, 1\n"
                           "*ELEMENT, TYPE=CPS4\n1, 1, 2, 3, 4\n");
  const Mesh flat = read_abaqus_inp(plane, "flat.inp", ".", as_shell);
  REQUIRE(flat.element_type() == ElementType::Shell4);
  REQUIRE(flat.dim() == 3);
  REQUIRE(flat.element_measure(0) == Approx(2.0));
  std::istringstream triangles("*NODE\n1, 0, 0\n2, 2, 0\n3, 2, 1\n*ELEMENT, TYPE=CPS3\n1, 1, 2, 3\n");
  REQUIRE_THROWS_WITH(read_abaqus_inp(triangles, "tri.inp", ".", as_shell),
                      ContainsSubstring("read as a shell"));
}

TEST_CASE("CalculiX shell decks: S4 cells, a section per thickness, the pressure's sign and "
          "results at the shell's own nodes",
          "[shell][io][cross-validation]") {
  ensure_directory("results/_test_tmp");
  ShellMeshSpec spec;
  spec.n1 = 4;
  spec.n2 = 2;
  FemModel model(make_structured_shell_mesh(spec), default_material(), 0.01, StressState::Shell,
                 IntegrationOptions());
  model.assign_thickness(0.02, {0, 1});
  DisplacementConstraint edge;
  edge.region.name = "edge";
  Selector s;
  s.kind = SelectorKind::Box;
  s.xmax = 0.0;
  edge.region.members.push_back(s);
  for (int k = 0; k < 6; ++k) edge.set(k, true);
  model.constraints() = {edge};
  LoadCaseSpec lc;
  lc.name = "pressure";
  PressureLoadSpec p;
  p.region.name = "all";
  p.region.members.push_back(Selector());
  p.pressure = 250.0;
  lc.pressures.push_back(p);
  PointLoadSpec m;
  m.region.name = "corner";
  Selector c;
  c.kind = SelectorKind::Box;
  c.xmin = 1.0;
  c.ymin = 1.0;
  m.region.members.push_back(c);
  m.moment = Vector3(0.0, 3.0, 0.0);
  lc.point_loads.push_back(m);
  model.load_case_specs() = {lc};
  model.finalize();
  REQUIRE(calculix_element_type(model) == "S4");
  const std::vector<std::string> decks = write_calculix_decks(model, "results/_test_tmp/shell", "unit");
  REQUIRE(decks.size() == 1);
  std::ifstream in(decks.front());
  std::stringstream buffer;
  buffer << in.rdbuf();
  const std::string text = buffer.str();
  REQUIRE_THAT(text, ContainsSubstring("*ELEMENT, TYPE=S4, ELSET=EALL"));
  // Two sections: the thick cells 1 and 2, the rest.
  REQUIRE_THAT(text, ContainsSubstring("*ELSET, ELSET=S1\n1, 2\n"));
  REQUIRE_THAT(text, ContainsSubstring("*ELSET, ELSET=S2\n3, 4, 5, 6, 7, 8\n"));
  REQUIRE_THAT(text, ContainsSubstring("*SHELL SECTION, ELSET=S1, MATERIAL=MAT1\n0.02\n"));
  REQUIRE_THAT(text, ContainsSubstring("*SHELL SECTION, ELSET=S2, MATERIAL=MAT1\n0.01\n"));
  REQUIRE(text.find("*SOLID SECTION") == std::string::npos);
  // CalculiX's shell P acts along the element normal, SparLab's against it.
  REQUIRE_THAT(text, ContainsSubstring("*DLOAD\n1, P, -250\n"));
  // The moment about y on the corner node (15 in CalculiX's numbering).
  REQUIRE_THAT(text, ContainsSubstring("15, 5, 3\n"));
  // Rotations are held on the edge: DOFs 4 to 6.
  REQUIRE_THAT(text, ContainsSubstring("1, 6, 6, 0\n"));
  REQUIRE_THAT(text, ContainsSubstring("*NODE FILE, OUTPUT=2D\nU\n"));
  std::remove(decks.front().c_str());
}
