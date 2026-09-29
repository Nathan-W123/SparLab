// Unit tests for the two-node Timoshenko beam (Beam2) and its sections.

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/elements/Beam2.hpp"
#include "sparlab/elements/BeamSection.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/io/CalculixWriter.hpp"
#include "sparlab/io/Config.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/io/Json.hpp"
#include "sparlab/io/MeshReader.hpp"
#include "sparlab/io/ResultWriter.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <random>
#include <sstream>

using namespace sparlab;
using Catch::Approx;
using Catch::Matchers::ContainsSubstring;

namespace {

constexpr Scalar kPi = 3.14159265358979323846;
constexpr Scalar kE = 70.0e9;
constexpr Scalar kNu = 0.3;
const Scalar kG = kE / (2.0 * (1.0 + kNu));

/// A general section with distinct properties in every direction.
BeamSection general_section(Scalar ky = 0.8, Scalar kz = 0.7) {
  BeamSection s;
  s.name = "test";
  s.area = 3.0e-3;
  s.iy = 4.0e-6;
  s.iz = 1.5e-6;
  s.torsion = 2.0e-6;
  s.shear_y = ky;
  s.shear_z = kz;
  return resolve(s, kNu);
}

Matrix material() {
  Matrix d = Matrix::Zero(2, 2);
  d(0, 0) = kE;
  d(1, 1) = kG;
  return d;
}

/// The six rigid-body motions of the nodes x (3 x 2), 12 x 6.
Matrix rigid_modes(const Matrix& x) {
  Matrix r = Matrix::Zero(12, 6);
  for (int a = 0; a < 2; ++a) {
    const Vector3 p = x.col(a);
    for (int k = 0; k < 3; ++k) r(6 * a + k, k) = 1.0;
    // Rotation about axis k: u = e_k x p, theta = e_k.
    for (int k = 0; k < 3; ++k) {
      const Vector3 u = Vector3::Unit(k).cross(p);
      r.block<3, 1>(6 * a, 3 + k) = u;
      r(6 * a + 3 + k, 3 + k) = 1.0;
    }
  }
  return r;
}

/// Przemieniecki's stiffness of the Timoshenko beam along x in its own axes.
Matrix closed_form_stiffness(const BeamSection& s, Scalar length) {
  const Scalar l = length;
  const Scalar py = 12.0 * kE * s.iz / (s.shear_y * kG * s.area * l * l);
  const Scalar pz = 12.0 * kE * s.iy / (s.shear_z * kG * s.area * l * l);
  Matrix k = Matrix::Zero(12, 12);
  const Scalar ax = kE * s.area / l;
  const Scalar tq = kG * s.torsion / l;
  k(0, 0) = k(6, 6) = ax;
  k(0, 6) = k(6, 0) = -ax;
  k(3, 3) = k(9, 9) = tq;
  k(3, 9) = k(9, 3) = -tq;
  // x-y plane: (v, theta_z) at DOFs 1, 5, 7, 11.
  {
    const Scalar c = kE * s.iz / ((1.0 + py) * l * l * l);
    const int i[4] = {1, 5, 7, 11};
    const Scalar m[4][4] = {{12.0, 6.0 * l, -12.0, 6.0 * l},
                            {6.0 * l, (4.0 + py) * l * l, -6.0 * l, (2.0 - py) * l * l},
                            {-12.0, -6.0 * l, 12.0, -6.0 * l},
                            {6.0 * l, (2.0 - py) * l * l, -6.0 * l, (4.0 + py) * l * l}};
    for (int a = 0; a < 4; ++a) {
      for (int b = 0; b < 4; ++b) k(i[a], i[b]) = c * m[a][b];
    }
  }
  // x-z plane: (w, theta_y) at DOFs 2, 4, 8, 10 (theta_y = -dw/dx).
  {
    const Scalar c = kE * s.iy / ((1.0 + pz) * l * l * l);
    const int i[4] = {2, 4, 8, 10};
    const Scalar m[4][4] = {{12.0, -6.0 * l, -12.0, -6.0 * l},
                            {-6.0 * l, (4.0 + pz) * l * l, 6.0 * l, (2.0 - pz) * l * l},
                            {-12.0, 6.0 * l, 12.0, 6.0 * l},
                            {-6.0 * l, (2.0 - pz) * l * l, 6.0 * l, (4.0 + pz) * l * l}};
    for (int a = 0; a < 4; ++a) {
      for (int b = 0; b < 4; ++b) k(i[a], i[b]) = c * m[a][b];
    }
  }
  return k;
}

}  // namespace

TEST_CASE("the beam stiffness is Przemieniecki's closed form and has the six rigid-body "
          "motions as its null space",
          "[beam][element]") {
  const Beam2Element element;
  const IntegrationOptions opts;
  const BeamSection s = general_section();
  // Along x: the local and global axes coincide.
  for (const Scalar l : {0.05, 0.4, 3.0}) {
    INFO("L = " << l);
    const Matrix g = Beam2Element::geometry(Vector3::Zero(), Vector3(l, 0.0, 0.0), s, kE, kG);
    const Matrix k = element.stiffness(g, material(), 1.0, opts);
    const Matrix exact = closed_form_stiffness(s, l);
    REQUIRE((k - exact).cwiseAbs().maxCoeff() <= 1e-12 * exact.cwiseAbs().maxCoeff());
  }
  // Oblique, with an orientation vector: symmetric, the rigid-body motions free,
  // no other zero-energy mode.
  BeamSection oriented = s;
  oriented.orientation = Vector3(0.3, -1.0, 0.4);
  const Vector3 x0(0.2, -0.1, 0.5);
  const Vector3 x1(1.1, 0.6, -0.3);
  const Matrix g = Beam2Element::geometry(x0, x1, oriented, kE, kG);
  const Matrix k = element.stiffness(g, material(), 1.0, opts);
  REQUIRE((k - k.transpose()).cwiseAbs().maxCoeff() <= 1e-15 * k.cwiseAbs().maxCoeff());
  Matrix x(3, 2);
  x << x0, x1;
  const Matrix r = rigid_modes(x);
  REQUIRE((k * r).norm() <= 1e-12 * k.norm() * r.norm());
  Eigen::SelfAdjointEigenSolver<Matrix> eig(k);
  int zero = 0;
  for (Eigen::Index i = 0; i < 12; ++i) {
    REQUIRE(eig.eigenvalues()(i) >= -1e-12 * eig.eigenvalues().maxCoeff());
    zero += eig.eigenvalues()(i) < 1e-12 * eig.eigenvalues().maxCoeff();
  }
  REQUIRE(zero == 6);
}

TEST_CASE("without shear deformation the beam is the Euler-Bernoulli element",
          "[beam][element]") {
  BeamSection s;
  s.name = "eb";
  s.area = 3.0e-3;
  s.iy = 4.0e-6;
  s.iz = 1.5e-6;
  s.torsion = 2.0e-6;
  s.shear_deformation = false;
  const BeamSection eb = resolve(s, kNu);
  const Scalar l = 0.7;
  const Matrix g = Beam2Element::geometry(Vector3::Zero(), Vector3(l, 0.0, 0.0), eb, kE, kG);
  const Matrix k = Beam2Element().stiffness(g, material(), 1.0, IntegrationOptions());
  // Phi = 0: the Hermite stiffness 12 EI / L^3, 6 EI / L^2, 4 EI / L, 2 EI / L.
  REQUIRE(k(1, 1) == Approx(12.0 * kE * eb.iz / (l * l * l)).epsilon(1e-13));
  REQUIRE(k(1, 5) == Approx(6.0 * kE * eb.iz / (l * l)).epsilon(1e-13));
  REQUIRE(k(5, 5) == Approx(4.0 * kE * eb.iz / l).epsilon(1e-13));
  REQUIRE(k(5, 11) == Approx(2.0 * kE * eb.iz / l).epsilon(1e-13));
  REQUIRE(k(2, 2) == Approx(12.0 * kE * eb.iy / (l * l * l)).epsilon(1e-13));
  REQUIRE(k(4, 4) == Approx(4.0 * kE * eb.iy / l).epsilon(1e-13));
  // A section without shear deformation takes no shear coefficients.
  BeamSection bad = s;
  bad.shear_y = 0.8;
  REQUIRE_THROWS_WITH(resolve(bad, kNu), ContainsSubstring("takes no shear coefficients"));
}

TEST_CASE("the beam mass holds rho A L in each direction and the rotary inertia of the "
          "section",
          "[beam][element]") {
  const Beam2Element element;
  const BeamSection s = general_section();
  const Scalar rho = 2700.0;
  const Vector3 x0(0.2, -0.1, 0.5);
  const Vector3 x1(1.1, 0.6, -0.3);
  const Scalar l = (x1 - x0).norm();
  const Matrix g = Beam2Element::geometry(x0, x1, s, kE, kG);
  const Matrix m = element.consistent_mass(g, rho, 1.0, IntegrationOptions());
  REQUIRE((m - m.transpose()).cwiseAbs().maxCoeff() <= 1e-15 * m.cwiseAbs().maxCoeff());
  Matrix x(3, 2);
  x << x0, x1;
  const Matrix r = rigid_modes(x);
  for (int k = 0; k < 3; ++k) {
    REQUIRE(r.col(k).dot(m * r.col(k)) == Approx(rho * s.area * l).epsilon(1e-13));
  }
  // A rigid rotation about the local axes through node 0 carries the kinetic
  // energy of the rod about that point plus the section's rotary inertia.
  const BeamFrame f = Beam2Element::frame(g);
  const auto rotation_about = [&](const Vector3& axis) {
    Vector q = Vector::Zero(12);
    for (int a = 0; a < 2; ++a) {
      q.segment<3>(6 * a) = axis.cross(Vector3(x.col(a)) - x0);
      q.segment<3>(6 * a + 3) = axis;
    }
    return q;
  };
  const Vector twist = rotation_about(f.x_axis());
  REQUIRE(twist.dot(m * twist) == Approx(rho * s.polar() * l).epsilon(1e-13));
  const Vector about_z = rotation_about(f.z_axis());
  REQUIRE(about_z.dot(m * about_z) ==
          Approx(rho * (s.area * l * l * l / 3.0 + s.iz * l)).epsilon(1e-13));
  const Vector about_y = rotation_about(f.y_axis());
  REQUIRE(about_y.dot(m * about_y) ==
          Approx(rho * (s.area * l * l * l / 3.0 + s.iy * l)).epsilon(1e-13));
  // Positive definite: every rotation has inertia.
  Eigen::SelfAdjointEigenSolver<Matrix> eig(m);
  REQUIRE(eig.eigenvalues().minCoeff() > 0.0);
}

TEST_CASE("the beam's lumped mass: rho A L / 2 per node and half the sections' inertia "
          "tensor, turned with the element",
          "[beam][element]") {
  const Beam2Element element;
  BeamSection s = general_section();
  s.orientation = Vector3(0.3, 1.0, 0.2);
  const Scalar rho = 2700.0;
  const Vector3 x0(0.2, -0.1, 0.5);
  const Vector3 x1(1.1, 0.6, -0.3);
  const Scalar l = (x1 - x0).norm();
  const Matrix g = Beam2Element::geometry(x0, x1, s, kE, kG);
  const Matrix m = element.lumped_mass(g, rho, 1.0, IntegrationOptions());
  REQUIRE(m.rows() == 12);
  REQUIRE((m - m.transpose()).cwiseAbs().maxCoeff() == 0.0);
  const BeamFrame f = Beam2Element::frame(g);
  const Scalar half = 0.5 * rho * l;
  for (int a = 0; a < 2; ++a) {
    REQUIRE((m.block<3, 3>(6 * a, 6 * a) - half * s.area * Matrix3::Identity())
                .cwiseAbs()
                .maxCoeff() <= 1e-15 * half * s.area);
    // Nothing couples the nodes, or a node's translations to its rotations.
    REQUIRE(m.block<3, 3>(6 * a, 6 * a + 3).cwiseAbs().maxCoeff() == 0.0);
    REQUIRE(m.block<6, 6>(6 * a, 6 * (1 - a)).cwiseAbs().maxCoeff() == 0.0);
    // In the element's axes the rotational block is diag(I_p, I_y, I_z).
    const Matrix3 local = f.rotation * m.block<3, 3>(6 * a + 3, 6 * a + 3) *
                          f.rotation.transpose();
    const Matrix3 expected = half * Vector3(s.polar(), s.iy, s.iz).asDiagonal();
    REQUIRE((local - expected).cwiseAbs().maxCoeff() <= 1e-14 * half * s.polar());
  }
  // The inclined element's rotational block is not diagonal in global axes:
  // no diagonal could give the twist its inertia rho I_p L here.
  REQUIRE(std::abs(m(3, 4)) > 1e-3 * m(3, 3));
  Eigen::SelfAdjointEigenSolver<Matrix> eig(m);
  REQUIRE(eig.eigenvalues().minCoeff() > 0.0);
}

TEST_CASE("the beam geometric stiffness: the axial force on the slopes, the rotations' "
          "gradients, and its derivative",
          "[beam][element][buckling]") {
  const Beam2Element element;
  const IntegrationOptions opts;
  const BeamSection s = general_section();
  const Scalar l = 0.8;
  const Matrix g = Beam2Element::geometry(Vector3::Zero(), Vector3(l, 0.0, 0.0), s, kE, kG);
  // A compression of 1e-4: N = -E A 1e-4.
  Vector u = Vector::Zero(12);
  u(6) = -1.0e-4 * l;
  const Scalar n = -kE * s.area * 1.0e-4;
  const Matrix kg = element.geometric_stiffness(g, material(), u, 1.0, 1.0, opts);
  REQUIRE((kg - kg.transpose()).cwiseAbs().maxCoeff() <= 1e-14 * kg.cwiseAbs().maxCoeff());
  // A rigid rotation about z: v = x, theta_z = 1 -> N int v'^2 = N L.
  Vector phi = Vector::Zero(12);
  phi(7) = l;
  phi(5) = phi(11) = 1.0;
  REQUIRE(phi.dot(kg * phi) == Approx(n * l).epsilon(1e-12));
  // A uniform twist rate: N (I_p / A) L / L^2.
  phi.setZero();
  phi(9) = 1.0;
  REQUIRE(phi.dot(kg * phi) == Approx(n * s.polar() / (s.area * l)).epsilon(1e-12));
  // A uniform stretch: N int u'^2 = N / L.
  phi.setZero();
  phi(6) = 1.0;
  REQUIRE(phi.dot(kg * phi) == Approx(n / l).epsilon(1e-12));
  // phi^T K_G(u) phi = g(phi)^T u on an oblique element.
  BeamSection oriented = s;
  oriented.orientation = Vector3(0.0, 0.3, 1.0);
  const Matrix go =
      Beam2Element::geometry(Vector3(0.1, 0.2, 0.3), Vector3(0.9, -0.4, 0.8), oriented, kE, kG);
  std::mt19937 rng(11u);
  std::uniform_real_distribution<Scalar> dist(-1.0, 1.0);
  for (int trial = 0; trial < 3; ++trial) {
    Vector ur(12);
    Vector pr(12);
    for (int i = 0; i < 12; ++i) {
      ur(i) = 1e-4 * dist(rng);
      pr(i) = dist(rng);
    }
    const Scalar direct = pr.dot(element.geometric_stiffness(go, material(), ur, 2.0, 1.0, opts) * pr);
    const Vector gd = element.geometric_stiffness_derivative(go, material(), pr, 2.0, 1.0, opts);
    REQUIRE(gd.dot(ur) == Approx(direct).epsilon(1e-12));
  }
}

TEST_CASE("a uniform line load: q L / 2 and q L^2 / 12 at the ends, whatever the shear "
          "parameter",
          "[beam][element][loads]") {
  const Beam2Element element;
  for (const Scalar kz : {0.8, 0.01}) {
    const BeamSection s = general_section(0.8, kz);
    const Scalar l = 0.6;
    const Matrix g = Beam2Element::geometry(Vector3::Zero(), Vector3(l, 0.0, 0.0), s, kE, kG);
    const Scalar q = -1500.0;
    const Vector f = element.line_load(g, Vector3(0.0, 0.0, q));
    INFO("k_z = " << kz);
    REQUIRE(f(2) == Approx(0.5 * q * l).epsilon(1e-13));
    REQUIRE(f(8) == Approx(0.5 * q * l).epsilon(1e-13));
    // theta_y = -dw/dx: a downward load's fixed-end moments.
    REQUIRE(f(4) == Approx(-q * l * l / 12.0).epsilon(1e-13));
    REQUIRE(f(10) == Approx(q * l * l / 12.0).epsilon(1e-13));
    for (const int i : {0, 1, 3, 5, 6, 7, 9, 11}) REQUIRE(std::abs(f(i)) <= 1e-12 * std::abs(q));
  }
}

TEST_CASE("beam axes: the default frame puts z' up, an orientation vector fixes y'",
          "[beam][element]") {
  const BeamSection s = general_section();
  const auto frame_of = [&](const Vector3& x1, const Vector3& orientation) {
    BeamSection o = s;
    o.orientation = orientation;
    return Beam2Element::frame(Beam2Element::geometry(Vector3::Zero(), x1, o, kE, kG));
  };
  // Along x: y' = Y, z' = Z; along y: y' = -X, z' = Z.
  BeamFrame f = frame_of(Vector3(2.0, 0.0, 0.0), Vector3::Zero());
  REQUIRE((f.y_axis() - Vector3::UnitY()).norm() <= 1e-15);
  REQUIRE((f.z_axis() - Vector3::UnitZ()).norm() <= 1e-15);
  f = frame_of(Vector3(0.0, 2.0, 0.0), Vector3::Zero());
  REQUIRE((f.y_axis() + Vector3::UnitX()).norm() <= 1e-15);
  REQUIRE((f.z_axis() - Vector3::UnitZ()).norm() <= 1e-15);
  // Vertical: z' along global X, so y' = z' x x' = -Y.
  f = frame_of(Vector3(0.0, 0.0, 3.0), Vector3::Zero());
  REQUIRE((f.z_axis() - Vector3::UnitX()).norm() <= 1e-15);
  REQUIRE((f.y_axis() + Vector3::UnitY()).norm() <= 1e-15);
  // Oblique: z' in the vertical plane through the axis, pointing up.
  f = frame_of(Vector3(1.0, 1.0, 1.0), Vector3::Zero());
  REQUIRE(f.z_axis().z() > 0.0);
  REQUIRE(std::abs(f.y_axis().z()) <= 1e-15);
  REQUIRE((f.rotation * f.rotation.transpose() - Matrix3::Identity()).norm() <= 1e-15);
  REQUIRE(f.rotation.determinant() == Approx(1.0).epsilon(1e-15));
  // An orientation vector: its part normal to the axis is y'.
  f = frame_of(Vector3(2.0, 0.0, 0.0), Vector3(0.5, 0.0, 1.0));
  REQUIRE((f.y_axis() - Vector3::UnitZ()).norm() <= 1e-15);
  REQUIRE((f.z_axis() + Vector3::UnitY()).norm() <= 1e-15);
  // Along the axis it cannot fix y'; nodes that coincide give no axis.
  REQUIRE_THROWS_WITH(frame_of(Vector3(2.0, 0.0, 0.0), Vector3(1.0, 0.0, 0.0)),
                      ContainsSubstring("lies along the element's axis"));
  REQUIRE_THROWS_WITH(frame_of(Vector3::Zero(), Vector3::Zero()),
                      ContainsSubstring("zero length"));
}

TEST_CASE("beam sections: the shapes' properties and Cowper's shear coefficients",
          "[beam][section]") {
  // A square: J = 0.1406 a^4 (Timoshenko and Goodier), I = a^4 / 12.
  BeamSection sq;
  sq.name = "square";
  sq.shape = BeamSectionShape::Rectangle;
  sq.width = sq.height = 0.02;
  const BeamSection r = resolve(sq, 0.3);
  REQUIRE(r.area == Approx(4.0e-4).epsilon(1e-15));
  REQUIRE(r.iy == Approx(std::pow(0.02, 4) / 12.0).epsilon(1e-14));
  REQUIRE(r.torsion / std::pow(0.02, 4) == Approx(0.140577).margin(1e-6));
  REQUIRE(r.shear_y == Approx(10.0 * 1.3 / (12.0 + 3.3)).epsilon(1e-15));
  REQUIRE(r.fibre_z == Approx(0.01).epsilon(1e-15));
  // A thin strip b >> h: J -> b h^3 / 3 (1 - 0.630 h / b).
  REQUIRE(rectangle_torsion_constant(1.0, 0.01) / (1.0 * 1e-6 / 3.0) ==
          Approx(1.0 - 0.630 * 0.01).margin(2e-5));
  // Width along y', height along z': I_y = b h^3 / 12.
  BeamSection rect = sq;
  rect.width = 0.01;
  rect.height = 0.03;
  const BeamSection rr = resolve(rect, 0.3);
  REQUIRE(rr.iy == Approx(0.01 * std::pow(0.03, 3) / 12.0).epsilon(1e-14));
  REQUIRE(rr.iz == Approx(0.03 * std::pow(0.01, 3) / 12.0).epsilon(1e-14));
  REQUIRE(rectangle_torsion_constant(0.01, 0.03) == rectangle_torsion_constant(0.03, 0.01));
  // Circle and tube; a thin tube's coefficient tends to 2 (1 + nu) / (4 + 3 nu).
  BeamSection c;
  c.name = "rod";
  c.shape = BeamSectionShape::Circle;
  c.radius = 0.01;
  const BeamSection rc = resolve(c, 0.3);
  REQUIRE(rc.torsion == Approx(0.5 * kPi * 1e-8).epsilon(1e-14));
  REQUIRE(rc.polar() == Approx(rc.torsion).epsilon(1e-14));
  REQUIRE(rc.shear_z == Approx(6.0 * 1.3 / 8.8).epsilon(1e-15));
  REQUIRE(cowper_tube(0.3, 0.0) == Approx(cowper_circle(0.3)).epsilon(1e-15));
  REQUIRE(cowper_tube(0.3, 1.0) == Approx(2.0 * 1.3 / 4.9).epsilon(1e-15));
  BeamSection t;
  t.name = "tube";
  t.shape = BeamSectionShape::Tube;
  t.radius = 0.02;
  t.inner_radius = 0.018;
  const BeamSection rt = resolve(t, 0.3);
  REQUIRE(rt.area == Approx(kPi * (4e-4 - 3.24e-4)).epsilon(1e-13));
  t.inner_radius = 0.02;
  REQUIRE_THROWS_WITH(resolve(t, 0.3), ContainsSubstring("smaller than the outer radius"));
  // A general section states its shear coefficients.
  BeamSection gen;
  gen.name = "gen";
  gen.area = 1e-3;
  gen.iy = gen.iz = 1e-7;
  gen.torsion = 2e-7;
  REQUIRE_THROWS_WITH(resolve(gen, 0.3), ContainsSubstring("states its shear coefficients"));
  gen.shear_y = gen.shear_z = 1.2;
  REQUIRE_THROWS_WITH(resolve(gen, 0.3), ContainsSubstring("must lie in (0, 1]"));
  REQUIRE(parse_beam_section_shape("tube") == BeamSectionShape::Tube);
  REQUIRE_THROWS_AS(parse_beam_section_shape("I"), ConfigError);
}

TEST_CASE("one beam element is exact for a cantilever under end loads, and its end forces "
          "balance",
          "[beam][element]") {
  // A cantilever along an oblique axis, clamped at node 0: K_11 u_1 = f_1.
  const Beam2Element element;
  const BeamSection s = general_section();
  const Vector3 x0(0.3, -0.2, 0.1);
  const Vector3 x1(1.2, 0.4, 0.9);
  const Scalar l = (x1 - x0).norm();
  const Matrix g = Beam2Element::geometry(x0, x1, s, kE, kG);
  const BeamFrame f = Beam2Element::frame(g);
  const Matrix k = element.stiffness(g, material(), 1.0, IntegrationOptions());
  const Matrix k11 = k.bottomRightCorner(6, 6);
  const Scalar p = 1000.0;
  const Scalar t = 50.0;
  // A shear along y' and a torque: w along y' = P L^3/(3 E I_z) + P L/(k_y G A).
  Vector load = Vector::Zero(6);
  load.head<3>() = p * f.y_axis();
  load.tail<3>() = t * f.x_axis();
  const Vector u1 = k11.ldlt().solve(load);
  const Scalar deflection = p * l * l * l / (3.0 * kE * s.iz) + p * l / (s.shear_y * kG * s.area);
  REQUIRE(u1.head<3>().dot(f.y_axis()) == Approx(deflection).epsilon(1e-11));
  REQUIRE(u1.tail<3>().dot(f.x_axis()) == Approx(t * l / (kG * s.torsion)).epsilon(1e-11));
  // The rotation at the tip about z': P L^2 / (2 E I_z).
  REQUIRE(u1.tail<3>().dot(f.z_axis()) == Approx(p * l * l / (2.0 * kE * s.iz)).epsilon(1e-11));
  // End forces: the root carries the shear P, the torque T and the moment
  // P L; the tip is free of moment.
  Vector ue = Vector::Zero(12);
  ue.tail<6>() = u1;
  const BeamEndForces ef = element.end_forces(g, material(), ue, Vector3::Zero());
  REQUIRE(ef.end(1) == Approx(p).epsilon(1e-10));
  REQUIRE(ef.start(1) == Approx(p).epsilon(1e-10));
  REQUIRE(ef.end(3) == Approx(t).epsilon(1e-10));
  REQUIRE(std::abs(ef.end(5)) <= 1e-9 * p * l);
  // M_z = E I_z theta_z' = P L at the root: the fibres at y' > 0 in
  // compression, the side the load bends the beam towards.
  REQUIRE(ef.start(5) == Approx(p * l).epsilon(1e-10));
}

namespace {

/// A cantilever 2 m along x in four elements, clamped at A: a tip force and
/// moment at B, and a uniform load along it.
const char* kFrameDeck = R"({
  "name": "frame_unit",
  "mesh": { "type": "frame",
            "points": [ { "name": "A", "position": [0, 0, 0] }, { "name": "B", "position": [2, 0, 0] } ],
            "members": [ { "name": "span", "from": "A", "to": "B", "elements": 4 } ] },
  "material": { "youngs_modulus": 70e9, "poisson_ratio": 0.3, "density": 2700 },
  "model": { "beam": { "sections": [
    { "name": "bar", "shape": "rectangle", "width": 0.05, "height": 0.1 } ] } },
  "boundary_conditions": [
    { "name": "clamp", "fix": ["x", "y", "z", "rx", "ry", "rz"], "region": { "group": "A" } } ],
  "load_cases": [
    { "name": "tip", "point_loads": [
        { "force": [0, 0, -1000], "moment": [300, 0, 0], "region": { "group": "B" } } ] },
    { "name": "uniform", "line_loads": [
        { "force_per_length": [0, 0, -2000], "region": { "group": "span" } } ] }
  ]
})";

std::string read_text(const std::string& path) {
  std::ifstream in(path);
  std::stringstream buffer;
  buffer << in.rdbuf();
  return buffer.str();
}

}  // namespace

TEST_CASE("beam decks: a frame with sections, moments and line loads solves exactly; what a "
          "beam refuses is refused",
          "[beam][io][config]") {
  const Configuration config = parse_configuration(json::parse(kFrameDeck, "frame.json"), "frame", true);
  REQUIRE(config.mesh_kind == MeshKind::Frame);
  REQUIRE(config.is_beam());
  REQUIRE(config.describe_mesh() == "frame of 1 member(s), 4 elements");
  REQUIRE(config.load_cases[1].line_loads.size() == 1);
  FemModel model = build_model(config);
  REQUIRE(model.mesh().element_type() == ElementType::Beam2);
  REQUIRE(model.dofs_per_node() == 6);
  Assembler assembler(model);
  StaticAnalysis analysis(model, assembler);
  const std::vector<StaticSolution> solutions = analysis.solve_all();
  // The rectangle is 0.05 wide along y' = y and 0.1 high along z' = z.
  const BeamSection s = model.section_of(0);
  const Scalar l = 2.0;
  const Scalar e = 70e9;
  const Scalar g = e / 2.6;
  const Scalar kga = s.shear_z * g * s.area;
  const Scalar eiy = e * 0.05 * 0.001 / 12.0;
  const Index tip = model.mesh().node_sets().at("B").front();
  const Vector& u_tip = solutions[0].displacement;
  REQUIRE(u_tip(6 * tip + 2) ==
          Approx(-1000.0 * (l * l * l / (3.0 * eiy) + l / kga)).epsilon(1e-11));
  REQUIRE(u_tip(6 * tip + 3) == Approx(300.0 * l / (g * s.torsion)).epsilon(1e-11));
  const Vector& u_line = solutions[1].displacement;
  REQUIRE(u_line(6 * tip + 2) ==
          Approx(-2000.0 * (l * l * l * l / (8.0 * eiy) + l * l / (2.0 * kga))).epsilon(1e-11));
  // The root's resultants: M_y = P L stretches the top fibres, and the
  // extreme-fibre stress is M c / I.
  const BeamField field = recover_beam_forces(model, assembler, u_tip, model.load_case_specs()[0]);
  const Index root_element = 0;
  REQUIRE(model.mesh().element_nodes(root_element)[0] == model.mesh().node_sets().at("A").front());
  REQUIRE(field.element[0].start(2) == Approx(-1000.0).epsilon(1e-10));
  REQUIRE(field.element[0].start(3) == Approx(300.0).epsilon(1e-10));
  REQUIRE(field.element[0].start(4) == Approx(1000.0 * l).epsilon(1e-10));
  REQUIRE(field.element_normal_stress(0) == Approx(1000.0 * l * 0.05 / (0.05 * 0.001 / 12.0)).epsilon(1e-10));

  // The result files: the end forces per element, and a mesh.json with the
  // sections and the distributed loads per element.
  ensure_directory("results/_test_tmp/beam_deck");
  ResultWriter writer("results/_test_tmp/beam_deck", config);
  writer.write_mesh(model);
  writer.write_beam_forces(model, "tip", field);
  writer.write_beam_vtk(model, "tip", u_tip, field);
  const json::Value mesh_json = json::parse_file("results/_test_tmp/beam_deck/mesh.json");
  const json::Value* beam = mesh_json.find("beam");
  REQUIRE(beam != nullptr);
  REQUIRE(beam->find("elements")->array_items().size() == 4);
  const json::Value& uniform = beam->find("load_cases")->array_items()[1];
  REQUIRE(uniform.find("distributed_N_per_m")->array_items()[2].array_items()[2].number_value() == -2000.0);
  REQUIRE(uniform.find("nodal_loads_node_fx_fy_fz_mx_my_mz")->array_items().empty());
  REQUIRE_THAT(read_text("results/_test_tmp/beam_deck/beam_tip.csv"),
               ContainsSubstring("N0[N],Qy0[N],Qz0[N],T0[Nm],My0[Nm],Mz0[Nm]"));
  REQUIRE_THAT(read_text("results/_test_tmp/beam_deck/fields_tip.vtk"),
               ContainsSubstring("CELL_TYPES 4\n3\n3\n3\n3\n"));

  // What a beam does not take is refused with the reason.
  const auto error = [](const std::string& text) {
    try {
      parse_configuration(json::parse(text, "deck.json"), "deck", true);
    } catch (const std::exception& e) {
      return std::string(e.what());
    }
    return std::string();
  };
  const auto edited = [](std::string text, const std::string& from, const std::string& to) {
    const std::size_t at = text.find(from);
    REQUIRE(at != std::string::npos);
    return text.replace(at, from.size(), to);
  };
  const std::string deck = kFrameDeck;
  const auto with = [&](const std::string& extra) {
    std::string text = deck;
    text.insert(text.rfind('}'), "," + extra);
    return error(text);
  };
  REQUIRE_THAT(with(R"("nonlinear": { "enabled": true })"), ContainsSubstring("small rotations"));
  REQUIRE_THAT(with(R"("topology": { "enabled": true })"), ContainsSubstring("section"));
  REQUIRE_THAT(error(edited(deck, R"("model": {)", R"("model": { "thickness": 0.1,)")),
               ContainsSubstring("does not apply to the beam mesh"));
  REQUIRE_THAT(error(edited(deck, R"("width": 0.05,)", R"("width": 0.05, "radius": 0.1,)")),
               ContainsSubstring("'model.beam.sections[0].radius' does not apply to a rectangle"));
  REQUIRE_THAT(error(edited(deck, R"("width": 0.05,)", R"("width": -0.05,)")),
               ContainsSubstring("the width must be positive"));
  REQUIRE_THAT(error(edited(deck, R"("to": "B", "elements": 4)", R"("to": "D", "elements": 4)")),
               ContainsSubstring("D"));
  REQUIRE_THAT(error(edited(deck, R"("line_loads": [)", R"("tractions": [ { "traction": [0, 0, 1], "region": { "group": "B" } } ], "line_loads": [)")),
               ContainsSubstring("a beam has no faces"));
  REQUIRE_THAT(error(edited(deck, R"({ "beam": { "sections": [
    { "name": "bar", "shape": "rectangle", "width": 0.05, "height": 0.1 } ] } })", R"({ })")),
               ContainsSubstring("'model.beam.sections' is required"));
  REQUIRE_THAT(error(edited(deck, R"("shape": "rectangle", "width": 0.05, "height": 0.1)",
                            R"("shape": "general", "area": 1e-3, "iy": 1e-6, "iz": 1e-6, "torsion": 2e-6)")),
               ContainsSubstring("states its shear coefficients"));
  // A section whose region selects nothing, and elements without a section.
  std::string partial = edited(deck, R"("name": "bar",)", R"("name": "bar", "region": { "box": { "xmax": 0.5 } },)");
  REQUIRE_THROWS_WITH(build_model(parse_configuration(json::parse(partial, "p.json"), "p", true)),
                      ContainsSubstring("3 of 4 beam elements have no cross-section"));
  std::string empty = edited(deck, R"("name": "bar",)", R"("name": "bar", "region": { "box": { "xmin": 5 } },)");
  REQUIRE_THROWS_WITH(build_model(parse_configuration(json::parse(empty, "e.json"), "e", true)),
                      ContainsSubstring("selected no element"));
  // Line loads belong to beams.
  const std::string solid = R"({
    "mesh": { "type": "structured_hex", "nx": 2, "ny": 1, "nz": 1, "lx": 1, "ly": 1, "lz": 1 },
    "material": { "youngs_modulus": 1e9, "poisson_ratio": 0.3 },
    "boundary_conditions": [ { "fix": ["x", "y", "z"], "region": { "box": { "xmax": 0 } } } ],
    "load_cases": [ { "name": "q", "line_loads": [ { "force_per_length": [0, 0, 1], "region": { "box": { "xmin": 1 } } } ] } ]
  })";
  REQUIRE_THAT(error(solid), ContainsSubstring("is a load per unit length along beams"));
  REQUIRE_THAT(error(edited(solid, R"("load_cases")", R"("model": { "stress_state": "beam" }, "load_cases")")),
               ContainsSubstring("needs a beam mesh"));
}

TEST_CASE("beam meshes from files: B31 elements and Gmsh lines become beams; trusses do not",
          "[beam][io][unstructured]") {
  const std::string inp =
      "*NODE\n1, 0, 0, 0\n2, 1, 0, 0\n3, 1, 1, 0\n4, 1, 1, 1\n"
      "*ELEMENT, TYPE=B31, ELSET=FRAME\n11, 1, 2\n12, 2, 3\n"
      "*ELEMENT, TYPE=B31H, ELSET=POST\n13, 3, 4\n"
      "*NSET, NSET=ROOT\n1\n"
      "*BEAM SECTION, ELSET=FRAME, MATERIAL=STEEL, SECTION=RECT\n0.1, 0.1\n0, 0, 1\n";
  std::istringstream in(inp);
  MeshReadReport report;
  const Mesh mesh = read_abaqus_inp(in, "frame.inp", ".", MeshReadOptions(), &report);
  REQUIRE(mesh.element_type() == ElementType::Beam2);
  REQUIRE(mesh.dim() == 3);
  REQUIRE(mesh.num_elements() == 3);
  REQUIRE(mesh.element_sets().at("FRAME") == std::vector<Index>({0, 1}));
  REQUIRE(mesh.element_sets().at("POST") == std::vector<Index>({2}));
  REQUIRE(mesh.node_sets().at("ROOT") == std::vector<Index>({0}));
  REQUIRE(mesh.element_measure(2) == Approx(1.0));
  const auto read_inp = [](const std::string& text, bool beam) {
    std::istringstream stream(text);
    MeshReadOptions options;
    options.beam = beam;
    return read_abaqus_inp(stream, "bad.inp", ".", options);
  };
  const std::string truss = "*NODE\n1, 0, 0, 0\n2, 1, 0, 0\n*ELEMENT, TYPE=T3D2\n1, 1, 2\n";
  REQUIRE_THROWS_WITH(read_inp(truss, false), ContainsSubstring("only points and lines"));
  REQUIRE_THROWS_WITH(read_inp(truss, true), ContainsSubstring("carries no bending"));
  REQUIRE_THROWS_WITH(read_inp("*NODE\n1, 0, 0, 0\n2, 1, 0, 0\n3, 2, 0, 0\n"
                               "*ELEMENT, TYPE=B31\n1, 1, 2\n*ELEMENT, TYPE=T3D2\n2, 2, 3\n",
                               false),
                      ContainsSubstring("mixes beams with other lines"));
  REQUIRE_THROWS_WITH(read_inp("*NODE\n1, 0, 0, 0\n2, 1, 0, 0\n3, 0.5, 0, 0\n"
                               "*ELEMENT, TYPE=B32\n1, 1, 2, 3\n",
                               true),
                      ContainsSubstring("B31"));

  // A Gmsh frame: 2-node lines in a physical curve, a physical point.
  const std::string msh =
      "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n"
      "$PhysicalNames\n2\n0 1 \"base\"\n1 2 \"legs\"\n$EndPhysicalNames\n"
      "$Nodes\n3\n1 0 0 0\n2 0 0 1\n3 1 0 1\n$EndNodes\n"
      "$Elements\n3\n1 15 2 1 1 1\n2 1 2 2 1 1 2\n3 1 2 2 2 2 3\n$EndElements\n";
  MeshReadOptions as_beam;
  as_beam.beam = true;
  std::istringstream gmsh(msh);
  const Mesh lines = read_gmsh(gmsh, "frame.msh", as_beam);
  REQUIRE(lines.element_type() == ElementType::Beam2);
  REQUIRE(lines.num_elements() == 2);
  REQUIRE(lines.element_sets().at("legs") == std::vector<Index>({0, 1}));
  REQUIRE(lines.node_sets().at("base") == std::vector<Index>({0}));
  std::istringstream plain(msh);
  REQUIRE_THROWS_WITH(read_gmsh(plain, "frame.msh", MeshReadOptions()),
                      ContainsSubstring("mesh.beam = true"));
  std::istringstream triangle("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n3\n1 0 0 0\n2 1 0 0\n"
                              "3 0 1 0\n$EndNodes\n$Elements\n1\n1 2 0 1 2 3\n$EndElements\n");
  REQUIRE_THROWS_WITH(read_gmsh(triangle, "tri.msh", as_beam),
                      ContainsSubstring("a frame is made of lines only"));
}

TEST_CASE("CalculiX beam decks: B31 with a RECT section per y' axis; other sections refused",
          "[beam][io][cross-validation]") {
  ensure_directory("results/_test_tmp");
  FrameMeshSpec spec;
  spec.points = {{"A", Vector3::Zero()}, {"B", Vector3(1.0, 0.0, 0.0)}, {"C", Vector3(1.0, 0.0, 1.0)}};
  FrameMember ab;
  ab.name = "AB";
  ab.from = "A";
  ab.to = "B";
  ab.elements = 2;
  FrameMember bc = ab;
  bc.name = "BC";
  bc.from = "B";
  bc.to = "C";
  spec.members = {ab, bc};
  const auto model_with = [&](const BeamSection& section) {
    FemModel model(make_frame_mesh(spec), IsotropicMaterial(210.0e9, 0.3, 7850.0), 1.0, StressState::Beam,
                   IntegrationOptions());
    std::vector<Index> all(4);
    std::iota(all.begin(), all.end(), Index{0});
    model.assign_section(section, all);
    DisplacementConstraint clamp;
    clamp.region.name = "A";
    clamp.region.members.emplace_back();
    clamp.region.members.back().kind = SelectorKind::Group;
    clamp.region.members.back().group = "A";
    for (int k = 0; k < 6; ++k) clamp.set(k, true);
    model.constraints() = {clamp};
    LoadCaseSpec lc;
    lc.name = "tip";
    PointLoadSpec p;
    p.region.name = "C";
    p.region.members.emplace_back();
    p.region.members.back().kind = SelectorKind::Group;
    p.region.members.back().group = "C";
    p.force = Vector3(0.0, 100.0, 0.0);
    p.moment = Vector3(0.0, 0.0, 7.0);
    lc.point_loads.push_back(p);
    model.load_case_specs() = {lc};
    model.finalize();
    return model;
  };
  BeamSection rect;
  rect.name = "rect";
  rect.shape = BeamSectionShape::Rectangle;
  rect.width = 0.02;
  rect.height = 0.04;
  const FemModel model = model_with(rect);
  REQUIRE(calculix_element_type(model) == "B31");
  const std::vector<std::string> decks = write_calculix_decks(model, "results/_test_tmp/beam", "unit");
  REQUIRE(decks.size() == 1);
  const std::string text = read_text(decks.front());
  // The points are nodes 1 to 3, the members' inner nodes follow.
  REQUIRE_THAT(text, ContainsSubstring("*ELEMENT, TYPE=B31, ELSET=EALL\n1, 1, 4\n2, 4, 2\n"));
  // The horizontal member's y' is y, the vertical one's -y (z' = x there):
  // one section each, the side along y' first.
  REQUIRE_THAT(text, ContainsSubstring("*BEAM SECTION, ELSET=B1, MATERIAL=MAT1, SECTION=RECT\n"
                                       "0.02, 0.040000000000000001\n0, 1, 0\n"));
  REQUIRE_THAT(text, ContainsSubstring("*BEAM SECTION, ELSET=B2, MATERIAL=MAT1, SECTION=RECT\n"
                                       "0.02, 0.040000000000000001\n0, -1, 0\n"));
  // The force on DOF 2 and the moment on DOF 6 of the tip C (node 3).
  REQUIRE_THAT(text, ContainsSubstring("\n3, 2, 100\n"));
  REQUIRE_THAT(text, ContainsSubstring("\n3, 6, 7\n"));
  REQUIRE_THAT(text, ContainsSubstring("*NODE FILE, OUTPUT=2D\nU\n"));
  std::remove(decks.front().c_str());
  BeamSection tube;
  tube.name = "tube";
  tube.shape = BeamSectionShape::Tube;
  tube.radius = 0.02;
  tube.inner_radius = 0.015;
  REQUIRE_THROWS_WITH(write_calculix_decks(model_with(tube), "results/_test_tmp/beam", "unit"),
                      ContainsSubstring("only rectangular sections are exported"));
  BeamSection bernoulli = rect;
  bernoulli.shear_deformation = false;
  REQUIRE_THROWS_WITH(write_calculix_decks(model_with(bernoulli), "results/_test_tmp/beam", "unit"),
                      ContainsSubstring("no shear deformation"));
}
