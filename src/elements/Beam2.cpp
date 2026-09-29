#include "sparlab/elements/Beam2.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/elements/Quadrature.hpp"

#include <Eigen/Geometry>

#include <cmath>
#include <sstream>

namespace sparlab {
namespace {

constexpr int kDofs = 12;
constexpr int kGeometryRows = 14;
using Operator = Eigen::Matrix<Scalar, 6, kDofs>;
using Matrix12 = Eigen::Matrix<Scalar, kDofs, kDofs>;
using Vector12 = Eigen::Matrix<Scalar, kDofs, 1>;
using Vector6 = Eigen::Matrix<Scalar, 6, 1>;

/// Everything the kernels need, read from the geometry matrix.
struct BeamData {
  BeamFrame frame;
  Scalar area = 0.0;
  Scalar iy = 0.0;
  Scalar iz = 0.0;
  Scalar torsion = 0.0;
  Scalar ky = 0.0;
  Scalar kz = 0.0;
  Scalar e = 0.0;
  Scalar g = 0.0;
  /// Phi_y = 12 E I_z / (k_y G A L^2) and Phi_z = 12 E I_y / (k_z G A L^2),
  /// 0 without shear deformation.
  Scalar phi_y = 0.0;
  Scalar phi_z = 0.0;
};

BeamData beam_data(const Matrix& geometry) {
  BeamData b;
  b.frame = Beam2Element::frame(geometry);
  b.area = geometry(6, 0);
  b.iy = geometry(7, 0);
  b.iz = geometry(8, 0);
  b.torsion = geometry(9, 0);
  b.ky = geometry(10, 0);
  b.kz = geometry(11, 0);
  b.e = geometry(12, 0);
  b.g = geometry(13, 0);
  if (!(b.area > 0.0) || !(b.iy > 0.0) || !(b.iz > 0.0) || !(b.torsion > 0.0) ||
      !(b.ky >= 0.0) || !(b.kz >= 0.0) || !(b.e > 0.0) || !(b.g > 0.0)) {
    throw ModelError("Beam2 needs a positive area, second moments, torsion constant and "
                     "moduli, and non-negative shear coefficients (rows 6-13 of its geometry)");
  }
  const Scalar l2 = b.frame.length * b.frame.length;
  b.phi_y = b.ky > 0.0 ? 12.0 * b.e * b.iz / (b.ky * b.g * b.area * l2) : 0.0;
  b.phi_z = b.kz > 0.0 ? 12.0 * b.e * b.iy / (b.kz * b.g * b.area * l2) : 0.0;
  return b;
}

/// The constitutive matrix must be the geometry's material.
void check_material(const Matrix& d, const BeamData& b) {
  if (d.rows() != 2 || d.cols() != 2) {
    throw ModelError("Beam2 takes the 2 x 2 constitutive matrix diag(E, G) of the beam "
                     "stress state");
  }
  if (std::abs(d(0, 0) - b.e) > 1e-12 * b.e || std::abs(d(1, 1) - b.g) > 1e-12 * b.g) {
    throw ModelError("Beam2: the constitutive matrix diag(E, G) differs from the moduli in "
                     "the element's geometry");
  }
}

/// The four deflection and four rotation functions of one bending plane and
/// their xi-derivatives: v = sum n_a q_a, theta = sum m_a q_a over
/// (v0, theta0, v1, theta1).
struct Bending {
  Eigen::Vector4d n, dn, m, dm;
};

Bending bending(Scalar xi, Scalar length, Scalar phi) {
  const Scalar mu = 1.0 / (1.0 + phi);
  const Scalar x2 = xi * xi;
  const Scalar x3 = x2 * xi;
  const Scalar l = length;
  Bending b;
  b.n << mu * (1.0 - 3.0 * x2 + 2.0 * x3 + phi * (1.0 - xi)),
      mu * l * (xi - 2.0 * x2 + x3 + 0.5 * phi * (xi - x2)),
      mu * (3.0 * x2 - 2.0 * x3 + phi * xi),
      mu * l * (-x2 + x3 - 0.5 * phi * (xi - x2));
  b.dn << mu * (-6.0 * xi + 6.0 * x2 - phi),
      mu * l * (1.0 - 4.0 * xi + 3.0 * x2 + 0.5 * phi * (1.0 - 2.0 * xi)),
      mu * (6.0 * xi - 6.0 * x2 + phi),
      mu * l * (-2.0 * xi + 3.0 * x2 - 0.5 * phi * (1.0 - 2.0 * xi));
  b.m << 6.0 * mu / l * (x2 - xi),
      mu * (1.0 - 4.0 * xi + 3.0 * x2 + phi * (1.0 - xi)),
      6.0 * mu / l * (xi - x2),
      mu * (-2.0 * xi + 3.0 * x2 + phi * xi);
  b.dm << 6.0 * mu / l * (2.0 * xi - 1.0),
      mu * (-4.0 + 6.0 * xi - phi),
      6.0 * mu / l * (1.0 - 2.0 * xi),
      mu * (-2.0 + 6.0 * xi + phi);
  return b;
}

/// The interpolation H (u, v, w, theta_x, theta_y, theta_z in local axes in
/// terms of the 12 local DOFs) and its derivative along x' at xi in [0, 1].
struct Interpolation {
  Operator h = Operator::Zero();
  Operator dh = Operator::Zero();
};

Interpolation interpolation(Scalar xi, const BeamData& b) {
  const Scalar length = b.frame.length;
  Interpolation in;
  const Scalar dx = 1.0 / length;  // d/dx' = (1 / L) d/dxi
  // Axial displacement and twist: linear.
  in.h(0, 0) = 1.0 - xi;
  in.h(0, 6) = xi;
  in.dh(0, 0) = -dx;
  in.dh(0, 6) = dx;
  in.h(3, 3) = 1.0 - xi;
  in.h(3, 9) = xi;
  in.dh(3, 3) = -dx;
  in.dh(3, 9) = dx;
  // The x'-y' plane: (v, theta_z) over (v0, theta_z0, v1, theta_z1).
  const Bending y = bending(xi, length, b.phi_y);
  const int vy[4] = {1, 5, 7, 11};
  for (int a = 0; a < 4; ++a) {
    in.h(1, vy[a]) = y.n(a);
    in.dh(1, vy[a]) = y.dn(a) * dx;
    in.h(5, vy[a]) = y.m(a);
    in.dh(5, vy[a]) = y.dm(a) * dx;
  }
  // The x'-z' plane: (w, psi) over (w0, psi0, w1, psi1) with psi = -theta_y,
  // so the rotation DOFs enter with a minus sign, and theta_y = -psi.
  const Bending z = bending(xi, length, b.phi_z);
  const int wz[4] = {2, 4, 8, 10};
  const Scalar sign[4] = {1.0, -1.0, 1.0, -1.0};
  for (int a = 0; a < 4; ++a) {
    in.h(2, wz[a]) = sign[a] * z.n(a);
    in.dh(2, wz[a]) = sign[a] * z.dn(a) * dx;
    in.h(4, wz[a]) = -sign[a] * z.m(a);
    in.dh(4, wz[a]) = -sign[a] * z.dm(a) * dx;
  }
  return in;
}

/// The generalised strains (u', gamma_y, gamma_z, theta_x', theta_y',
/// theta_z') on the local DOFs.
Operator strain_rows(const Interpolation& in) {
  Operator s;
  s.row(0) = in.dh.row(0);
  s.row(1) = in.dh.row(1) - in.h.row(5);
  s.row(2) = in.dh.row(2) + in.h.row(4);
  s.row(3) = in.dh.row(3);
  s.row(4) = in.dh.row(4);
  s.row(5) = in.dh.row(5);
  return s;
}

/// The 12 x 12 rotation of the DOFs to local axes (local = T global).
Matrix12 transformation(const Matrix3& r) {
  Matrix12 t = Matrix12::Zero();
  for (int k = 0; k < 4; ++k) t.block<3, 3>(3 * k, 3 * k) = r;
  return t;
}

/// Four Gauss points on [-1, 1], mapped to [0, 1]: exact for the degree-6
/// mass integrand, and for every other one here.
const std::vector<QuadraturePoint1D>& rule() { return gauss_legendre_line(4); }

/// The section stiffnesses diag(EA, k_y GA, k_z GA, GJ, EI_y, EI_z); a
/// direction without shear deformation has none (its shear strain is zero).
Vector6 section_stiffness(const BeamData& b) {
  Vector6 dg;
  dg << b.e * b.area, b.ky * b.g * b.area, b.kz * b.g * b.area, b.g * b.torsion, b.e * b.iy,
      b.e * b.iz;
  return dg;
}

Matrix12 local_stiffness(const BeamData& b) {
  const Scalar l = b.frame.length;
  const Vector6 dg = section_stiffness(b);
  Matrix12 k = Matrix12::Zero();
  for (const auto& gp : rule()) {
    const Operator bs = strain_rows(interpolation(0.5 * (gp.xi + 1.0), b));
    k.noalias() += (0.5 * gp.weight * l) * (bs.transpose() * dg.asDiagonal() * bs);
  }
  return 0.5 * (k + k.transpose());
}

/// int H'^T diag(1, 1, 1, I_p/A, I_y/A, I_z/A) H' dx': the geometric
/// stiffness per unit axial force.
Matrix12 local_geometric(const BeamData& b) {
  const Scalar l = b.frame.length;
  Vector6 w;
  w << 1.0, 1.0, 1.0, (b.iy + b.iz) / b.area, b.iy / b.area, b.iz / b.area;
  Matrix12 kg = Matrix12::Zero();
  for (const auto& gp : rule()) {
    const Interpolation in = interpolation(0.5 * (gp.xi + 1.0), b);
    kg.noalias() += (0.5 * gp.weight * l) * (in.dh.transpose() * w.asDiagonal() * in.dh);
  }
  return 0.5 * (kg + kg.transpose());
}

/// Consistent local nodal forces of a uniform line load q (local components).
Vector12 local_line_load(const BeamData& b, const Vector3& q_local) {
  const Scalar l = b.frame.length;
  Vector12 f = Vector12::Zero();
  for (const auto& gp : rule()) {
    const Interpolation in = interpolation(0.5 * (gp.xi + 1.0), b);
    f.noalias() += (0.5 * gp.weight * l) * (in.h.topRows(3).transpose() * q_local);
  }
  return f;
}

Vector12 element_dofs(const Vector& ue) {
  if (ue.size() != kDofs) {
    std::ostringstream os;
    os << "Beam2 takes 12 element DOFs, received " << ue.size();
    throw ModelError(os.str());
  }
  return ue;
}

}  // namespace

BeamFrame Beam2Element::frame(const Matrix& geometry) {
  if (geometry.cols() != 2 || geometry.rows() != kGeometryRows) {
    std::ostringstream os;
    os << "Beam2 expects a 14 x 2 geometry matrix (coordinates, orientation, section, "
          "moduli), received "
       << geometry.rows() << " x " << geometry.cols();
    throw MeshError(os.str());
  }
  const Vector3 x0 = geometry.block<3, 1>(0, 0);
  const Vector3 x1 = geometry.block<3, 1>(0, 1);
  BeamFrame f;
  f.length = (x1 - x0).norm();
  if (!(f.length > 0.0) || !std::isfinite(f.length)) {
    throw MeshError("Beam2 element of zero length: its two nodes coincide");
  }
  const Vector3 ex = (x1 - x0) / f.length;
  const Vector3 v = geometry.block<3, 1>(3, 0);
  Vector3 ey;
  Vector3 ez;
  if (v.norm() > 0.0) {
    const Vector3 p = v - v.dot(ex) * ex;
    if (!(p.norm() > 1.0e-6 * v.norm())) {
      throw MeshError("Beam2 orientation vector lies along the element's axis: it cannot "
                      "fix the section's y' axis");
    }
    ey = p.normalized();
    ez = ex.cross(ey);
  } else {
    // z' the projection of global Z ("up"), or of global X for an element
    // within 0.1 degree of vertical.
    static const Scalar kParallel = std::cos(0.1 * 3.14159265358979323846 / 180.0);
    const Vector3 up = std::abs(ex.z()) > kParallel ? Vector3::UnitX() : Vector3::UnitZ();
    ez = (up - up.dot(ex) * ex).normalized();
    ey = ez.cross(ex);
  }
  f.rotation.row(0) = ex.transpose();
  f.rotation.row(1) = ey.transpose();
  f.rotation.row(2) = ez.transpose();
  return f;
}

Matrix Beam2Element::geometry(const Vector3& x0, const Vector3& x1, const BeamSection& section,
                              Scalar youngs_modulus, Scalar shear_modulus) {
  Matrix g(kGeometryRows, 2);
  for (int c = 0; c < 2; ++c) {
    g.block<3, 1>(0, c) = c == 0 ? x0 : x1;
    g.block<3, 1>(3, c) = section.orientation;
    g(6, c) = section.area;
    g(7, c) = section.iy;
    g(8, c) = section.iz;
    g(9, c) = section.torsion;
    g(10, c) = section.shear_deformation ? section.shear_y : 0.0;
    g(11, c) = section.shear_deformation ? section.shear_z : 0.0;
    g(12, c) = youngs_modulus;
    g(13, c) = shear_modulus;
  }
  return g;
}

Matrix Beam2Element::stiffness(const Matrix& geometry, const Matrix& d, Scalar,
                               const IntegrationOptions&) const {
  const BeamData b = beam_data(geometry);
  check_material(d, b);
  const Matrix12 t = transformation(b.frame.rotation);
  return t.transpose() * local_stiffness(b) * t;
}

Matrix Beam2Element::consistent_mass(const Matrix& geometry, Scalar density, Scalar,
                                     const IntegrationOptions&) const {
  const BeamData b = beam_data(geometry);
  const Scalar l = b.frame.length;
  Vector6 rho;
  rho << b.area, b.area, b.area, b.iy + b.iz, b.iy, b.iz;
  rho *= density;
  Matrix12 m = Matrix12::Zero();
  for (const auto& gp : rule()) {
    const Interpolation in = interpolation(0.5 * (gp.xi + 1.0), b);
    m.noalias() += (0.5 * gp.weight * l) * (in.h.transpose() * rho.asDiagonal() * in.h);
  }
  const Matrix12 t = transformation(b.frame.rotation);
  return t.transpose() * (0.5 * (m + m.transpose())) * t;
}

Matrix Beam2Element::lumped_mass(const Matrix& geometry, Scalar density, Scalar,
                                 const IntegrationOptions&) const {
  const BeamData b = beam_data(geometry);
  const Scalar half = 0.5 * density * b.frame.length;
  const Matrix3& r = b.frame.rotation;
  const Matrix3 turned = r.transpose() * Vector3(b.iy + b.iz, b.iy, b.iz).asDiagonal() * r;
  const Matrix3 inertia = (0.5 * half) * (turned + turned.transpose());
  Matrix m = Matrix::Zero(kDofs, kDofs);
  for (int a = 0; a < 2; ++a) {
    m.block<3, 3>(6 * a, 6 * a) = (half * b.area) * Matrix3::Identity();
    m.block<3, 3>(6 * a + 3, 6 * a + 3) = inertia;
  }
  return m;
}

StrainOperator Beam2Element::strain_operator(const Matrix& geometry,
                                             const NaturalPoint& point) const {
  const BeamData b = beam_data(geometry);
  StrainOperator op;
  op.b = strain_rows(interpolation(0.5 * (point.xi + 1.0), b)) *
         transformation(b.frame.rotation);
  op.detJ = 0.5 * b.frame.length;
  return op;
}

Vector Beam2Element::shape_functions(const NaturalPoint& point) const {
  Vector n(2);
  n << 0.5 * (1.0 - point.xi), 0.5 * (1.0 + point.xi);
  return n;
}

std::vector<NaturalPoint> Beam2Element::stress_evaluation_points(
    const IntegrationOptions&) const {
  std::vector<NaturalPoint> out;
  for (const auto& gp : gauss_legendre_line(2)) out.push_back(NaturalPoint{gp.xi, 0.0, 0.0});
  return out;
}

std::vector<IntegrationPoint> Beam2Element::integration_rule(const IntegrationOptions&) const {
  std::vector<IntegrationPoint> out;
  for (const auto& gp : gauss_legendre_line(2)) {
    IntegrationPoint ip;
    ip.point = NaturalPoint{gp.xi, 0.0, 0.0};
    ip.weight = gp.weight;
    out.push_back(ip);
  }
  return out;
}

Matrix Beam2Element::geometric_stiffness(const Matrix& geometry, const Matrix& d,
                                         const Vector& ue, Scalar stress_scale, Scalar,
                                         const IntegrationOptions&) const {
  const BeamData b = beam_data(geometry);
  check_material(d, b);
  const Matrix12 t = transformation(b.frame.rotation);
  const Vector12 local = t * element_dofs(ue);
  const Scalar axial = stress_scale * b.e * b.area * (local(6) - local(0)) / b.frame.length;
  return t.transpose() * (axial * local_geometric(b)) * t;
}

Vector Beam2Element::geometric_stiffness_derivative(const Matrix& geometry, const Matrix& d,
                                                    const Vector& phi, Scalar stress_scale,
                                                    Scalar, const IntegrationOptions&) const {
  const BeamData b = beam_data(geometry);
  check_material(d, b);
  const Matrix12 t = transformation(b.frame.rotation);
  const Vector12 p = t * element_dofs(phi);
  // phi^T K_G(u) phi = s (EA / L) (u_1 - u_0) (p^T G p), linear in the local
  // axial displacements u_0 (DOF 0) and u_1 (DOF 6).
  const Scalar factor =
      stress_scale * b.e * b.area / b.frame.length * p.dot(local_geometric(b) * p);
  Vector12 g_local = Vector12::Zero();
  g_local(0) = -factor;
  g_local(6) = factor;
  return t.transpose() * g_local;
}

Vector Beam2Element::boundary_traction(const Matrix&, int, const Vector3&, Scalar,
                                       const IntegrationOptions&) const {
  throw ConfigError("a beam has no faces to take a traction: load it with line loads "
                    "(force per unit length) or point loads");
}

Vector Beam2Element::line_load(const Matrix& geometry, const Vector3& q) const {
  const BeamData b = beam_data(geometry);
  const Matrix12 t = transformation(b.frame.rotation);
  return t.transpose() * local_line_load(b, b.frame.rotation * q);
}

BeamEndForces Beam2Element::end_forces(const Matrix& geometry, const Matrix& d,
                                       const Vector& ue, const Vector3& q) const {
  const BeamData b = beam_data(geometry);
  check_material(d, b);
  const Matrix12 t = transformation(b.frame.rotation);
  const Vector12 f =
      local_stiffness(b) * (t * element_dofs(ue)) - local_line_load(b, b.frame.rotation * q);
  BeamEndForces out;
  out.frame = b.frame;
  out.start = -f.head<6>();
  out.end = f.tail<6>();
  return out;
}

}  // namespace sparlab
