#include "sparlab/elements/Shell4.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/elements/Quad4.hpp"
#include "sparlab/elements/Quadrature.hpp"

#include <Eigen/Geometry>
#include <Eigen/LU>

#include <array>
#include <cmath>
#include <sstream>

namespace sparlab {
namespace {

constexpr int kNodes = 4;
constexpr int kDofs = 24;
using Row = Eigen::Matrix<Scalar, 1, kDofs>;
using Operator3 = Eigen::Matrix<Scalar, 3, kDofs>;
using Operator6 = Eigen::Matrix<Scalar, 6, kDofs>;
using Matrix6 = Eigen::Matrix<Scalar, 6, 6>;

/// Natural coordinates of the corners, in node order.
constexpr std::array<std::array<Scalar, 2>, 4> kCorners{{{-1.0, -1.0}, {1.0, -1.0},
                                                         {1.0, 1.0}, {-1.0, 1.0}}};

struct ShellGeometry {
  Eigen::Matrix<Scalar, 3, 4> x;  ///< nodal coordinates
  Eigen::Matrix<Scalar, 3, 4> v;  ///< nodal directors, unit
};

Vector3 surface_normal(const Eigen::Matrix<Scalar, 3, 4>& x, Scalar r, Scalar s) {
  const Eigen::Matrix<Scalar, 4, 2> dn = quad4_shape_gradients_natural(r, s);
  const Vector3 gr = x * dn.col(0);
  const Vector3 gs = x * dn.col(1);
  return gr.cross(gs);
}

ShellGeometry shell_geometry(const Matrix& geometry) {
  if (geometry.cols() != kNodes || (geometry.rows() != 6 && geometry.rows() != 3)) {
    std::ostringstream os;
    os << "Shell4 expects a 6 x 4 geometry matrix (coordinates, then directors) or a 3 x 4 "
          "coordinate matrix, received "
       << geometry.rows() << " x " << geometry.cols();
    throw MeshError(os.str());
  }
  ShellGeometry g;
  g.x = geometry.topRows(3);
  for (int k = 0; k < kNodes; ++k) {
    const Vector3 v = geometry.rows() == 6
                          ? Vector3(geometry.col(k).tail(3))
                          : surface_normal(g.x, kCorners[static_cast<std::size_t>(k)][0],
                                           kCorners[static_cast<std::size_t>(k)][1]);
    const Scalar len = v.norm();
    if (!(len > 0.0) || !std::isfinite(len)) {
      std::ostringstream os;
      os << "Shell4 node " << k << " has no director: the element is collapsed at that "
            "corner (two edges parallel)";
      throw MeshError(os.str());
    }
    g.v.col(k) = v / len;
  }
  return g;
}

/// Psi with Psi theta = theta x V.
Matrix3 cross_with(const Vector3& v) {
  Matrix3 psi;
  psi << 0.0, v.z(), -v.y(),
         -v.z(), 0.0, v.x(),
         v.y(), -v.x(), 0.0;
  return psi;
}

/// e1 of the local frame for the normal e3: the projection of global x, or of
/// global z when x lies within 0.1 degree of e3.
Vector3 local_e1(const Vector3& e3) {
  static const Scalar kParallel = std::sin(0.1 * 3.14159265358979323846 / 180.0);
  Vector3 e1 = Vector3::UnitX() - e3.x() * e3;
  if (e1.norm() < kParallel) e1 = Vector3::UnitZ() - e3.z() * e3;
  return e1.normalized();
}

/// The quantities of one point (r, s, zeta).
struct ShellPoint {
  Eigen::Vector4d n;
  Eigen::Matrix<Scalar, 4, 2> dn;
  Vector3 gr, gs, gz;
  Scalar det = 0.0;
  Matrix3 frame;  ///< columns e1, e2, e3
  Matrix3 t;      ///< T = E^T G^-T: t(a, i) = e_a . g^i
};

ShellPoint evaluate(const ShellGeometry& g, Scalar thickness, Scalar r, Scalar s, Scalar z) {
  ShellPoint p;
  p.n = quad4_shape_functions(r, s);
  p.dn = quad4_shape_gradients_natural(r, s);
  const Scalar h = 0.5 * thickness;
  const Eigen::Matrix<Scalar, 3, 4> xz = g.x + (z * h) * g.v;
  p.gr = xz * p.dn.col(0);
  p.gs = xz * p.dn.col(1);
  p.gz = h * (g.v * p.n);
  Matrix3 base;
  base.col(0) = p.gr;
  base.col(1) = p.gs;
  base.col(2) = p.gz;
  p.det = base.determinant();
  if (!(p.det > 0.0)) {
    std::ostringstream os;
    os << "Shell4 volume Jacobian is " << p.det << " m^3 at (r, s, zeta) = (" << r << ", "
       << s << ", " << z << "): the element is folded, collapsed, or its directors point "
       << "against its normal";
    throw MeshError(os.str());
  }
  const Vector3 e3 = p.gz.normalized();
  const Vector3 e1 = local_e1(e3);
  p.frame.col(0) = e1;
  p.frame.col(1) = e3.cross(e1);
  p.frame.col(2) = e3;
  p.t = p.frame.transpose() * base.inverse().transpose();
  return p;
}

/// Adds factor * d(a . u_,j)/dq to `row`: j = 0 (r), 1 (s), 2 (zeta).
void add_row(Row& row, const ShellGeometry& g, const ShellPoint& p, Scalar thickness,
             Scalar z, const Vector3& a, int j, Scalar factor) {
  const Scalar h = 0.5 * thickness;
  for (int k = 0; k < kNodes; ++k) {
    const Vector3 va = g.v.col(k).cross(a);  // a . (theta x V) = (V x a) . theta
    if (j < 2) {
      const Scalar dnk = p.dn(k, j);
      row.segment<3>(6 * k) += (factor * dnk) * a.transpose();
      row.segment<3>(6 * k + 3) += (factor * z * h * dnk) * va.transpose();
    } else {
      row.segment<3>(6 * k + 3) += (factor * h * p.n(k)) * va.transpose();
    }
  }
}

/// The covariant strain rows at a point: rr, ss, rs, rzeta, szeta (tensor
/// components).
struct CovariantRows {
  Row rr, ss, rs, rz, sz;
};

CovariantRows covariant_rows(const ShellGeometry& g, const ShellPoint& p, Scalar thickness,
                             Scalar z) {
  CovariantRows c;
  c.rr.setZero();
  c.ss.setZero();
  c.rs.setZero();
  c.rz.setZero();
  c.sz.setZero();
  add_row(c.rr, g, p, thickness, z, p.gr, 0, 1.0);
  add_row(c.ss, g, p, thickness, z, p.gs, 1, 1.0);
  add_row(c.rs, g, p, thickness, z, p.gr, 1, 0.5);
  add_row(c.rs, g, p, thickness, z, p.gs, 0, 0.5);
  add_row(c.rz, g, p, thickness, z, p.gr, 2, 0.5);
  add_row(c.rz, g, p, thickness, z, p.gz, 0, 0.5);
  add_row(c.sz, g, p, thickness, z, p.gs, 2, 0.5);
  add_row(c.sz, g, p, thickness, z, p.gz, 1, 0.5);
  return c;
}

/// The local strain operator (11, 22, 33, 12, 23, 13; 33 zero) at a point,
/// its transverse shears interpolated from the tying points (MITC4).
Operator6 local_operator(const ShellGeometry& g, const ShellPoint& p, Scalar thickness,
                         Scalar r, Scalar s, Scalar z) {
  CovariantRows c = covariant_rows(g, p, thickness, z);
  // Tying points: A (0, 1), C (0, -1) for r-zeta; D (1, 0), B (-1, 0) for s-zeta.
  const CovariantRows a = covariant_rows(g, evaluate(g, thickness, 0.0, 1.0, z), thickness, z);
  const CovariantRows cc =
      covariant_rows(g, evaluate(g, thickness, 0.0, -1.0, z), thickness, z);
  const CovariantRows d = covariant_rows(g, evaluate(g, thickness, 1.0, 0.0, z), thickness, z);
  const CovariantRows b =
      covariant_rows(g, evaluate(g, thickness, -1.0, 0.0, z), thickness, z);
  c.rz = 0.5 * (1.0 + s) * a.rz + 0.5 * (1.0 - s) * cc.rz;
  c.sz = 0.5 * (1.0 + r) * d.sz + 0.5 * (1.0 - r) * b.sz;

  // eps~(i, j) as rows; (2, 2) not used.
  const std::array<std::array<const Row*, 3>, 3> e{{{&c.rr, &c.rs, &c.rz},
                                                    {&c.rs, &c.ss, &c.sz},
                                                    {&c.rz, &c.sz, nullptr}}};
  const auto component = [&](int aa, int bb) {
    Row out = Row::Zero();
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        const Row* ij = e[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        if (ij == nullptr) continue;
        const Scalar coef = p.t(aa, i) * p.t(bb, j);
        if (coef != 0.0) out += coef * (*ij);
      }
    }
    return out;
  };
  Operator6 op = Operator6::Zero();
  op.row(0) = component(0, 0);
  op.row(1) = component(1, 1);
  op.row(3) = 2.0 * component(0, 1);
  op.row(4) = 2.0 * component(1, 2);
  op.row(5) = 2.0 * component(0, 2);
  return op;
}

Matrix3 plane_stress(const Matrix& d) {
  if (d.rows() != 3 || d.cols() != 3) {
    std::ostringstream os;
    os << "Shell4 expects the 3 x 3 plane-stress constitutive matrix, received " << d.rows()
       << " x " << d.cols();
    throw ModelError(os.str());
  }
  return d;
}

/// The local 6 x 6 material: plane stress in (11, 22, 12), k G in shear.
Matrix6 shell_material(const Matrix3& d) {
  Matrix6 c = Matrix6::Zero();
  const std::array<int, 3> map{0, 1, 3};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      c(map[static_cast<std::size_t>(i)], map[static_cast<std::size_t>(j)]) = d(i, j);
    }
  }
  c(4, 4) = Shell4Element::kShearFactor * d(2, 2);
  c(5, 5) = Shell4Element::kShearFactor * d(2, 2);
  return c;
}

void check_thickness(Scalar thickness) {
  if (!(thickness > 0.0)) {
    std::ostringstream os;
    os << "shell thickness must be positive (got " << thickness << " m)";
    throw ConfigError(os.str());
  }
}

/// The displacement gradient operators u_,r, u_,s, u_,zeta (3 x 24 each).
std::array<Operator3, 3> gradient_operators(const ShellGeometry& g, const ShellPoint& p,
                                            Scalar thickness, Scalar z) {
  const Scalar h = 0.5 * thickness;
  std::array<Operator3, 3> out;
  for (Operator3& o : out) o.setZero();
  for (int k = 0; k < kNodes; ++k) {
    const Matrix3 psi = cross_with(g.v.col(k));
    for (int j = 0; j < 2; ++j) {
      out[static_cast<std::size_t>(j)].block<3, 3>(0, 6 * k) =
          p.dn(k, j) * Matrix3::Identity();
      out[static_cast<std::size_t>(j)].block<3, 3>(0, 6 * k + 3) = (z * h * p.dn(k, j)) * psi;
    }
    out[2].block<3, 3>(0, 6 * k + 3) = (h * p.n(k)) * psi;
  }
  return out;
}

/// D_a = sum_i T(a, i) u_,i for the in-plane directions a = 0, 1.
std::array<Operator3, 2> in_plane_gradients(const ShellGeometry& g, const ShellPoint& p,
                                            Scalar thickness, Scalar z) {
  const std::array<Operator3, 3> grad = gradient_operators(g, p, thickness, z);
  std::array<Operator3, 2> out;
  for (int a = 0; a < 2; ++a) {
    out[static_cast<std::size_t>(a)] =
        p.t(a, 0) * grad[0] + p.t(a, 1) * grad[1] + p.t(a, 2) * grad[2];
  }
  return out;
}

/// The drilling constraint at (r, s) of the mid-surface: the rotation about
/// the normal n less the in-plane rotation of the mid-surface,
/// n . sum N_k theta_k - 1/2 n . sum_a (a^a x u_,a), with a^a the in-plane
/// dual basis of g_r, g_s. Both equal n . omega under a rigid rotation omega,
/// on any geometry, so the rigid-body motions stay free.
Row drilling_row(const ShellGeometry& g, Scalar r, Scalar s) {
  const Eigen::Vector4d nn = quad4_shape_functions(r, s);
  const Eigen::Matrix<Scalar, 4, 2> dn = quad4_shape_gradients_natural(r, s);
  const Vector3 gr = g.x * dn.col(0);
  const Vector3 gs = g.x * dn.col(1);
  const Vector3 n = gr.cross(gs).normalized();
  Eigen::Matrix2d metric;
  metric << gr.dot(gr), gr.dot(gs), gs.dot(gr), gs.dot(gs);
  const Eigen::Matrix2d inv = metric.inverse();
  const Vector3 ar = inv(0, 0) * gr + inv(0, 1) * gs;
  const Vector3 as = inv(1, 0) * gr + inv(1, 1) * gs;
  // n . (a x u_,a) = (n x a) . u_,a = u_,a . (n x a)
  const Vector3 nr = n.cross(ar);
  const Vector3 ns = n.cross(as);
  Row c = Row::Zero();
  for (int k = 0; k < kNodes; ++k) {
    c.segment<3>(6 * k + 3) = nn(k) * n.transpose();
    c.segment<3>(6 * k) = -0.5 * (dn(k, 0) * nr + dn(k, 1) * ns).transpose();
  }
  return c;
}

Scalar plane_von_mises(Scalar s11, Scalar s22, Scalar s12) {
  return std::sqrt(std::max(0.0, s11 * s11 - s11 * s22 + s22 * s22 + 3.0 * s12 * s12));
}

}  // namespace

Shell4Element::Shell4Element(Scalar drilling_factor) : drilling_factor_(drilling_factor) {
  if (!(drilling_factor >= 0.0) || !std::isfinite(drilling_factor)) {
    std::ostringstream os;
    os << "the shell drilling stiffness factor must be finite and >= 0 (got "
       << drilling_factor << ")";
    throw ConfigError(os.str());
  }
}

Vector Shell4Element::shape_functions(const NaturalPoint& point) const {
  return quad4_shape_functions(point.xi, point.eta);
}

std::vector<NaturalPoint> Shell4Element::stress_evaluation_points(
    const IntegrationOptions& opts) const {
  std::vector<NaturalPoint> out;
  for (const IntegrationPoint& ip : integration_rule(opts)) out.push_back(ip.point);
  return out;
}

std::vector<IntegrationPoint> Shell4Element::integration_rule(const IntegrationOptions&) const {
  std::vector<IntegrationPoint> out;
  for (const auto& gp : gauss_legendre_square(2)) {
    for (const auto& zp : gauss_legendre_line(2)) {
      IntegrationPoint ip;
      ip.point.xi = gp.xi;
      ip.point.eta = gp.eta;
      ip.point.zeta = zp.xi;
      ip.weight = gp.weight * zp.weight;
      out.push_back(ip);
    }
  }
  return out;
}

StrainOperator Shell4Element::strain_operator(const Matrix&, const NaturalPoint&) const {
  throw ModelError("the shell strain operator depends on the thickness: call "
                   "Shell4Element::strain_operator_at");
}

StrainOperator Shell4Element::strain_operator_at(const Matrix& geometry, Scalar thickness,
                                                 const NaturalPoint& point) const {
  check_thickness(thickness);
  const ShellGeometry g = shell_geometry(geometry);
  const ShellPoint p = evaluate(g, thickness, point.xi, point.eta, point.zeta);
  StrainOperator op;
  op.b = local_operator(g, p, thickness, point.xi, point.eta, point.zeta);
  op.detJ = p.det;
  return op;
}

Matrix Shell4Element::stiffness(const Matrix& geometry, const Matrix& d_in, Scalar thickness,
                                const IntegrationOptions&) const {
  check_thickness(thickness);
  const ShellGeometry g = shell_geometry(geometry);
  const Matrix3 d = plane_stress(d_in);
  const Matrix6 c = shell_material(d);
  Matrix ke = Matrix::Zero(kDofs, kDofs);
  for (const auto& gp : gauss_legendre_square(2)) {
    for (const auto& zp : gauss_legendre_line(2)) {
      const ShellPoint p = evaluate(g, thickness, gp.xi, gp.eta, zp.xi);
      const Operator6 b = local_operator(g, p, thickness, gp.xi, gp.eta, zp.xi);
      ke.noalias() += (gp.weight * zp.weight * p.det) * (b.transpose() * c * b);
    }
  }
  // The rotation about the normal moves no point of the shell: a penalty ties
  // it to the in-plane rotation of the mid-surface (the class comment).
  const Scalar kd = drilling_factor_ * d(2, 2) * thickness;
  if (kd > 0.0) {
    for (const auto& gp : gauss_legendre_square(2)) {
      const Row drill = drilling_row(g, gp.xi, gp.eta);
      const Scalar area = surface_normal(g.x, gp.xi, gp.eta).norm();
      ke.noalias() += (kd * gp.weight * area) * (drill.transpose() * drill);
    }
  }
  return 0.5 * (ke + ke.transpose());
}

Matrix Shell4Element::consistent_mass(const Matrix& geometry, Scalar density, Scalar thickness,
                                      const IntegrationOptions& opts) const {
  check_thickness(thickness);
  const ShellGeometry g = shell_geometry(geometry);
  const Scalar h = 0.5 * thickness;
  Matrix me = Matrix::Zero(kDofs, kDofs);
  for (const auto& gp : gauss_legendre_square(opts.mass_points)) {
    for (const auto& zp : gauss_legendre_line(3)) {
      const ShellPoint p = evaluate(g, thickness, gp.xi, gp.eta, zp.xi);
      Operator3 nu = Operator3::Zero();
      for (int k = 0; k < kNodes; ++k) {
        nu.block<3, 3>(0, 6 * k) = p.n(k) * Matrix3::Identity();
        nu.block<3, 3>(0, 6 * k + 3) = (zp.xi * h * p.n(k)) * cross_with(g.v.col(k));
      }
      me.noalias() += (density * gp.weight * zp.weight * p.det) * (nu.transpose() * nu);
    }
  }
  return 0.5 * (me + me.transpose());
}

Matrix Shell4Element::geometric_stiffness(const Matrix& geometry, const Matrix& d_in,
                                          const Vector& ue, Scalar stress_scale,
                                          Scalar thickness, const IntegrationOptions&) const {
  check_thickness(thickness);
  if (ue.size() != kDofs) throw ModelError("Shell4 geometric stiffness needs 24 displacements");
  const ShellGeometry g = shell_geometry(geometry);
  const Matrix6 c = shell_material(plane_stress(d_in));
  Matrix kg = Matrix::Zero(kDofs, kDofs);
  for (const auto& gp : gauss_legendre_square(2)) {
    for (const auto& zp : gauss_legendre_line(2)) {
      const ShellPoint p = evaluate(g, thickness, gp.xi, gp.eta, zp.xi);
      const Operator6 b = local_operator(g, p, thickness, gp.xi, gp.eta, zp.xi);
      const Eigen::Matrix<Scalar, 6, 1> sigma = stress_scale * (c * (b * ue));
      const std::array<Operator3, 2> da = in_plane_gradients(g, p, thickness, zp.xi);
      const Scalar w = gp.weight * zp.weight * p.det;
      kg.noalias() += (w * sigma(0)) * (da[0].transpose() * da[0]);
      kg.noalias() += (w * sigma(1)) * (da[1].transpose() * da[1]);
      kg.noalias() += (w * sigma(3)) * (da[0].transpose() * da[1] + da[1].transpose() * da[0]);
    }
  }
  return 0.5 * (kg + kg.transpose());
}

Vector Shell4Element::geometric_stiffness_derivative(const Matrix& geometry, const Matrix& d_in,
                                                     const Vector& phi, Scalar stress_scale,
                                                     Scalar thickness,
                                                     const IntegrationOptions&) const {
  check_thickness(thickness);
  if (phi.size() != kDofs) throw ModelError("Shell4 needs a 24-component mode");
  const ShellGeometry g = shell_geometry(geometry);
  const Matrix3 d = plane_stress(d_in);
  Vector out = Vector::Zero(kDofs);
  for (const auto& gp : gauss_legendre_square(2)) {
    for (const auto& zp : gauss_legendre_line(2)) {
      const ShellPoint p = evaluate(g, thickness, gp.xi, gp.eta, zp.xi);
      const Operator6 b = local_operator(g, p, thickness, gp.xi, gp.eta, zp.xi);
      const std::array<Operator3, 2> da = in_plane_gradients(g, p, thickness, zp.xi);
      const Vector3 g0 = da[0] * phi;
      const Vector3 g1 = da[1] * phi;
      const Vector3 hat(g0.dot(g0), g1.dot(g1), 2.0 * g0.dot(g1));
      Operator3 b_in;
      b_in.row(0) = b.row(0);
      b_in.row(1) = b.row(1);
      b_in.row(2) = b.row(3);
      out.noalias() += (stress_scale * gp.weight * zp.weight * p.det) *
                       (b_in.transpose() * (d * hat));
    }
  }
  return out;
}

Vector Shell4Element::boundary_traction(const Matrix& geometry, int local_face,
                                        const Vector3& traction, Scalar thickness,
                                        const IntegrationOptions& opts) const {
  check_thickness(thickness);
  const std::vector<int>& edge = face_nodes(local_face);
  const Vector3 xa = geometry.block<3, 1>(0, edge[0]);
  const Vector3 xb = geometry.block<3, 1>(0, edge[1]);
  const Scalar half_length = 0.5 * (xb - xa).norm();
  Vector f = Vector::Zero(kDofs);
  for (const auto& gp : gauss_legendre_line(opts.edge_points)) {
    const Scalar na = 0.5 * (1.0 - gp.xi);
    const Scalar nb = 0.5 * (1.0 + gp.xi);
    const Scalar w = gp.weight * half_length * thickness;
    f.segment<3>(6 * edge[0]) += (w * na) * traction;
    f.segment<3>(6 * edge[1]) += (w * nb) * traction;
  }
  return f;
}

Vector Shell4Element::pressure_load(const Matrix& geometry, Scalar pressure,
                                    const IntegrationOptions&) const {
  const ShellGeometry g = shell_geometry(geometry);
  Vector f = Vector::Zero(kDofs);
  for (const auto& gp : gauss_legendre_square(2)) {
    const Eigen::Vector4d n = quad4_shape_functions(gp.xi, gp.eta);
    const Vector3 area_vector = surface_normal(g.x, gp.xi, gp.eta);
    for (int k = 0; k < kNodes; ++k) {
      f.segment<3>(6 * k) -= (pressure * gp.weight * n(k)) * area_vector;
    }
  }
  return f;
}

ShellResultants Shell4Element::resultants(const Matrix& geometry, const Matrix& d_in,
                                          const Vector& ue, Scalar thickness, Scalar r,
                                          Scalar s) const {
  check_thickness(thickness);
  if (ue.size() != kDofs) throw ModelError("Shell4 resultants need 24 displacements");
  const ShellGeometry g = shell_geometry(geometry);
  const Matrix6 c = shell_material(plane_stress(d_in));
  const Scalar h = 0.5 * thickness;
  const auto stress_at = [&](Scalar z) -> Eigen::Matrix<Scalar, 6, 1> {
    const ShellPoint p = evaluate(g, thickness, r, s, z);
    return c * (local_operator(g, p, thickness, r, s, z) * ue);
  };
  ShellResultants out;
  const ShellPoint mid = evaluate(g, thickness, r, s, 0.0);
  out.e1 = mid.frame.col(0);
  out.e2 = mid.frame.col(1);
  out.e3 = mid.frame.col(2);
  for (const auto& zp : gauss_legendre_line(2)) {
    const Eigen::Matrix<Scalar, 6, 1> sigma = stress_at(zp.xi);
    const Vector3 in_plane(sigma(0), sigma(1), sigma(3));
    out.membrane += (h * zp.weight) * in_plane;
    out.moment += (h * h * zp.weight * zp.xi) * in_plane;
    out.shear += (h * zp.weight) * Eigen::Vector2d(sigma(5), sigma(4));
  }
  const Eigen::Matrix<Scalar, 6, 1> top = stress_at(1.0);
  const Eigen::Matrix<Scalar, 6, 1> bottom = stress_at(-1.0);
  const Eigen::Matrix<Scalar, 6, 1> middle = stress_at(0.0);
  out.stress_top = Vector3(top(0), top(1), top(3));
  out.stress_bottom = Vector3(bottom(0), bottom(1), bottom(3));
  out.von_mises_top = plane_von_mises(top(0), top(1), top(3));
  out.von_mises_bottom = plane_von_mises(bottom(0), bottom(1), bottom(3));
  // The transverse shear stress of a homogeneous section is parabolic, zero
  // on the faces and 3 Q / (2 t) on the mid-surface.
  const Scalar tau13 = 1.5 * out.shear(0) / thickness;
  const Scalar tau23 = 1.5 * out.shear(1) / thickness;
  out.von_mises_mid = std::sqrt(std::max(
      0.0, middle(0) * middle(0) - middle(0) * middle(1) + middle(1) * middle(1) +
               3.0 * (middle(3) * middle(3) + tau13 * tau13 + tau23 * tau23)));
  return out;
}

Scalar Shell4Element::area(const Matrix& geometry) {
  const ShellGeometry g = shell_geometry(geometry);
  Scalar a = 0.0;
  for (const auto& gp : gauss_legendre_square(2)) {
    a += gp.weight * surface_normal(g.x, gp.xi, gp.eta).norm();
  }
  return a;
}

Vector3 Shell4Element::normal(const Matrix& geometry, Scalar r, Scalar s) {
  const ShellGeometry g = shell_geometry(geometry);
  return surface_normal(g.x, r, s).normalized();
}

}  // namespace sparlab
