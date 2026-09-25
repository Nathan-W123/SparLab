#include "sparlab/elements/Quad4.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/elements/Quadrature.hpp"

#include <sstream>

namespace sparlab {
namespace {

/// Nodal coordinates in the fixed-row layout the kernels below are written
/// for. The kernels keep their original fixed-size Eigen types deliberately:
/// the products are then evaluated by exactly the same code paths as before
/// the interface became dimension-generic, which is what keeps every 2-D
/// result bit-for-bit reproducible.
using Quad4Coords = Eigen::Matrix<Scalar, 2, Eigen::Dynamic>;

Quad4Coords quad4_coords(const Matrix& coords) {
  if (coords.rows() != 2 || coords.cols() != 4) {
    std::ostringstream os;
    os << "Quad4 expects a 2 x 4 nodal coordinate matrix, received " << coords.rows()
       << " x " << coords.cols();
    throw MeshError(os.str());
  }
  return coords;
}

Matrix3 quad4_constitutive(const Matrix& d) {
  if (d.rows() != 3 || d.cols() != 3) {
    std::ostringstream os;
    os << "Quad4 expects a 3 x 3 constitutive matrix, received " << d.rows() << " x "
       << d.cols();
    throw ModelError(os.str());
  }
  return d;
}

/// Strain operator in the element's native fixed-size layout.
struct Quad4Kernel {
  Eigen::Matrix<Scalar, 3, Eigen::Dynamic> b;  ///< 3 x 8
  Scalar detJ = 0.0;
};

Quad4Kernel quad4_kernel(const Quad4Coords& coords, const NaturalPoint& point) {
  const Eigen::Matrix<Scalar, 4, 2> dn_dxi =
      quad4_shape_gradients_natural(point.xi, point.eta);

  // J = sum_a x_a (dN_a/dxi)^T  =>  J(i,j) = dx_i/dxi_j
  const Matrix2 jac = coords * dn_dxi;
  // The 2x2 determinant and inverse are written out: the assembly loop runs
  // this millions of times and neither needs Eigen's LU machinery.
  const Scalar det = jac(0, 0) * jac(1, 1) - jac(0, 1) * jac(1, 0);
  if (!(det > 0.0)) {
    std::ostringstream os;
    os << "Quad4 Jacobian determinant is " << det << " m^2 at (xi, eta) = (" << point.xi
       << ", " << point.eta
       << "); the element is inverted, collapsed or its nodes are not counter-clockwise";
    throw MeshError(os.str());
  }

  // dN/dx = dN/dxi * J^{-1}  (row a holds grad N_a).
  Matrix2 jinv;
  jinv(0, 0) = jac(1, 1) / det;
  jinv(0, 1) = -jac(0, 1) / det;
  jinv(1, 0) = -jac(1, 0) / det;
  jinv(1, 1) = jac(0, 0) / det;
  const Eigen::Matrix<Scalar, 4, 2> dn_dx = dn_dxi * jinv;

  Quad4Kernel op;
  op.detJ = det;
  op.b.setZero(3, 8);
  for (int a = 0; a < 4; ++a) {
    const Scalar dx = dn_dx(a, 0);
    const Scalar dy = dn_dx(a, 1);
    op.b(0, 2 * a + 0) = dx;   // eps_xx = du/dx
    op.b(1, 2 * a + 1) = dy;   // eps_yy = dv/dy
    op.b(2, 2 * a + 0) = dy;   // gamma_xy = du/dy + dv/dx
    op.b(2, 2 * a + 1) = dx;
  }
  return op;
}

}  // namespace

Eigen::Vector4d quad4_shape_functions(Scalar xi, Scalar eta) {
  Eigen::Vector4d n;
  n << 0.25 * (1.0 - xi) * (1.0 - eta),
       0.25 * (1.0 + xi) * (1.0 - eta),
       0.25 * (1.0 + xi) * (1.0 + eta),
       0.25 * (1.0 - xi) * (1.0 + eta);
  return n;
}

Eigen::Matrix<Scalar, 4, 2> quad4_shape_gradients_natural(Scalar xi, Scalar eta) {
  Eigen::Matrix<Scalar, 4, 2> g;
  // dN/dxi                          dN/deta
  g(0, 0) = -0.25 * (1.0 - eta);  g(0, 1) = -0.25 * (1.0 - xi);
  g(1, 0) =  0.25 * (1.0 - eta);  g(1, 1) = -0.25 * (1.0 + xi);
  g(2, 0) =  0.25 * (1.0 + eta);  g(2, 1) =  0.25 * (1.0 + xi);
  g(3, 0) = -0.25 * (1.0 + eta);  g(3, 1) =  0.25 * (1.0 - xi);
  return g;
}

Vector Quad4Element::shape_functions(const NaturalPoint& point) const {
  return quad4_shape_functions(point.xi, point.eta);
}

StrainOperator Quad4Element::strain_operator(const Matrix& coords,
                                             const NaturalPoint& point) const {
  const Quad4Kernel kernel = quad4_kernel(quad4_coords(coords), point);
  StrainOperator op;
  op.b = kernel.b;
  op.detJ = kernel.detJ;
  return op;
}

Matrix Quad4Element::stiffness(const Matrix& coords_in, const Matrix& d_in,
                               Scalar thickness, const IntegrationOptions& opts) const {
  if (!(thickness > 0.0)) {
    std::ostringstream os;
    os << "element thickness must be positive (got " << thickness << " m)";
    throw ConfigError(os.str());
  }
  const Quad4Coords coords = quad4_coords(coords_in);
  const Matrix3 d = quad4_constitutive(d_in);
  Matrix ke = Matrix::Zero(num_dofs(), num_dofs());
  for (const auto& gp : gauss_legendre_square(opts.stiffness_points)) {
    NaturalPoint p;
    p.xi = gp.xi;
    p.eta = gp.eta;
    const Quad4Kernel op = quad4_kernel(coords, p);
    ke.noalias() += (thickness * op.detJ * gp.weight) * (op.b.transpose() * d * op.b);
  }
  // Enforce exact symmetry: the integrand is symmetric, so any asymmetry is
  // pure round-off. Symmetrising keeps Cholesky and the eigen solver happy.
  return 0.5 * (ke + ke.transpose());
}

Matrix Quad4Element::consistent_mass(const Matrix& coords_in, Scalar density,
                                     Scalar thickness,
                                     const IntegrationOptions& opts) const {
  if (!(thickness > 0.0)) {
    std::ostringstream os;
    os << "element thickness must be positive (got " << thickness << " m)";
    throw ConfigError(os.str());
  }
  const Quad4Coords coords = quad4_coords(coords_in);
  Matrix me = Matrix::Zero(num_dofs(), num_dofs());
  for (const auto& gp : gauss_legendre_square(opts.mass_points)) {
    NaturalPoint p;
    p.xi = gp.xi;
    p.eta = gp.eta;
    const Eigen::Vector4d n = quad4_shape_functions(p.xi, p.eta);
    const Eigen::Matrix<Scalar, 4, 2> dn_dxi =
        quad4_shape_gradients_natural(p.xi, p.eta);
    const Matrix2 jac = coords * dn_dxi;
    const Scalar det = jac(0, 0) * jac(1, 1) - jac(0, 1) * jac(1, 0);
    if (!(det > 0.0)) {
      std::ostringstream os;
      os << "Quad4 Jacobian determinant is " << det
         << " m^2 during mass integration; element geometry is invalid";
      throw MeshError(os.str());
    }
    // N_mat is 2 x 8 with the node-major DOF ordering.
    Eigen::Matrix<Scalar, 2, 8> nmat = Eigen::Matrix<Scalar, 2, 8>::Zero();
    for (int a = 0; a < 4; ++a) {
      nmat(0, 2 * a + 0) = n(a);
      nmat(1, 2 * a + 1) = n(a);
    }
    me.noalias() +=
        (density * thickness * det * gp.weight) * (nmat.transpose() * nmat);
  }
  return 0.5 * (me + me.transpose());
}

std::vector<NaturalPoint> Quad4Element::stress_evaluation_points(
    const IntegrationOptions& opts) const {
  std::vector<NaturalPoint> points;
  for (const auto& gp : gauss_legendre_square(opts.stiffness_points)) {
    NaturalPoint p;
    p.xi = gp.xi;
    p.eta = gp.eta;
    points.push_back(p);
  }
  return points;
}

Vector Quad4Element::boundary_traction(const Matrix& coords_in, int local_face,
                                       const Vector3& traction, Scalar thickness,
                                       const IntegrationOptions& opts) const {
  const std::vector<int>& en = face_nodes(local_face);
  if (traction.z() != 0.0) {
    std::ostringstream os;
    os << "a plane element cannot carry an out-of-plane traction (t_z = " << traction.z()
       << " Pa); use a 3-D mesh or drop the z component";
    throw ConfigError(os.str());
  }
  const Quad4Coords coords = quad4_coords(coords_in);
  const Vector2 xa = coords.col(en[0]);
  const Vector2 xb = coords.col(en[1]);

  Vector fe = Vector::Zero(num_dofs());
  // Parametrise the edge with s in [-1, 1]; for a straight Q4 edge the
  // Jacobian is the constant half-length |x_b - x_a| / 2.
  const Scalar half_length = 0.5 * (xb - xa).norm();
  if (!(half_length > 0.0)) throw MeshError("zero-length element edge in traction load");

  for (const auto& gp : gauss_legendre_line(opts.edge_points)) {
    const Scalar na = 0.5 * (1.0 - gp.xi);
    const Scalar nb = 0.5 * (1.0 + gp.xi);
    const Scalar scale = thickness * half_length * gp.weight;
    fe(2 * en[0] + 0) += scale * na * traction.x();
    fe(2 * en[0] + 1) += scale * na * traction.y();
    fe(2 * en[1] + 0) += scale * nb * traction.x();
    fe(2 * en[1] + 1) += scale * nb * traction.y();
  }
  return fe;
}

std::vector<IntegrationPoint> Quad4Element::integration_rule(
    const IntegrationOptions& opts) const {
  std::vector<IntegrationPoint> rule;
  for (const auto& gp : gauss_legendre_square(opts.stiffness_points)) {
    IntegrationPoint ip;
    ip.point.xi = gp.xi;
    ip.point.eta = gp.eta;
    ip.weight = gp.weight;
    rule.push_back(ip);
  }
  return rule;
}

}  // namespace sparlab
