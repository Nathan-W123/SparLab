#include "sparlab/elements/Quad4.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/elements/Quadrature.hpp"

#include <sstream>

namespace sparlab {

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

StrainOperator Quad4Element::strain_operator(
    const Eigen::Matrix<Scalar, 2, Eigen::Dynamic>& coords,
    const NaturalPoint& point) const {
  if (coords.cols() != 4) {
    std::ostringstream os;
    os << "Quad4 expects 4 nodal coordinate columns, received " << coords.cols();
    throw MeshError(os.str());
  }

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

  StrainOperator op;
  op.detJ = det;
  op.b.setZero(kVoigt, num_dofs());
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

Matrix Quad4Element::stiffness(const Eigen::Matrix<Scalar, 2, Eigen::Dynamic>& coords,
                               const Matrix3& d, Scalar thickness,
                               const IntegrationOptions& opts) const {
  if (!(thickness > 0.0)) {
    std::ostringstream os;
    os << "element thickness must be positive (got " << thickness << " m)";
    throw ConfigError(os.str());
  }
  Matrix ke = Matrix::Zero(num_dofs(), num_dofs());
  for (const auto& gp : gauss_legendre_square(opts.stiffness_points)) {
    NaturalPoint p;
    p.xi = gp.xi;
    p.eta = gp.eta;
    const StrainOperator op = strain_operator(coords, p);
    ke.noalias() += (thickness * op.detJ * gp.weight) * (op.b.transpose() * d * op.b);
  }
  // Enforce exact symmetry: the integrand is symmetric, so any asymmetry is
  // pure round-off. Symmetrising keeps Cholesky and the eigen solver happy.
  return 0.5 * (ke + ke.transpose());
}

Matrix Quad4Element::consistent_mass(
    const Eigen::Matrix<Scalar, 2, Eigen::Dynamic>& coords, Scalar density,
    Scalar thickness, const IntegrationOptions& opts) const {
  if (!(thickness > 0.0)) {
    std::ostringstream os;
    os << "element thickness must be positive (got " << thickness << " m)";
    throw ConfigError(os.str());
  }
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

std::array<int, 2> Quad4Element::edge_nodes(int local_edge) const {
  static const std::array<std::array<int, 2>, 4> table = {
      {{0, 1}, {1, 2}, {2, 3}, {3, 0}}};
  if (local_edge < 0 || local_edge > 3) {
    std::ostringstream os;
    os << "Quad4 local edge index " << local_edge << " is outside [0, 3]";
    throw MeshError(os.str());
  }
  return table[static_cast<std::size_t>(local_edge)];
}

Vector Quad4Element::edge_traction(const Eigen::Matrix<Scalar, 2, Eigen::Dynamic>& coords,
                                  int local_edge, const Vector2& traction,
                                  Scalar thickness,
                                  const IntegrationOptions& opts) const {
  const std::array<int, 2> en = edge_nodes(local_edge);
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

std::unique_ptr<Element> make_element(ElementType type) {
  switch (type) {
    case ElementType::Quad4: return std::make_unique<Quad4Element>();
  }
  throw ConfigError("no element implementation registered for the requested type");
}

}  // namespace sparlab
