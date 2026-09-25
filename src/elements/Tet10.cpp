#include "sparlab/elements/Tet10.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/elements/Quadrature.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <limits>
#include <sstream>

namespace sparlab {
namespace {

using Tet10Coords = Eigen::Matrix<Scalar, 3, 10>;
using Tet10Gradients = Eigen::Matrix<Scalar, 10, 3>;
using Tet10StrainMatrix = Eigen::Matrix<Scalar, 6, 30>;

/// Natural coordinates of the ten nodes (corners, then the edge midpoints of
/// 0-1, 1-2, 2-0, 0-3, 1-3, 2-3).
constexpr Scalar kNodes[10][3] = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                                  {0.0, 0.0, 1.0}, {0.5, 0.0, 0.0}, {0.5, 0.5, 0.0},
                                  {0.0, 0.5, 0.0}, {0.0, 0.0, 0.5}, {0.5, 0.0, 0.5},
                                  {0.0, 0.5, 0.5}};

/// Corner pairs of the six edge nodes 4..9.
constexpr int kEdges[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};

Tet10Coords tet10_coords(const Matrix& coords) {
  if (coords.rows() != 3 || coords.cols() != 10) {
    std::ostringstream os;
    os << "Tet10 expects a 3 x 10 nodal coordinate matrix, received " << coords.rows()
       << " x " << coords.cols();
    throw MeshError(os.str());
  }
  return coords;
}

void require_unit_thickness(Scalar thickness) {
  if (thickness != 1.0) {
    std::ostringstream os;
    os << "a solid (Tet10) element has no thickness; received " << thickness
       << " m. Leave model.thickness at its default of 1 for a 3-D mesh";
    throw ConfigError(os.str());
  }
}

/// Isoparametric map at one natural point: Cartesian gradients and det J.
struct Tet10Mapping {
  Tet10Gradients dn_dx;
  Scalar det = 0.0;
};

Tet10Mapping tet10_mapping(const Tet10Coords& x, Scalar xi, Scalar eta, Scalar zeta,
                           const char* context) {
  const Tet10Gradients dn = tet10_shape_gradients_natural(xi, eta, zeta);
  // J(i, j) = dx_i / dxi_j.
  const Matrix3 jac = x * dn;
  const Scalar det = jac.determinant();
  if (!(det > 0.0)) {
    std::ostringstream os;
    os << "Tet10 Jacobian determinant is " << det << " m^3 at natural point (" << xi
       << ", " << eta << ", " << zeta << ")" << (context[0] != '\0' ? " " : "") << context
       << "; the element is inverted, or an edge node sits so far off the straight edge "
          "that the curved element folds. Check that the corners satisfy "
          "(x1 - x0) x (x2 - x0) . (x3 - x0) > 0 and that edge nodes lie near their "
          "edge midpoints";
    throw MeshError(os.str());
  }
  Tet10Mapping map;
  map.det = det;
  // dN/dx_i = sum_j dN/dxi_j dxi_j/dx_i, i.e. dn * J^{-1}.
  map.dn_dx = dn * jac.inverse();
  return map;
}

Tet10StrainMatrix tet10_strain_matrix(const Tet10Gradients& g) {
  Tet10StrainMatrix b = Tet10StrainMatrix::Zero();
  for (int a = 0; a < 10; ++a) {
    const Scalar dx = g(a, 0);
    const Scalar dy = g(a, 1);
    const Scalar dz = g(a, 2);
    const int c = 3 * a;
    b(0, c + 0) = dx;  // eps_xx
    b(1, c + 1) = dy;  // eps_yy
    b(2, c + 2) = dz;  // eps_zz
    b(3, c + 0) = dy;  // gamma_xy
    b(3, c + 1) = dx;
    b(4, c + 1) = dz;  // gamma_yz
    b(4, c + 2) = dy;
    b(5, c + 0) = dz;  // gamma_zx
    b(5, c + 2) = dx;
  }
  return b;
}

/// Six-node triangle shape functions and their (r, s) derivatives, in the
/// face node order corners 0-2, then the edges 0-1, 1-2, 2-0.
void tri6_shape(Scalar r, Scalar s, Eigen::Matrix<Scalar, 6, 1>& n,
                Eigen::Matrix<Scalar, 6, 2>& dn) {
  const Scalar l0 = 1.0 - r - s;
  n << l0 * (2.0 * l0 - 1.0), r * (2.0 * r - 1.0), s * (2.0 * s - 1.0), 4.0 * l0 * r,
      4.0 * r * s, 4.0 * s * l0;
  // dL0/dr = dL0/ds = -1.
  dn(0, 0) = -(4.0 * l0 - 1.0);
  dn(0, 1) = -(4.0 * l0 - 1.0);
  dn(1, 0) = 4.0 * r - 1.0;
  dn(1, 1) = 0.0;
  dn(2, 0) = 0.0;
  dn(2, 1) = 4.0 * s - 1.0;
  dn(3, 0) = 4.0 * (l0 - r);
  dn(3, 1) = -4.0 * r;
  dn(4, 0) = 4.0 * s;
  dn(4, 1) = 4.0 * r;
  dn(5, 0) = -4.0 * s;
  dn(5, 1) = 4.0 * (l0 - s);
}

/// The 16-point collapsed Gauss rule on the face: degree 6, enough for the
/// quadratic shape functions times a quartic area element of a curved face.
constexpr int kFacePoints = 4;

}  // namespace

Eigen::Matrix<Scalar, 10, 1> tet10_shape_functions(Scalar xi, Scalar eta, Scalar zeta) {
  const Scalar l[4] = {1.0 - xi - eta - zeta, xi, eta, zeta};
  Eigen::Matrix<Scalar, 10, 1> n;
  for (int i = 0; i < 4; ++i) n(i) = l[i] * (2.0 * l[i] - 1.0);
  for (int k = 0; k < 6; ++k) n(4 + k) = 4.0 * l[kEdges[k][0]] * l[kEdges[k][1]];
  return n;
}

Eigen::Matrix<Scalar, 10, 3> tet10_shape_gradients_natural(Scalar xi, Scalar eta,
                                                           Scalar zeta) {
  const Scalar l[4] = {1.0 - xi - eta - zeta, xi, eta, zeta};
  // dL_i / d(xi, eta, zeta).
  static const Scalar dl[4][3] = {{-1.0, -1.0, -1.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                                  {0.0, 0.0, 1.0}};
  Eigen::Matrix<Scalar, 10, 3> g;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 3; ++j) g(i, j) = (4.0 * l[i] - 1.0) * dl[i][j];
  }
  for (int k = 0; k < 6; ++k) {
    const int a = kEdges[k][0];
    const int b = kEdges[k][1];
    for (int j = 0; j < 3; ++j) g(4 + k, j) = 4.0 * (l[a] * dl[b][j] + l[b] * dl[a][j]);
  }
  return g;
}

Scalar tet10_volume(const Matrix& coords_in, Scalar* min_det) {
  const Tet10Coords x = tet10_coords(coords_in);
  Scalar volume = 0.0;
  for (const QuadraturePoint3D& gp : collapsed_gauss_tetrahedron(3)) {
    const Matrix3 jac = x * tet10_shape_gradients_natural(gp.xi, gp.eta, gp.zeta);
    volume += gp.weight * jac.determinant();
  }
  if (min_det != nullptr) {
    Scalar lowest = std::numeric_limits<Scalar>::max();
    for (const QuadraturePoint3D& gp : tetrahedron_rule_4()) {
      const Matrix3 jac = x * tet10_shape_gradients_natural(gp.xi, gp.eta, gp.zeta);
      lowest = std::min(lowest, jac.determinant());
    }
    *min_det = lowest;
  }
  return volume;
}

Scalar tet10_jacobian_ratio(const Matrix& coords_in) {
  const Tet10Coords x = tet10_coords(coords_in);
  Scalar lowest = std::numeric_limits<Scalar>::max();
  Scalar highest = -std::numeric_limits<Scalar>::max();
  for (const auto& node : kNodes) {
    const Matrix3 jac = x * tet10_shape_gradients_natural(node[0], node[1], node[2]);
    const Scalar det = jac.determinant();
    lowest = std::min(lowest, det);
    highest = std::max(highest, det);
  }
  if (!(highest > 0.0)) return -1.0;
  return lowest / highest;
}

Scalar tri6_face_area(const Matrix& face) {
  if (face.rows() != 3 || face.cols() != 6) {
    throw MeshError("tri6_face_area expects a 3 x 6 coordinate matrix");
  }
  Scalar area = 0.0;
  Eigen::Matrix<Scalar, 6, 1> n;
  Eigen::Matrix<Scalar, 6, 2> dn;
  for (const QuadraturePoint2D& gp : collapsed_gauss_triangle(kFacePoints)) {
    tri6_shape(gp.xi, gp.eta, n, dn);
    const Vector3 xr = face * dn.col(0);
    const Vector3 xs = face * dn.col(1);
    area += gp.weight * xr.cross(xs).norm();
  }
  return area;
}

NaturalPoint Tet10Element::reference_centroid() const {
  NaturalPoint p;
  p.xi = 0.25;
  p.eta = 0.25;
  p.zeta = 0.25;
  return p;
}

Vector Tet10Element::shape_functions(const NaturalPoint& point) const {
  return tet10_shape_functions(point.xi, point.eta, point.zeta);
}

StrainOperator Tet10Element::strain_operator(const Matrix& coords,
                                             const NaturalPoint& point) const {
  const Tet10Mapping map =
      tet10_mapping(tet10_coords(coords), point.xi, point.eta, point.zeta, "");
  StrainOperator op;
  op.detJ = map.det;
  op.b = tet10_strain_matrix(map.dn_dx);
  return op;
}

Matrix Tet10Element::stiffness(const Matrix& coords_in, const Matrix& d_in, Scalar thickness,
                               const IntegrationOptions& /*opts*/) const {
  require_unit_thickness(thickness);
  if (d_in.rows() != 6 || d_in.cols() != 6) {
    std::ostringstream os;
    os << "Tet10 expects a 6 x 6 constitutive matrix, received " << d_in.rows() << " x "
       << d_in.cols();
    throw ModelError(os.str());
  }
  const Tet10Coords x = tet10_coords(coords_in);
  const Matrix6 d = d_in;
  Eigen::Matrix<Scalar, 30, 30> ke = Eigen::Matrix<Scalar, 30, 30>::Zero();
  for (const QuadraturePoint3D& gp : tetrahedron_rule_4()) {
    const Tet10Mapping map =
        tet10_mapping(x, gp.xi, gp.eta, gp.zeta, "during stiffness integration");
    const Tet10StrainMatrix b = tet10_strain_matrix(map.dn_dx);
    ke.noalias() += (map.det * gp.weight) * (b.transpose() * d * b);
  }
  return 0.5 * (ke + ke.transpose());
}

Matrix Tet10Element::consistent_mass(const Matrix& coords_in, Scalar density,
                                     Scalar thickness,
                                     const IntegrationOptions& /*opts*/) const {
  require_unit_thickness(thickness);
  const Tet10Coords x = tet10_coords(coords_in);
  // The scalar block int rho N N^T is shared by the three components.
  Eigen::Matrix<Scalar, 10, 10> block = Eigen::Matrix<Scalar, 10, 10>::Zero();
  for (const QuadraturePoint3D& gp : collapsed_gauss_tetrahedron(4)) {
    const Eigen::Matrix<Scalar, 10, 1> n = tet10_shape_functions(gp.xi, gp.eta, gp.zeta);
    const Tet10Mapping map =
        tet10_mapping(x, gp.xi, gp.eta, gp.zeta, "during mass integration");
    block.noalias() += (density * map.det * gp.weight) * (n * n.transpose());
  }
  Matrix me = Matrix::Zero(30, 30);
  for (int a = 0; a < 10; ++a) {
    for (int b = 0; b < 10; ++b) {
      const Scalar value = 0.5 * (block(a, b) + block(b, a));
      for (int k = 0; k < 3; ++k) me(3 * a + k, 3 * b + k) = value;
    }
  }
  return me;
}

std::vector<NaturalPoint> Tet10Element::stress_evaluation_points(
    const IntegrationOptions& /*opts*/) const {
  std::vector<NaturalPoint> points;
  for (const QuadraturePoint3D& gp : tetrahedron_rule_4()) {
    NaturalPoint p;
    p.xi = gp.xi;
    p.eta = gp.eta;
    p.zeta = gp.zeta;
    points.push_back(p);
  }
  return points;
}

std::vector<IntegrationPoint> Tet10Element::integration_rule(
    const IntegrationOptions& /*opts*/) const {
  std::vector<IntegrationPoint> rule;
  for (const QuadraturePoint3D& gp : tetrahedron_rule_4()) {
    IntegrationPoint ip;
    ip.point.xi = gp.xi;
    ip.point.eta = gp.eta;
    ip.point.zeta = gp.zeta;
    ip.weight = gp.weight;
    rule.push_back(ip);
  }
  return rule;
}

Vector Tet10Element::boundary_traction(const Matrix& coords_in, int local_face,
                                       const Vector3& traction, Scalar thickness,
                                       const IntegrationOptions& /*opts*/) const {
  require_unit_thickness(thickness);
  const Tet10Coords x = tet10_coords(coords_in);
  const std::vector<int>& fn = face_nodes(local_face);
  // The face is the 6-node triangle through its own nodes; on it the
  // element's shape functions reduce to the face's, and the others vanish.
  Eigen::Matrix<Scalar, 3, 6> xf;
  for (int a = 0; a < 6; ++a) xf.col(a) = x.col(fn[static_cast<std::size_t>(a)]);

  Vector fe = Vector::Zero(num_dofs());
  Scalar area = 0.0;
  Eigen::Matrix<Scalar, 6, 1> n;
  Eigen::Matrix<Scalar, 6, 2> dn;
  for (const QuadraturePoint2D& gp : collapsed_gauss_triangle(kFacePoints)) {
    tri6_shape(gp.xi, gp.eta, n, dn);
    const Scalar da = (xf * dn.col(0)).cross(xf * dn.col(1)).norm();
    area += gp.weight * da;
    for (int a = 0; a < 6; ++a) {
      const int node = fn[static_cast<std::size_t>(a)];
      const Scalar scale = gp.weight * da * n(a);
      fe(3 * node + 0) += scale * traction.x();
      fe(3 * node + 1) += scale * traction.y();
      fe(3 * node + 2) += scale * traction.z();
    }
  }
  if (!(area > 0.0)) throw MeshError("zero-area element face in traction load");
  return fe;
}

}  // namespace sparlab
