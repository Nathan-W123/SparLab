#include "sparlab/elements/Tet4.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <Eigen/Dense>

#include <sstream>

namespace sparlab {
namespace {

void check_coords(const Matrix& coords) {
  if (coords.rows() != 3 || coords.cols() != 4) {
    std::ostringstream os;
    os << "Tet4 expects a 3 x 4 nodal coordinate matrix, received " << coords.rows()
       << " x " << coords.cols();
    throw MeshError(os.str());
  }
}

void require_unit_thickness(Scalar thickness) {
  if (thickness != 1.0) {
    std::ostringstream os;
    os << "a solid (Tet4) element has no thickness; received " << thickness
       << " m. Leave model.thickness at its default of 1 for a 3-D mesh";
    throw ConfigError(os.str());
  }
}

Matrix strain_matrix(const Eigen::Matrix<Scalar, 4, 3>& g) {
  Matrix b = Matrix::Zero(6, 12);
  for (int a = 0; a < 4; ++a) {
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

}  // namespace

Eigen::Matrix<Scalar, 4, 3> tet4_shape_gradients(const Matrix& coords, Scalar* volume) {
  check_coords(coords);
  // J(i, j) = dx_i / dxi_j, columns are the edge vectors from node 0.
  Matrix3 jac;
  jac.col(0) = coords.col(1) - coords.col(0);
  jac.col(1) = coords.col(2) - coords.col(0);
  jac.col(2) = coords.col(3) - coords.col(0);
  const Scalar det = jac.determinant();
  if (!(det > 0.0)) {
    std::ostringstream os;
    os << "Tet4 signed volume is " << det / 6.0
       << " m^3; the element is inverted or flat, or its nodes are not ordered with "
          "(x1 - x0) x (x2 - x0) . (x3 - x0) > 0";
    throw MeshError(os.str());
  }
  // grad N_a = J^{-T} grad_xi N_a with grad_xi N = [-1 -1 -1; 1 0 0; 0 1 0; 0 0 1].
  const Matrix3 jinv = jac.inverse();
  Eigen::Matrix<Scalar, 4, 3> g;
  for (int j = 0; j < 3; ++j) {
    g(1 + j, 0) = jinv(j, 0);
    g(1 + j, 1) = jinv(j, 1);
    g(1 + j, 2) = jinv(j, 2);
  }
  g.row(0) = -(g.row(1) + g.row(2) + g.row(3));
  if (volume != nullptr) *volume = det / 6.0;
  return g;
}

NaturalPoint Tet4Element::reference_centroid() const {
  NaturalPoint p;
  p.xi = 0.25;
  p.eta = 0.25;
  p.zeta = 0.25;
  return p;
}

Vector Tet4Element::shape_functions(const NaturalPoint& point) const {
  Vector n(4);
  n << 1.0 - point.xi - point.eta - point.zeta, point.xi, point.eta, point.zeta;
  return n;
}

StrainOperator Tet4Element::strain_operator(const Matrix& coords,
                                            const NaturalPoint& /*point*/) const {
  Scalar volume = 0.0;
  const Eigen::Matrix<Scalar, 4, 3> g = tet4_shape_gradients(coords, &volume);
  StrainOperator op;
  op.b = strain_matrix(g);
  op.detJ = 6.0 * volume;
  return op;
}

Matrix Tet4Element::stiffness(const Matrix& coords, const Matrix& d, Scalar thickness,
                              const IntegrationOptions& /*opts*/) const {
  require_unit_thickness(thickness);
  if (d.rows() != 6 || d.cols() != 6) {
    std::ostringstream os;
    os << "Tet4 expects a 6 x 6 constitutive matrix, received " << d.rows() << " x "
       << d.cols();
    throw ModelError(os.str());
  }
  Scalar volume = 0.0;
  const Matrix b = strain_matrix(tet4_shape_gradients(coords, &volume));
  const Matrix ke = volume * (b.transpose() * d * b);
  return 0.5 * (ke + ke.transpose());
}

Matrix Tet4Element::consistent_mass(const Matrix& coords, Scalar density, Scalar thickness,
                                    const IntegrationOptions& /*opts*/) const {
  require_unit_thickness(thickness);
  Scalar volume = 0.0;
  tet4_shape_gradients(coords, &volume);
  const Scalar m = density * volume / 20.0;
  Matrix me = Matrix::Zero(12, 12);
  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      const Scalar value = (a == b ? 2.0 : 1.0) * m;
      for (int k = 0; k < 3; ++k) me(3 * a + k, 3 * b + k) = value;
    }
  }
  return me;
}

std::vector<NaturalPoint> Tet4Element::stress_evaluation_points(
    const IntegrationOptions& /*opts*/) const {
  return {reference_centroid()};
}

Vector Tet4Element::boundary_traction(const Matrix& coords, int local_face,
                                      const Vector3& traction, Scalar thickness,
                                      const IntegrationOptions& /*opts*/) const {
  require_unit_thickness(thickness);
  check_coords(coords);
  const std::vector<int>& fn = face_nodes(local_face);
  const Vector3 a = coords.col(fn[0]);
  const Vector3 b = coords.col(fn[1]);
  const Vector3 c = coords.col(fn[2]);
  const Scalar area = 0.5 * (b - a).cross(c - a).norm();
  if (!(area > 0.0)) throw MeshError("zero-area element face in traction load");
  Vector fe = Vector::Zero(num_dofs());
  const Scalar share = area / 3.0;
  for (int node : fn) {
    fe(3 * node + 0) += share * traction.x();
    fe(3 * node + 1) += share * traction.y();
    fe(3 * node + 2) += share * traction.z();
  }
  return fe;
}

std::vector<IntegrationPoint> Tet4Element::integration_rule(
    const IntegrationOptions& /*opts*/) const {
  // One point at the centroid; the reference tetrahedron has volume 1/6 and
  // det J = 6V, so w det J = V.
  IntegrationPoint ip;
  ip.point = reference_centroid();
  ip.weight = 1.0 / 6.0;
  return {ip};
}

}  // namespace sparlab
