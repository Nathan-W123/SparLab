#include "sparlab/elements/Tri3.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <sstream>

namespace sparlab {
namespace {

void check_coords(const Matrix& coords) {
  if (coords.rows() != 2 || coords.cols() != 3) {
    std::ostringstream os;
    os << "Tri3 expects a 2 x 3 nodal coordinate matrix, received " << coords.rows()
       << " x " << coords.cols();
    throw MeshError(os.str());
  }
}

void check_thickness(Scalar thickness) {
  if (!(thickness > 0.0)) {
    std::ostringstream os;
    os << "element thickness must be positive (got " << thickness << " m)";
    throw ConfigError(os.str());
  }
}

Matrix strain_matrix(const Eigen::Matrix<Scalar, 3, 2>& grad) {
  Matrix b = Matrix::Zero(3, 6);
  for (int a = 0; a < 3; ++a) {
    const Scalar dx = grad(a, 0);
    const Scalar dy = grad(a, 1);
    b(0, 2 * a + 0) = dx;  // eps_xx = du/dx
    b(1, 2 * a + 1) = dy;  // eps_yy = dv/dy
    b(2, 2 * a + 0) = dy;  // gamma_xy = du/dy + dv/dx
    b(2, 2 * a + 1) = dx;
  }
  return b;
}

}  // namespace

Eigen::Matrix<Scalar, 3, 2> tri3_shape_gradients(const Matrix& coords, Scalar* area) {
  check_coords(coords);
  const Scalar x1 = coords(0, 0), y1 = coords(1, 0);
  const Scalar x2 = coords(0, 1), y2 = coords(1, 1);
  const Scalar x3 = coords(0, 2), y3 = coords(1, 2);
  const Scalar twice_area = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
  if (!(twice_area > 0.0)) {
    std::ostringstream os;
    os << "Tri3 signed area is " << 0.5 * twice_area
       << " m^2; the element is inverted, collapsed or its nodes are not "
          "counter-clockwise";
    throw MeshError(os.str());
  }
  Eigen::Matrix<Scalar, 3, 2> g;
  g(0, 0) = (y2 - y3) / twice_area;  g(0, 1) = (x3 - x2) / twice_area;
  g(1, 0) = (y3 - y1) / twice_area;  g(1, 1) = (x1 - x3) / twice_area;
  g(2, 0) = (y1 - y2) / twice_area;  g(2, 1) = (x2 - x1) / twice_area;
  if (area != nullptr) *area = 0.5 * twice_area;
  return g;
}

NaturalPoint Tri3Element::reference_centroid() const {
  NaturalPoint p;
  p.xi = 1.0 / 3.0;
  p.eta = 1.0 / 3.0;
  return p;
}

Vector Tri3Element::shape_functions(const NaturalPoint& point) const {
  Vector n(3);
  n << 1.0 - point.xi - point.eta, point.xi, point.eta;
  return n;
}

StrainOperator Tri3Element::strain_operator(const Matrix& coords,
                                            const NaturalPoint& /*point*/) const {
  Scalar area = 0.0;
  const Eigen::Matrix<Scalar, 3, 2> g = tri3_shape_gradients(coords, &area);
  StrainOperator op;
  op.b = strain_matrix(g);
  op.detJ = 2.0 * area;
  return op;
}

Matrix Tri3Element::stiffness(const Matrix& coords, const Matrix& d, Scalar thickness,
                              const IntegrationOptions& /*opts*/) const {
  check_thickness(thickness);
  if (d.rows() != 3 || d.cols() != 3) {
    std::ostringstream os;
    os << "Tri3 expects a 3 x 3 constitutive matrix, received " << d.rows() << " x "
       << d.cols();
    throw ModelError(os.str());
  }
  Scalar area = 0.0;
  const Matrix b = strain_matrix(tri3_shape_gradients(coords, &area));
  const Matrix ke = (thickness * area) * (b.transpose() * d * b);
  return 0.5 * (ke + ke.transpose());
}

Matrix Tri3Element::consistent_mass(const Matrix& coords, Scalar density, Scalar thickness,
                                    const IntegrationOptions& /*opts*/) const {
  check_thickness(thickness);
  Scalar area = 0.0;
  tri3_shape_gradients(coords, &area);
  const Scalar m = density * thickness * area / 12.0;
  Matrix me = Matrix::Zero(6, 6);
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      const Scalar value = (a == b ? 2.0 : 1.0) * m;
      me(2 * a + 0, 2 * b + 0) = value;
      me(2 * a + 1, 2 * b + 1) = value;
    }
  }
  return me;
}

std::vector<NaturalPoint> Tri3Element::stress_evaluation_points(
    const IntegrationOptions& /*opts*/) const {
  return {reference_centroid()};
}

Vector Tri3Element::boundary_traction(const Matrix& coords, int local_face,
                                      const Vector3& traction, Scalar thickness,
                                      const IntegrationOptions& /*opts*/) const {
  const std::vector<int>& en = face_nodes(local_face);
  if (traction.z() != 0.0) {
    std::ostringstream os;
    os << "a plane element cannot carry an out-of-plane traction (t_z = " << traction.z()
       << " Pa); use a 3-D mesh or drop the z component";
    throw ConfigError(os.str());
  }
  check_coords(coords);
  check_thickness(thickness);
  const Scalar length = (coords.col(en[1]) - coords.col(en[0])).norm();
  if (!(length > 0.0)) throw MeshError("zero-length element edge in traction load");
  Vector fe = Vector::Zero(num_dofs());
  const Scalar share = 0.5 * thickness * length;
  for (int node : en) {
    fe(2 * node + 0) += share * traction.x();
    fe(2 * node + 1) += share * traction.y();
  }
  return fe;
}

std::vector<IntegrationPoint> Tri3Element::integration_rule(
    const IntegrationOptions& /*opts*/) const {
  // One point at the centroid; the reference triangle has area 1/2 and
  // det J = 2A, so w det J = A.
  IntegrationPoint ip;
  ip.point = reference_centroid();
  ip.weight = 0.5;
  return {ip};
}

}  // namespace sparlab
