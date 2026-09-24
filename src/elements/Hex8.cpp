#include "sparlab/elements/Hex8.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/elements/Quad4.hpp"
#include "sparlab/elements/Quadrature.hpp"

#include <Eigen/Dense>

#include <sstream>

namespace sparlab {
namespace {

/// Natural coordinates of the eight corners in VTK order.
constexpr Scalar kCorner[8][3] = {{-1.0, -1.0, -1.0}, {1.0, -1.0, -1.0}, {1.0, 1.0, -1.0},
                                  {-1.0, 1.0, -1.0},  {-1.0, -1.0, 1.0}, {1.0, -1.0, 1.0},
                                  {1.0, 1.0, 1.0},    {-1.0, 1.0, 1.0}};

using Hex8Coords = Eigen::Matrix<Scalar, 3, 8>;

Hex8Coords hex8_coords(const Matrix& coords) {
  if (coords.rows() != 3 || coords.cols() != 8) {
    std::ostringstream os;
    os << "Hex8 expects a 3 x 8 nodal coordinate matrix, received " << coords.rows()
       << " x " << coords.cols();
    throw MeshError(os.str());
  }
  return coords;
}

void require_unit_thickness(Scalar thickness) {
  if (thickness != 1.0) {
    std::ostringstream os;
    os << "a solid (Hex8) element has no thickness; received " << thickness
       << " m. Leave model.thickness at its default of 1 for a 3-D mesh";
    throw ConfigError(os.str());
  }
}

/// Jacobian, its determinant and the Cartesian shape gradients at one point.
struct Hex8Mapping {
  Eigen::Matrix<Scalar, 8, 3> dn_dx;
  Scalar det = 0.0;
};

Hex8Mapping hex8_mapping(const Hex8Coords& coords, Scalar xi, Scalar eta, Scalar zeta,
                         const char* context) {
  const Eigen::Matrix<Scalar, 8, 3> dn_dxi = hex8_shape_gradients_natural(xi, eta, zeta);
  // J(i,j) = dx_i / dxi_j
  const Matrix3 jac = coords * dn_dxi;
  // Explicit 3x3 determinant and adjugate: this runs at every Gauss point of
  // every element and needs no LU machinery.
  const Scalar c00 = jac(1, 1) * jac(2, 2) - jac(1, 2) * jac(2, 1);
  const Scalar c01 = jac(1, 2) * jac(2, 0) - jac(1, 0) * jac(2, 2);
  const Scalar c02 = jac(1, 0) * jac(2, 1) - jac(1, 1) * jac(2, 0);
  const Scalar det = jac(0, 0) * c00 + jac(0, 1) * c01 + jac(0, 2) * c02;
  if (!(det > 0.0)) {
    std::ostringstream os;
    os << "Hex8 Jacobian determinant is " << det << " m^3 at (xi, eta, zeta) = (" << xi
       << ", " << eta << ", " << zeta << ") " << context
       << "; the element is inverted, folded or its nodes do not follow the VTK "
          "hexahedron ordering";
    throw MeshError(os.str());
  }
  Matrix3 jinv;
  jinv(0, 0) = c00 / det;
  jinv(1, 0) = c01 / det;
  jinv(2, 0) = c02 / det;
  jinv(0, 1) = (jac(0, 2) * jac(2, 1) - jac(0, 1) * jac(2, 2)) / det;
  jinv(1, 1) = (jac(0, 0) * jac(2, 2) - jac(0, 2) * jac(2, 0)) / det;
  jinv(2, 1) = (jac(0, 1) * jac(2, 0) - jac(0, 0) * jac(2, 1)) / det;
  jinv(0, 2) = (jac(0, 1) * jac(1, 2) - jac(0, 2) * jac(1, 1)) / det;
  jinv(1, 2) = (jac(0, 2) * jac(1, 0) - jac(0, 0) * jac(1, 2)) / det;
  jinv(2, 2) = (jac(0, 0) * jac(1, 1) - jac(0, 1) * jac(1, 0)) / det;

  Hex8Mapping map;
  map.det = det;
  // dN/dx = dN/dxi * J^{-1}  (row a holds grad N_a).
  map.dn_dx = dn_dxi * jinv;
  return map;
}

Eigen::Matrix<Scalar, 6, 24> hex8_strain_matrix(const Eigen::Matrix<Scalar, 8, 3>& dn_dx) {
  Eigen::Matrix<Scalar, 6, 24> b = Eigen::Matrix<Scalar, 6, 24>::Zero();
  for (int a = 0; a < 8; ++a) {
    const Scalar dx = dn_dx(a, 0);
    const Scalar dy = dn_dx(a, 1);
    const Scalar dz = dn_dx(a, 2);
    const int c = 3 * a;
    b(0, c + 0) = dx;  // eps_xx
    b(1, c + 1) = dy;  // eps_yy
    b(2, c + 2) = dz;  // eps_zz
    b(3, c + 0) = dy;  // gamma_xy = du/dy + dv/dx
    b(3, c + 1) = dx;
    b(4, c + 1) = dz;  // gamma_yz = dv/dz + dw/dy
    b(4, c + 2) = dy;
    b(5, c + 0) = dz;  // gamma_zx = du/dz + dw/dx
    b(5, c + 2) = dx;
  }
  return b;
}

}  // namespace

Eigen::Matrix<Scalar, 8, 1> hex8_shape_functions(Scalar xi, Scalar eta, Scalar zeta) {
  Eigen::Matrix<Scalar, 8, 1> n;
  for (int a = 0; a < 8; ++a) {
    n(a) = 0.125 * (1.0 + xi * kCorner[a][0]) * (1.0 + eta * kCorner[a][1]) *
           (1.0 + zeta * kCorner[a][2]);
  }
  return n;
}

Eigen::Matrix<Scalar, 8, 3> hex8_shape_gradients_natural(Scalar xi, Scalar eta,
                                                         Scalar zeta) {
  Eigen::Matrix<Scalar, 8, 3> g;
  for (int a = 0; a < 8; ++a) {
    const Scalar fx = 1.0 + xi * kCorner[a][0];
    const Scalar fy = 1.0 + eta * kCorner[a][1];
    const Scalar fz = 1.0 + zeta * kCorner[a][2];
    g(a, 0) = 0.125 * kCorner[a][0] * fy * fz;
    g(a, 1) = 0.125 * kCorner[a][1] * fx * fz;
    g(a, 2) = 0.125 * kCorner[a][2] * fx * fy;
  }
  return g;
}

Scalar hex8_volume(const Matrix& coords_in, Scalar* min_det) {
  const Hex8Coords coords = hex8_coords(coords_in);
  Scalar volume = 0.0;
  Scalar lowest = std::numeric_limits<Scalar>::max();
  for (const auto& gp : gauss_legendre_cube(2)) {
    const Eigen::Matrix<Scalar, 8, 3> dn_dxi =
        hex8_shape_gradients_natural(gp.xi, gp.eta, gp.zeta);
    const Matrix3 jac = coords * dn_dxi;
    const Scalar det = jac.determinant();
    lowest = std::min(lowest, det);
    volume += gp.weight * det;
  }
  if (min_det != nullptr) *min_det = lowest;
  return volume;
}

Scalar hex8_face_area(const Matrix& face) {
  if (face.rows() != 3 || face.cols() != 4) {
    throw MeshError("hex8_face_area expects a 3 x 4 coordinate matrix");
  }
  Scalar area = 0.0;
  for (const auto& gp : gauss_legendre_square(2)) {
    const Eigen::Matrix<Scalar, 4, 2> dn = quad4_shape_gradients_natural(gp.xi, gp.eta);
    const Vector3 xs = face * dn.col(0);
    const Vector3 xt = face * dn.col(1);
    area += gp.weight * xs.cross(xt).norm();
  }
  return area;
}

Vector Hex8Element::shape_functions(const NaturalPoint& point) const {
  return hex8_shape_functions(point.xi, point.eta, point.zeta);
}

StrainOperator Hex8Element::strain_operator(const Matrix& coords,
                                            const NaturalPoint& point) const {
  const Hex8Mapping map =
      hex8_mapping(hex8_coords(coords), point.xi, point.eta, point.zeta, "");
  StrainOperator op;
  op.detJ = map.det;
  op.b = hex8_strain_matrix(map.dn_dx);
  return op;
}

Matrix Hex8Element::stiffness(const Matrix& coords_in, const Matrix& d_in,
                              Scalar thickness, const IntegrationOptions& opts) const {
  require_unit_thickness(thickness);
  if (d_in.rows() != 6 || d_in.cols() != 6) {
    std::ostringstream os;
    os << "Hex8 expects a 6 x 6 constitutive matrix, received " << d_in.rows() << " x "
       << d_in.cols();
    throw ModelError(os.str());
  }
  const Hex8Coords coords = hex8_coords(coords_in);
  const Matrix6 d = d_in;
  Eigen::Matrix<Scalar, 24, 24> ke = Eigen::Matrix<Scalar, 24, 24>::Zero();
  for (const auto& gp : gauss_legendre_cube(opts.stiffness_points)) {
    const Hex8Mapping map = hex8_mapping(coords, gp.xi, gp.eta, gp.zeta,
                                         "during stiffness integration");
    const Eigen::Matrix<Scalar, 6, 24> b = hex8_strain_matrix(map.dn_dx);
    ke.noalias() += (map.det * gp.weight) * (b.transpose() * d * b);
  }
  // Symmetrise: the integrand is symmetric, so any asymmetry is round-off.
  return 0.5 * (ke + ke.transpose());
}

Matrix Hex8Element::consistent_mass(const Matrix& coords_in, Scalar density,
                                    Scalar thickness,
                                    const IntegrationOptions& opts) const {
  require_unit_thickness(thickness);
  const Hex8Coords coords = hex8_coords(coords_in);
  Eigen::Matrix<Scalar, 24, 24> me = Eigen::Matrix<Scalar, 24, 24>::Zero();
  for (const auto& gp : gauss_legendre_cube(opts.mass_points)) {
    const Eigen::Matrix<Scalar, 8, 1> n = hex8_shape_functions(gp.xi, gp.eta, gp.zeta);
    const Hex8Mapping map =
        hex8_mapping(coords, gp.xi, gp.eta, gp.zeta, "during mass integration");
    // N_mat is 3 x 24 with the node-major DOF ordering.
    Eigen::Matrix<Scalar, 3, 24> nmat = Eigen::Matrix<Scalar, 3, 24>::Zero();
    for (int a = 0; a < 8; ++a) {
      nmat(0, 3 * a + 0) = n(a);
      nmat(1, 3 * a + 1) = n(a);
      nmat(2, 3 * a + 2) = n(a);
    }
    me.noalias() += (density * map.det * gp.weight) * (nmat.transpose() * nmat);
  }
  return 0.5 * (me + me.transpose());
}

std::vector<NaturalPoint> Hex8Element::stress_evaluation_points(
    const IntegrationOptions& opts) const {
  std::vector<NaturalPoint> points;
  for (const auto& gp : gauss_legendre_cube(opts.stiffness_points)) {
    NaturalPoint p;
    p.xi = gp.xi;
    p.eta = gp.eta;
    p.zeta = gp.zeta;
    points.push_back(p);
  }
  return points;
}

Vector Hex8Element::boundary_traction(const Matrix& coords_in, int local_face,
                                      const Vector3& traction, Scalar thickness,
                                      const IntegrationOptions& opts) const {
  require_unit_thickness(thickness);
  const std::vector<int>& fn = face_nodes(local_face);
  const Hex8Coords coords = hex8_coords(coords_in);

  // The face is a bilinear quadrilateral parametrised by (s, t) in [-1, 1]^2
  // through the Q4 shape functions in the face's own node order; its area
  // element is |x_s x x_t| ds dt.
  Eigen::Matrix<Scalar, 3, 4> xf;
  for (int a = 0; a < 4; ++a) xf.col(a) = coords.col(fn[static_cast<std::size_t>(a)]);

  Vector fe = Vector::Zero(num_dofs());
  Scalar area = 0.0;
  for (const auto& gp : gauss_legendre_square(opts.edge_points)) {
    const Eigen::Vector4d n = quad4_shape_functions(gp.xi, gp.eta);
    const Eigen::Matrix<Scalar, 4, 2> dn = quad4_shape_gradients_natural(gp.xi, gp.eta);
    const Vector3 xs = xf * dn.col(0);
    const Vector3 xt = xf * dn.col(1);
    const Scalar da = xs.cross(xt).norm();
    area += gp.weight * da;
    for (int a = 0; a < 4; ++a) {
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
