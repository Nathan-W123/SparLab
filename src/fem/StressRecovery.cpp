#include "sparlab/fem/StressRecovery.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <Eigen/Eigenvalues>

#include <cmath>
#include <sstream>

namespace sparlab {
namespace {

void gather_element_displacement(const Mesh& mesh, Index e, const Vector& u, Vector& ue) {
  const int npe = mesh.nodes_per_elem();
  const int dim = mesh.dim();
  ue.resize(npe * dim);
  const Index* nodes = mesh.element_nodes(e);
  for (int a = 0; a < npe; ++a) {
    for (int k = 0; k < dim; ++k) ue(dim * a + k) = u(nodes[a] * dim + k);
  }
}

void check_displacement(const FemModel& model, const Vector& u) {
  if (u.size() != model.dofs().num_dofs()) {
    std::ostringstream os;
    os << "displacement vector has length " << u.size() << " but the model has "
       << model.dofs().num_dofs() << " DOFs";
    throw ModelError(os.str());
  }
}

}  // namespace

Scalar von_mises(const Vector& s, StressState state, Scalar poisson) {
  if (s.size() == 3) {
    const Scalar sx = s(0);
    const Scalar sy = s(1);
    const Scalar sxy = s(2);
    const Scalar sz = (state == StressState::PlaneStrain) ? poisson * (sx + sy) : 0.0;
    const Scalar j = 0.5 * ((sx - sy) * (sx - sy) + (sy - sz) * (sy - sz) +
                            (sz - sx) * (sz - sx)) +
                     3.0 * sxy * sxy;
    return std::sqrt(std::max(j, 0.0));
  }
  if (s.size() == 6) {
    const Scalar sx = s(0);
    const Scalar sy = s(1);
    const Scalar sz = s(2);
    const Scalar j = 0.5 * ((sx - sy) * (sx - sy) + (sy - sz) * (sy - sz) +
                            (sz - sx) * (sz - sx)) +
                     3.0 * (s(3) * s(3) + s(4) * s(4) + s(5) * s(5));
    return std::sqrt(std::max(j, 0.0));
  }
  std::ostringstream os;
  os << "von_mises expects a Voigt vector of 3 or 6 components, received " << s.size();
  throw ModelError(os.str());
}

void principal_stresses(const Vector& s, Scalar& s1, Scalar& s2) {
  if (s.size() != 3) {
    throw ModelError("principal_stresses expects a 3-component (plane) Voigt vector");
  }
  const Scalar mean = 0.5 * (s(0) + s(1));
  const Scalar dev = 0.5 * (s(0) - s(1));
  const Scalar r = std::sqrt(dev * dev + s(2) * s(2));
  s1 = mean + r;
  s2 = mean - r;
}

Vector3 principal_stresses_3d(const Vector& s) {
  if (s.size() != 6) {
    throw ModelError("principal_stresses_3d expects a 6-component Voigt vector");
  }
  Matrix3 sigma;
  sigma << s(0), s(3), s(5),
           s(3), s(1), s(4),
           s(5), s(4), s(2);
  Eigen::SelfAdjointEigenSolver<Matrix3> es(sigma, Eigen::EigenvaluesOnly);
  if (es.info() != Eigen::Success) {
    throw SolverError("the 3 x 3 principal-stress eigensolve failed");
  }
  // Eigen returns ascending order; report descending.
  return es.eigenvalues().reverse();
}

Vector element_strain_at(const FemModel& model, Index element, const NaturalPoint& point,
                         const Vector& displacement) {
  check_displacement(model, displacement);
  Vector ue;
  gather_element_displacement(model.mesh(), element, displacement, ue);
  const StrainOperator op =
      model.element().strain_operator(model.mesh().element_coordinates(element), point);
  return op.b * ue;
}

Vector element_stress_at(const FemModel& model, Index element, const NaturalPoint& point,
                         const Vector& displacement, Scalar stiffness_scale) {
  const Vector strain = element_strain_at(model, element, point, displacement);
  return stiffness_scale * (model.constitutive() * strain);
}

StressField recover_stresses(const FemModel& model, const Assembler& assembler,
                             const Vector& displacement,
                             const Vector* stiffness_scale) {
  check_displacement(model, displacement);
  const Mesh& mesh = model.mesh();
  const Index ne = mesh.num_elements();
  const Index nn = mesh.num_nodes();
  const int nv = model.element().num_voigt();
  const bool solid = mesh.dim() == 3;
  if (stiffness_scale && stiffness_scale->size() != ne) {
    std::ostringstream os;
    os << "stiffness scale vector has length " << stiffness_scale->size()
       << " but the mesh has " << ne << " elements";
    throw ModelError(os.str());
  }

  StressField field;
  field.element_strain.setZero(nv, ne);
  field.element_stress.setZero(nv, ne);
  field.element_solid_stress.setZero(nv, ne);
  field.element_von_mises.setZero(ne);
  field.element_solid_von_mises.setZero(ne);
  field.element_principal_max.setZero(ne);
  if (solid) field.element_principal_mid.setZero(ne);
  field.element_principal_min.setZero(ne);
  field.element_strain_energy.setZero(ne);
  field.nodal_stress.setZero(nv, nn);
  field.nodal_von_mises.setZero(nn);

  Vector nodal_weight = Vector::Zero(nn);
  const std::vector<NaturalPoint> points =
      model.element().stress_evaluation_points(model.integration());
  if (points.empty()) throw ModelError("element reports no stress evaluation points");

  Vector ue;
  for (Index e = 0; e < ne; ++e) {
    gather_element_displacement(mesh, e, displacement, ue);
    const Matrix coords = mesh.element_coordinates(e);
    const Scalar s = stiffness_scale ? (*stiffness_scale)(e) : 1.0;

    Vector strain_avg = Vector::Zero(nv);
    for (const NaturalPoint& p : points) {
      const StrainOperator op = model.element().strain_operator(coords, p);
      strain_avg += op.b * ue;
    }
    strain_avg /= static_cast<Scalar>(points.size());

    const Vector solid_stress = model.constitutive() * strain_avg;
    const Vector macro_stress = s * solid_stress;

    field.element_strain.col(e) = strain_avg;
    field.element_stress.col(e) = macro_stress;
    field.element_solid_stress.col(e) = solid_stress;
    field.element_von_mises(e) =
        von_mises(macro_stress, model.stress_state(), model.material().poisson_ratio());
    field.element_solid_von_mises(e) =
        von_mises(solid_stress, model.stress_state(), model.material().poisson_ratio());
    if (solid) {
      const Vector3 principal = principal_stresses_3d(macro_stress);
      field.element_principal_max(e) = principal(0);
      field.element_principal_mid(e) = principal(1);
      field.element_principal_min(e) = principal(2);
    } else {
      Scalar s1 = 0.0;
      Scalar s2 = 0.0;
      principal_stresses(macro_stress, s1, s2);
      field.element_principal_max(e) = s1;
      field.element_principal_min(e) = s2;
    }

    // Exact element strain energy from the element stiffness matrix.
    const Matrix& ke = assembler.element_stiffness(e);
    field.element_strain_energy(e) = 0.5 * s * ue.dot(ke * ue);

    // Measure-weighted scatter to nodes for smooth plotting.
    const Scalar w = mesh.element_measure(e);
    const Index* nodes = mesh.element_nodes(e);
    for (int a = 0; a < mesh.nodes_per_elem(); ++a) {
      field.nodal_stress.col(nodes[a]) += w * macro_stress;
      nodal_weight(nodes[a]) += w;
    }
  }

  for (Index n = 0; n < nn; ++n) {
    if (nodal_weight(n) > 0.0) field.nodal_stress.col(n) /= nodal_weight(n);
    field.nodal_von_mises(n) = von_mises(field.nodal_stress.col(n), model.stress_state(),
                                         model.material().poisson_ratio());
  }

  return field;
}

}  // namespace sparlab
