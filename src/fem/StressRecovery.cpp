#include "sparlab/fem/StressRecovery.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <cmath>
#include <sstream>

namespace sparlab {
namespace {

void gather_element_displacement(const Mesh& mesh, Index e, const Vector& u, Vector& ue) {
  const int npe = mesh.nodes_per_elem();
  ue.resize(npe * kDofsPerNode);
  const Index* nodes = mesh.element_nodes(e);
  for (int a = 0; a < npe; ++a) {
    ue(kDofsPerNode * a + 0) = u(nodes[a] * kDofsPerNode + 0);
    ue(kDofsPerNode * a + 1) = u(nodes[a] * kDofsPerNode + 1);
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

Scalar von_mises(const Vector3& s, StressState state, Scalar poisson) {
  const Scalar sx = s(0);
  const Scalar sy = s(1);
  const Scalar sxy = s(2);
  const Scalar sz = (state == StressState::PlaneStrain) ? poisson * (sx + sy) : 0.0;
  const Scalar j = 0.5 * ((sx - sy) * (sx - sy) + (sy - sz) * (sy - sz) +
                          (sz - sx) * (sz - sx)) +
                   3.0 * sxy * sxy;
  return std::sqrt(std::max(j, 0.0));
}

void principal_stresses(const Vector3& s, Scalar& s1, Scalar& s2) {
  const Scalar mean = 0.5 * (s(0) + s(1));
  const Scalar dev = 0.5 * (s(0) - s(1));
  const Scalar r = std::sqrt(dev * dev + s(2) * s(2));
  s1 = mean + r;
  s2 = mean - r;
}

Vector3 element_strain_at(const FemModel& model, Index element, const NaturalPoint& point,
                          const Vector& displacement) {
  check_displacement(model, displacement);
  Vector ue;
  gather_element_displacement(model.mesh(), element, displacement, ue);
  const StrainOperator op =
      model.element().strain_operator(model.mesh().element_coordinates(element), point);
  return op.b * ue;
}

Vector3 element_stress_at(const FemModel& model, Index element, const NaturalPoint& point,
                          const Vector& displacement, Scalar stiffness_scale) {
  const Vector3 strain = element_strain_at(model, element, point, displacement);
  return stiffness_scale * (model.constitutive() * strain);
}

StressField recover_stresses(const FemModel& model, const Assembler& assembler,
                             const Vector& displacement,
                             const Vector* stiffness_scale) {
  check_displacement(model, displacement);
  const Mesh& mesh = model.mesh();
  const Index ne = mesh.num_elements();
  const Index nn = mesh.num_nodes();
  if (stiffness_scale && stiffness_scale->size() != ne) {
    std::ostringstream os;
    os << "stiffness scale vector has length " << stiffness_scale->size()
       << " but the mesh has " << ne << " elements";
    throw ModelError(os.str());
  }

  StressField field;
  field.element_strain.setZero(kVoigt, ne);
  field.element_stress.setZero(kVoigt, ne);
  field.element_solid_stress.setZero(kVoigt, ne);
  field.element_von_mises.setZero(ne);
  field.element_solid_von_mises.setZero(ne);
  field.element_principal_max.setZero(ne);
  field.element_principal_min.setZero(ne);
  field.element_strain_energy.setZero(ne);
  field.nodal_stress.setZero(kVoigt, nn);
  field.nodal_von_mises.setZero(nn);

  Vector nodal_weight = Vector::Zero(nn);
  const std::vector<NaturalPoint> points =
      model.element().stress_evaluation_points(model.integration());
  if (points.empty()) throw ModelError("element reports no stress evaluation points");

  Vector ue;
  for (Index e = 0; e < ne; ++e) {
    gather_element_displacement(mesh, e, displacement, ue);
    const auto coords = mesh.element_coordinates(e);
    const Scalar s = stiffness_scale ? (*stiffness_scale)(e) : 1.0;

    Vector3 strain_avg = Vector3::Zero();
    for (const NaturalPoint& p : points) {
      const StrainOperator op = model.element().strain_operator(coords, p);
      strain_avg += op.b * ue;
    }
    strain_avg /= static_cast<Scalar>(points.size());

    const Vector3 solid_stress = model.constitutive() * strain_avg;
    const Vector3 macro_stress = s * solid_stress;

    field.element_strain.col(e) = strain_avg;
    field.element_stress.col(e) = macro_stress;
    field.element_solid_stress.col(e) = solid_stress;
    field.element_von_mises(e) =
        von_mises(macro_stress, model.stress_state(), model.material().poisson_ratio());
    field.element_solid_von_mises(e) =
        von_mises(solid_stress, model.stress_state(), model.material().poisson_ratio());
    Scalar s1 = 0.0;
    Scalar s2 = 0.0;
    principal_stresses(macro_stress, s1, s2);
    field.element_principal_max(e) = s1;
    field.element_principal_min(e) = s2;

    // Exact element strain energy from the element stiffness matrix.
    const Matrix& ke = assembler.element_stiffness(e);
    field.element_strain_energy(e) = 0.5 * s * ue.dot(ke * ue);

    // Area-weighted scatter to nodes for smooth plotting.
    const Scalar w = mesh.element_area(e);
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
