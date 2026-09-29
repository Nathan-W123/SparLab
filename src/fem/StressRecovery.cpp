#include "sparlab/fem/StressRecovery.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/Loads.hpp"

#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace sparlab {
namespace {

void gather_element_displacement(const FemModel& model, Index e, const Vector& u, Vector& ue) {
  const Mesh& mesh = model.mesh();
  model.dofs().gather(mesh.element_nodes(e), mesh.nodes_per_elem(), u, ue);
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

Scalar von_mises_plane(Scalar sx, Scalar sy, Scalar sxy, Scalar sz) {
  const Scalar j = 0.5 * ((sx - sy) * (sx - sy) + (sy - sz) * (sy - sz) + (sz - sx) * (sz - sx)) +
                   3.0 * sxy * sxy;
  return std::sqrt(std::max(j, 0.0));
}

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
  gather_element_displacement(model, element, displacement, ue);
  const StrainOperator op =
      model.element().strain_operator(model.mesh().element_coordinates(element), point);
  return op.b * ue;
}

Vector element_stress_at(const FemModel& model, Index element, const NaturalPoint& point,
                         const Vector& displacement, Scalar stiffness_scale,
                         const Vector* temperature) {
  Vector strain = element_strain_at(model, element, point, displacement);
  if (temperature != nullptr) strain -= element_thermal_strain(model, element, point, *temperature);
  return stiffness_scale * (model.constitutive_of(element) * strain);
}

ShellField recover_shell_resultants(const FemModel& model, const Assembler& assembler,
                                    const Vector& displacement) {
  check_displacement(model, displacement);
  const auto* shell = dynamic_cast<const Shell4Element*>(&model.element());
  if (shell == nullptr) {
    throw ModelError("recover_shell_resultants needs a shell model; a continuum model's "
                     "stresses come from recover_stresses");
  }
  const Mesh& mesh = model.mesh();
  const Index ne = mesh.num_elements();
  const Index nn = mesh.num_nodes();
  ShellField field;
  field.element.resize(static_cast<std::size_t>(ne));
  field.element_von_mises.setZero(ne);
  field.element_strain_energy.setZero(ne);
  field.nodal_von_mises.setZero(nn);
  Vector nodal_weight = Vector::Zero(nn);
  Vector ue;
  for (Index e = 0; e < ne; ++e) {
    gather_element_displacement(model, e, displacement, ue);
    const ShellResultants r = shell->resultants(model.element_geometry(e),
                                                model.constitutive_of(e), ue,
                                                model.thickness_of(e), 0.0, 0.0);
    field.element[static_cast<std::size_t>(e)] = r;
    const Scalar vm = std::max({r.von_mises_top, r.von_mises_bottom, r.von_mises_mid});
    field.element_von_mises(e) = vm;
    field.element_strain_energy(e) = 0.5 * ue.dot(assembler.element_stiffness(e) * ue);
    const Scalar w = mesh.element_measure(e);
    const Index* nodes = mesh.element_nodes(e);
    for (int a = 0; a < mesh.nodes_per_elem(); ++a) {
      field.nodal_von_mises(nodes[a]) += w * vm;
      nodal_weight(nodes[a]) += w;
    }
  }
  for (Index n = 0; n < nn; ++n) {
    if (nodal_weight(n) > 0.0) field.nodal_von_mises(n) /= nodal_weight(n);
  }
  return field;
}

std::vector<Vector3> beam_distributed_loads(const FemModel& model,
                                            const LoadCaseSpec& load_case) {
  if (!model.is_beam()) throw ModelError("beam_distributed_loads needs a beam model");
  const Mesh& mesh = model.mesh();
  std::vector<Vector3> q = beam_line_loads(model, load_case);
  const Index ne = mesh.num_elements();
  std::vector<Vector3> density(static_cast<std::size_t>(ne), Vector3::Zero());
  for (Index e = 0; e < ne; ++e) {
    density[static_cast<std::size_t>(e)] = model.material_of(e).density() * load_case.gravity;
  }
  for (const BodyForceSpec& body : load_case.body_forces) {
    if (body.whole_model) {
      for (Vector3& b : density) b += body.force_density;
    } else {
      for (Index e : body.region.select_elements(mesh)) {
        density[static_cast<std::size_t>(e)] += body.force_density;
      }
    }
  }
  for (Index e = 0; e < ne; ++e) {
    const Vector3& b = density[static_cast<std::size_t>(e)];
    if (b.squaredNorm() > 0.0) q[static_cast<std::size_t>(e)] += model.section_of(e).area * b;
  }
  return q;
}

BeamField recover_beam_forces(const FemModel& model, const Assembler& assembler,
                              const Vector& displacement, const LoadCaseSpec& load_case) {
  check_displacement(model, displacement);
  const auto* beam = dynamic_cast<const Beam2Element*>(&model.element());
  if (beam == nullptr) {
    throw ModelError("recover_beam_forces needs a beam model; a continuum model's stresses "
                     "come from recover_stresses");
  }
  const Mesh& mesh = model.mesh();
  const Index ne = mesh.num_elements();
  const Index nn = mesh.num_nodes();
  const std::vector<Vector3> q = beam_distributed_loads(model, load_case);
  BeamField field;
  field.element.resize(static_cast<std::size_t>(ne));
  field.element_normal_stress.setZero(ne);
  field.element_strain_energy.setZero(ne);
  field.nodal_normal_stress.setZero(nn);
  Vector nodal_weight = Vector::Zero(nn);
  Vector ue;
  for (Index e = 0; e < ne; ++e) {
    gather_element_displacement(model, e, displacement, ue);
    const BeamSection s = model.section_of(e);
    const BeamEndForces f = beam->end_forces(model.element_geometry(e), model.constitutive_of(e),
                                             ue, q[static_cast<std::size_t>(e)]);
    field.element[static_cast<std::size_t>(e)] = f;
    Scalar sigma = std::numeric_limits<Scalar>::quiet_NaN();
    if (s.fibre_y > 0.0 || s.fibre_z > 0.0) {
      sigma = 0.0;
      for (const auto* r : {&f.start, &f.end}) {
        const Scalar n = std::abs((*r)(0)) / s.area;
        const Scalar bending =
            s.round() ? std::hypot((*r)(4), (*r)(5)) * s.fibre_y / s.iy
                      : std::abs((*r)(4)) * s.fibre_z / s.iy + std::abs((*r)(5)) * s.fibre_y / s.iz;
        sigma = std::max(sigma, n + bending);
      }
    }
    field.element_normal_stress(e) = sigma;
    field.element_strain_energy(e) = 0.5 * ue.dot(assembler.element_stiffness(e) * ue);
    if (std::isfinite(sigma)) {
      const Scalar w = mesh.element_measure(e);
      const Index* nodes = mesh.element_nodes(e);
      for (int a = 0; a < 2; ++a) {
        field.nodal_normal_stress(nodes[a]) += w * sigma;
        nodal_weight(nodes[a]) += w;
      }
    }
  }
  for (Index n = 0; n < nn; ++n) {
    if (nodal_weight(n) > 0.0) field.nodal_normal_stress(n) /= nodal_weight(n);
  }
  return field;
}

StressField recover_stresses(const FemModel& model, const Assembler& assembler,
                             const Vector& displacement,
                             const Vector* stiffness_scale,
                             const Vector* temperature) {
  check_displacement(model, displacement);
  if (model.dofs_per_node() != model.dim()) {
    throw ModelError("recover_stresses is the continuum recovery; a shell or beam model "
                     "reports its resultants through recover_shell_resultants or "
                     "recover_beam_forces");
  }
  if (temperature != nullptr && temperature->size() != model.mesh().num_nodes()) {
    throw ModelError("recover_stresses: the temperature field does not match the mesh");
  }
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
  const bool plane_strain = model.stress_state() == StressState::PlaneStrain;
  if (plane_strain) field.element_sigma_zz.setZero(ne);
  if (temperature != nullptr) field.element_temperature.setZero(ne);
  Vector nodal_sigma_zz = plane_strain ? Vector::Zero(nn) : Vector();

  Vector nodal_weight = Vector::Zero(nn);
  const std::vector<NaturalPoint> points =
      model.element().stress_evaluation_points(model.integration());
  if (points.empty()) throw ModelError("element reports no stress evaluation points");
  const std::vector<IntegrationPoint> rule = model.element().integration_rule(model.integration());
  const Scalar t = mesh.dim() == 2 ? model.thickness() : 1.0;

  Vector ue;
  for (Index e = 0; e < ne; ++e) {
    gather_element_displacement(model, e, displacement, ue);
    const Matrix coords = mesh.element_coordinates(e);
    const Scalar s = stiffness_scale ? (*stiffness_scale)(e) : 1.0;
    const IsotropicMaterial& material = model.material_of(e);
    const Matrix& d = model.constitutive_of(e);
    const bool thermal = temperature != nullptr && material.thermal_expansion() != 0.0;

    Vector strain_avg = Vector::Zero(nv);
    Vector eps0_avg = Vector::Zero(nv);
    Scalar dt_avg = 0.0;
    for (const NaturalPoint& p : points) {
      const StrainOperator op = model.element().strain_operator(coords, p);
      strain_avg += op.b * ue;
      if (temperature != nullptr) {
        const Scalar dt = element_temperature_change(model, e, p, *temperature);
        dt_avg += dt;
        if (thermal) eps0_avg += material.thermal_strain(model.stress_state(), dt);
      }
    }
    strain_avg /= static_cast<Scalar>(points.size());
    eps0_avg /= static_cast<Scalar>(points.size());
    dt_avg /= static_cast<Scalar>(points.size());

    const Vector solid_stress = thermal ? Vector(d * (strain_avg - eps0_avg))
                                        : Vector(d * strain_avg);
    const Vector macro_stress = s * solid_stress;

    field.element_strain.col(e) = strain_avg;
    field.element_stress.col(e) = macro_stress;
    field.element_solid_stress.col(e) = solid_stress;
    if (temperature != nullptr) {
      field.element_temperature(e) = dt_avg + material.reference_temperature();
    }
    if (plane_strain && thermal) {
      // sigma_zz = nu (sxx + syy) - E alpha dT, times the stiffness factor.
      const Scalar solid_zz =
          material.plane_strain_sigma_zz(solid_stress(0), solid_stress(1), dt_avg);
      field.element_sigma_zz(e) = s * solid_zz;
      field.element_von_mises(e) =
          von_mises_plane(macro_stress(0), macro_stress(1), macro_stress(2), s * solid_zz);
      field.element_solid_von_mises(e) =
          von_mises_plane(solid_stress(0), solid_stress(1), solid_stress(2), solid_zz);
    } else {
      if (plane_strain) {
        field.element_sigma_zz(e) =
            material.poisson_ratio() * (macro_stress(0) + macro_stress(1));
      }
      field.element_von_mises(e) =
          von_mises(macro_stress, model.stress_state(), material.poisson_ratio());
      field.element_solid_von_mises(e) =
          von_mises(solid_stress, model.stress_state(), material.poisson_ratio());
    }
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

    if (thermal) {
      // Elastic strain energy 1/2 int (Bu - eps0)^T D (Bu - eps0) dV with the
      // stiffness rule, the quadrature of the thermal load itself.
      Scalar energy = 0.0;
      for (const IntegrationPoint& ip : rule) {
        const StrainOperator op = model.element().strain_operator(coords, ip.point);
        const Vector elastic =
            op.b * ue - element_thermal_strain(model, e, ip.point, *temperature);
        energy += 0.5 * t * ip.weight * op.detJ * elastic.dot(d * elastic);
      }
      field.element_strain_energy(e) = s * energy;
    } else {
      // Exact element strain energy from the element stiffness matrix.
      const Matrix& ke = assembler.element_stiffness(e);
      field.element_strain_energy(e) = 0.5 * s * ue.dot(ke * ue);
    }

    // Measure-weighted scatter to nodes for smooth plotting.
    const Scalar w = mesh.element_measure(e);
    const Index* nodes = mesh.element_nodes(e);
    for (int a = 0; a < mesh.nodes_per_elem(); ++a) {
      field.nodal_stress.col(nodes[a]) += w * macro_stress;
      if (plane_strain) nodal_sigma_zz(nodes[a]) += w * field.element_sigma_zz(e);
      nodal_weight(nodes[a]) += w;
    }
  }

  for (Index n = 0; n < nn; ++n) {
    if (nodal_weight(n) > 0.0) {
      field.nodal_stress.col(n) /= nodal_weight(n);
      if (plane_strain) nodal_sigma_zz(n) /= nodal_weight(n);
    }
    if (plane_strain && temperature != nullptr) {
      const auto s = field.nodal_stress.col(n);
      field.nodal_von_mises(n) = von_mises_plane(s(0), s(1), s(2), nodal_sigma_zz(n));
    } else {
      field.nodal_von_mises(n) = von_mises(field.nodal_stress.col(n), model.stress_state(),
                                           model.material().poisson_ratio());
    }
  }

  return field;
}

}  // namespace sparlab
