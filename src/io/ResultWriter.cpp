#include "sparlab/io/ResultWriter.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/core/Version.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/io/VtkWriter.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace sparlab {
namespace {

std::string sanitise(const std::string& name) {
  std::string out;
  out.reserve(name.size());
  for (char c : name) {
    const bool ok = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
                    (c >= '0' && c <= '9') || c == '-' || c == '_';
    out.push_back(ok ? c : '_');
  }
  if (out.empty()) out = "unnamed";
  return out;
}

json::Value vector2_json(const Vector2& v) {
  json::Value out = json::Value::make_array();
  out.push_back(json::Value::make_number(v.x()));
  out.push_back(json::Value::make_number(v.y()));
  return out;
}

json::Value equilibrium_json(const EquilibriumCheck& eq) {
  json::Value out = json::Value::make_object();
  out.set("applied_force_N", vector2_json(eq.applied_force));
  out.set("reaction_force_N", vector2_json(eq.reaction_force));
  out.set("force_residual_N", vector2_json(eq.force_residual));
  out.set("applied_moment_Nm", json::Value::make_number(eq.applied_moment));
  out.set("reaction_moment_Nm", json::Value::make_number(eq.reaction_moment));
  out.set("moment_residual_Nm", json::Value::make_number(eq.moment_residual));
  out.set("relative_force_error", json::Value::make_number(eq.relative_force_error));
  out.set("relative_moment_error", json::Value::make_number(eq.relative_moment_error));
  return out;
}

json::Value modal_json(const ModalResult& modal) {
  json::Value out = json::Value::make_object();
  out.set("num_modes", json::Value::make_number(modal.eigenvalues.size()));
  out.set("eigenvalues_1_per_s2", json::array_of(modal.eigenvalues));
  out.set("angular_frequencies_rad_per_s", json::array_of(modal.angular_frequencies));
  out.set("frequencies_hz", json::array_of(modal.frequencies_hz));
  out.set("eigenpair_residuals", json::array_of(modal.modal_residuals));
  if (modal.modal_mass_fraction.size() > 0) {
    out.set("low_density_kinetic_energy_fraction",
            json::array_of(modal.modal_mass_fraction));
  }
  out.set("total_mass_kg", json::Value::make_number(modal.total_mass));
  out.set("subspace_size", json::Value::make_number(modal.subspace_size));
  out.set("iterations", json::Value::make_number(modal.iterations));
  out.set("converged", json::Value::make_bool(modal.converged));
  out.set("final_relative_change", json::Value::make_number(modal.final_change));
  out.set("warnings", json::array_of(modal.warnings));
  return out;
}

json::Value timings_json(const TimingLedger& timings) {
  json::Value out = json::Value::make_object();
  for (const auto& kv : timings.totals()) {
    out.set(kv.first, json::Value::make_number(kv.second));
  }
  return out;
}

json::Value mesh_stats_json(const FemModel& model) {
  const Mesh& mesh = model.mesh();
  json::Value out = json::Value::make_object();
  out.set("element_type", json::Value::make_string(to_string(mesh.element_type())));
  out.set("num_nodes", json::Value::make_number(mesh.num_nodes()));
  out.set("num_elements", json::Value::make_number(mesh.num_elements()));
  out.set("num_dofs", json::Value::make_number(model.dofs().num_dofs()));
  out.set("num_free_dofs", json::Value::make_number(model.dofs().num_free()));
  out.set("num_prescribed_dofs",
          json::Value::make_number(model.dofs().num_constrained()));
  const Eigen::Vector4d bb = mesh.bounding_box();
  json::Value box = json::Value::make_array();
  for (int i = 0; i < 4; ++i) box.push_back(json::Value::make_number(bb(i)));
  out.set("bounding_box_xmin_ymin_xmax_ymax_m", box);
  out.set("thickness_m", json::Value::make_number(model.thickness()));
  out.set("domain_volume_m3", json::Value::make_number(model.domain_volume()));
  if (mesh.structured_info().has_value()) {
    const StructuredGridInfo& info = *mesh.structured_info();
    json::Value grid = json::Value::make_object();
    grid.set("nx", json::Value::make_number(info.nx));
    grid.set("ny", json::Value::make_number(info.ny));
    grid.set("lx_m", json::Value::make_number(info.lx));
    grid.set("ly_m", json::Value::make_number(info.ly));
    grid.set("uniform", json::Value::make_bool(info.uniform));
    out.set("structured_grid", grid);
  }
  return out;
}

json::Value material_json(const IsotropicMaterial& m, StressState state) {
  json::Value out = json::Value::make_object();
  out.set("name", json::Value::make_string(m.name()));
  out.set("youngs_modulus_Pa", json::Value::make_number(m.youngs_modulus()));
  out.set("poisson_ratio", json::Value::make_number(m.poisson_ratio()));
  out.set("density_kg_per_m3", json::Value::make_number(m.density()));
  out.set("shear_modulus_Pa", json::Value::make_number(m.shear_modulus()));
  out.set("stress_state", json::Value::make_string(
                              state == StressState::PlaneStress ? "plane_stress"
                                                                : "plane_strain"));
  return out;
}

json::Value diagnostics_json(const ModelDiagnostics& diag) {
  json::Value out = json::Value::make_object();
  out.set("well_posed", json::Value::make_bool(diag.well_posed()));
  out.set("num_element_groups", json::Value::make_number(diag.components.size()));
  out.set("problems", json::array_of(diag.problems));
  json::Value groups = json::Value::make_array();
  for (const MeshComponent& c : diag.components) {
    json::Value g = json::Value::make_object();
    g.set("num_elements", json::Value::make_number(c.elements.size()));
    g.set("num_nodes", json::Value::make_number(c.nodes.size()));
    g.set("prescribed_dofs", json::Value::make_number(c.prescribed_dofs));
    g.set("rigid_body_null_dimension",
          json::Value::make_number(c.rigid_null_dimension));
    groups.push_back(g);
  }
  out.set("element_groups", groups);
  return out;
}

}  // namespace

ResultWriter::ResultWriter(std::string directory, const Configuration& config)
    : directory_(std::move(directory)), config_(config) {
  ensure_directory(directory_);
  log::debug("result directory: ", directory_);
}

std::string ResultWriter::file(const std::string& name) const {
  return path_join(directory_, name);
}

void ResultWriter::write_json(const std::string& file_name,
                              const json::Value& value) const {
  const std::string path = file(file_name);
  std::ofstream out(path, std::ios::out | std::ios::trunc);
  if (!out) throw IoError("cannot open '" + path + "' for writing");
  out << json::dump(value, 2) << '\n';
  out.flush();
  if (!out) throw IoError("failed while writing '" + path + "'");
}

void ResultWriter::write_config() const {
  write_json("config.json", config_.document);
}

void ResultWriter::write_mesh(const FemModel& model) const {
  const Mesh& mesh = model.mesh();
  json::Value doc = json::Value::make_object();
  doc.set("case", json::Value::make_string(config_.name));
  doc.set("element_type", json::Value::make_string(to_string(mesh.element_type())));
  doc.set("nodes_per_element", json::Value::make_number(mesh.nodes_per_elem()));

  json::Value nodes = json::Value::make_array();
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    nodes.push_back(vector2_json(mesh.node(n)));
  }
  doc.set("nodes_m", nodes);

  json::Value elements = json::Value::make_array();
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    json::Value row = json::Value::make_array();
    const Index* en = mesh.element_nodes(e);
    for (int a = 0; a < mesh.nodes_per_elem(); ++a) {
      row.push_back(json::Value::make_number(en[a]));
    }
    elements.push_back(row);
  }
  doc.set("elements", elements);

  // Prescribed DOFs, for the boundary-condition figure.
  json::Value constraints = json::Value::make_array();
  for (Index d : model.dofs().constrained_dofs()) {
    json::Value entry = json::Value::make_object();
    entry.set("node", json::Value::make_number(d / kDofsPerNode));
    entry.set("component",
              json::Value::make_string(d % kDofsPerNode == 0 ? "x" : "y"));
    entry.set("value_m", json::Value::make_number(model.dofs().prescribed_value(d)));
    constraints.push_back(entry);
  }
  doc.set("prescribed_dofs", constraints);

  // Applied nodal forces per load case (non-zero entries only).
  json::Value cases = json::Value::make_array();
  const std::vector<Vector>& loads = model.load_vectors();
  for (std::size_t l = 0; l < loads.size(); ++l) {
    json::Value entry = json::Value::make_object();
    entry.set("name", json::Value::make_string(model.load_case_specs()[l].name));
    entry.set("weight",
              json::Value::make_number(model.load_case_specs()[l].weight));
    json::Value forces = json::Value::make_array();
    for (Index n = 0; n < mesh.num_nodes(); ++n) {
      const Scalar fx = loads[l](n * kDofsPerNode + 0);
      const Scalar fy = loads[l](n * kDofsPerNode + 1);
      if (fx == 0.0 && fy == 0.0) continue;
      json::Value f = json::Value::make_object();
      f.set("node", json::Value::make_number(n));
      f.set("fx_N", json::Value::make_number(fx));
      f.set("fy_N", json::Value::make_number(fy));
      forces.push_back(f);
    }
    entry.set("nodal_forces", forces);
    cases.push_back(entry);
  }
  doc.set("load_cases", cases);

  write_json("mesh.json", doc);
}

void ResultWriter::write_displacement(const Mesh& mesh, const std::string& load_case,
                                      const Vector& displacement) const {
  CsvWriter csv(file("displacement_" + sanitise(load_case) + ".csv"),
                {"node", "x[m]", "y[m]", "ux[m]", "uy[m]", "umag[m]"});
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    const Vector2 x = mesh.node(n);
    const Scalar ux = displacement(n * kDofsPerNode + 0);
    const Scalar uy = displacement(n * kDofsPerNode + 1);
    csv.row(n, {x.x(), x.y(), ux, uy, std::hypot(ux, uy)});
  }
  csv.close();
}

void ResultWriter::write_stress(const Mesh& mesh, const std::string& load_case,
                                const StressField& field, const Vector* density) const {
  CsvWriter csv(file("stress_" + sanitise(load_case) + ".csv"),
                {"element", "cx[m]", "cy[m]", "area[m2]", "density[-]", "exx[-]",
                 "eyy[-]", "gxy[-]", "sxx[Pa]", "syy[Pa]", "sxy[Pa]", "von_mises[Pa]",
                 "solid_von_mises[Pa]", "principal_max[Pa]", "principal_min[Pa]",
                 "strain_energy[J]"});
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    const Vector2 c = mesh.element_centroid(e);
    csv.row(e, {c.x(), c.y(), mesh.element_area(e), density ? (*density)(e) : 1.0,
                field.element_strain(0, e), field.element_strain(1, e),
                field.element_strain(2, e), field.element_stress(0, e),
                field.element_stress(1, e), field.element_stress(2, e),
                field.element_von_mises(e), field.element_solid_von_mises(e),
                field.element_principal_max(e), field.element_principal_min(e),
                field.element_strain_energy(e)});
  }
  csv.close();
}

void ResultWriter::write_reactions(const Mesh& mesh, const DofManager& dofs,
                                   const std::string& load_case,
                                   const Vector& reactions) const {
  CsvWriter csv(file("reactions_" + sanitise(load_case) + ".csv"),
                {"node", "x[m]", "y[m]", "rx[N]", "ry[N]", "rmag[N]"});
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    const bool cx = dofs.is_constrained(n * kDofsPerNode + 0);
    const bool cy = dofs.is_constrained(n * kDofsPerNode + 1);
    if (!cx && !cy) continue;
    const Vector2 x = mesh.node(n);
    const Scalar rx = reactions(n * kDofsPerNode + 0);
    const Scalar ry = reactions(n * kDofsPerNode + 1);
    csv.row(n, {x.x(), x.y(), rx, ry, std::hypot(rx, ry)});
  }
  csv.close();
}

void ResultWriter::write_static_vtk(const Mesh& mesh, const std::string& load_case,
                                    const Vector& displacement,
                                    const StressField& field, const Vector* density,
                                    const Vector* stiffness_factor) const {
  VtkWriter writer(mesh, "SparLab static solution: " + config_.name + " / " + load_case);
  Vector magnitude(mesh.num_nodes());
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    magnitude(n) = std::hypot(displacement(n * kDofsPerNode + 0),
                              displacement(n * kDofsPerNode + 1));
  }
  writer.add_point_vectors("displacement", displacement);
  writer.add_point_scalars("displacement_magnitude", magnitude);
  writer.add_point_scalars("nodal_von_mises", field.nodal_von_mises);
  writer.add_cell_scalars("sigma_xx", field.element_stress.row(0).transpose());
  writer.add_cell_scalars("sigma_yy", field.element_stress.row(1).transpose());
  writer.add_cell_scalars("sigma_xy", field.element_stress.row(2).transpose());
  writer.add_cell_scalars("von_mises", field.element_von_mises);
  writer.add_cell_scalars("principal_max", field.element_principal_max);
  writer.add_cell_scalars("principal_min", field.element_principal_min);
  writer.add_cell_scalars("strain_energy", field.element_strain_energy);
  if (density) writer.add_cell_scalars("density", *density);
  if (stiffness_factor) writer.add_cell_scalars("stiffness_factor", *stiffness_factor);
  writer.write(file("fields_" + sanitise(load_case) + ".vtk"));
}

void ResultWriter::write_modal(const Mesh& mesh, const ModalResult& modal,
                               const std::string& tag) const {
  const std::string suffix = tag.empty() ? "" : "_" + sanitise(tag);
  {
    CsvWriter csv(file("modes" + suffix + ".csv"),
                  {"mode", "eigenvalue[1/s2]", "omega[rad/s]", "frequency[Hz]",
                   "eigenpair_residual[-]", "low_density_ke_fraction[-]"});
    for (Eigen::Index i = 0; i < modal.eigenvalues.size(); ++i) {
      const Scalar frac =
          modal.modal_mass_fraction.size() > i ? modal.modal_mass_fraction(i) : 0.0;
      csv.row(static_cast<Index>(i),
              {modal.eigenvalues(i), modal.angular_frequencies(i),
               modal.frequencies_hz(i), modal.modal_residuals(i), frac});
    }
    csv.close();
  }

  if (!config_.output.write_mode_shapes) return;

  std::vector<std::string> header{"node", "x[m]", "y[m]"};
  for (Eigen::Index i = 0; i < modal.mode_shapes.cols(); ++i) {
    std::ostringstream ux;
    std::ostringstream uy;
    ux << "ux_mode" << i << "[m]";
    uy << "uy_mode" << i << "[m]";
    header.push_back(ux.str());
    header.push_back(uy.str());
  }
  CsvWriter csv(file("mode_shapes" + suffix + ".csv"), header);
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    std::vector<Scalar> row{mesh.node(n).x(), mesh.node(n).y()};
    for (Eigen::Index i = 0; i < modal.mode_shapes.cols(); ++i) {
      row.push_back(modal.mode_shapes(n * kDofsPerNode + 0, i));
      row.push_back(modal.mode_shapes(n * kDofsPerNode + 1, i));
    }
    csv.row(n, row);
  }
  csv.close();

  if (config_.output.write_vtk) {
    for (Eigen::Index i = 0; i < modal.mode_shapes.cols(); ++i) {
      std::ostringstream title;
      title << "SparLab mode " << i << " at " << modal.frequencies_hz(i) << " Hz";
      VtkWriter writer(mesh, title.str());
      writer.add_point_vectors("mode_shape", modal.mode_shapes.col(i));
      Vector magnitude(mesh.num_nodes());
      for (Index n = 0; n < mesh.num_nodes(); ++n) {
        magnitude(n) = std::hypot(modal.mode_shapes(n * kDofsPerNode + 0, i),
                                  modal.mode_shapes(n * kDofsPerNode + 1, i));
      }
      writer.add_point_scalars("mode_shape_magnitude", magnitude);
      std::ostringstream name;
      name << "mode" << suffix << "_" << i << ".vtk";
      writer.write(file(name.str()));
    }
  }
}

void ResultWriter::write_history(const TopologyOptimizationResult& result) const {
  CsvWriter csv(file("history.csv"),
                {"iteration", "penalty[-]", "compliance[J]", "volume[m3]",
                 "volume_fraction[-]", "max_design_change[-]", "lagrange_multiplier[-]",
                 "grey_level[-]", "oc_bisections[-]", "volume_converged[-]",
                 "seconds[s]"});
  for (const TopologyIteration& it : result.history) {
    csv.row(it.iteration,
            {it.penalty, it.compliance, it.volume, it.volume_fraction, it.max_change,
             it.lambda, it.gray_level, static_cast<Scalar>(it.bisections),
             it.volume_converged ? 1.0 : 0.0, it.seconds});
  }
  csv.close();
}

void ResultWriter::write_density(const Mesh& mesh, const DesignDomain& domain,
                                 const TopologyOptimizationResult& result) const {
  CsvWriter csv(file("density_final.csv"),
                {"element", "cx[m]", "cy[m]", "area[m2]", "volume[m3]", "design_x[-]",
                 "physical_density[-]", "stiffness_factor[-]", "passive_tag[-]",
                 "strain_energy[J]"});
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    const Vector2 c = mesh.element_centroid(e);
    csv.row(e, {c.x(), c.y(), mesh.element_area(e), domain.element_volumes()(e),
                result.design(e), result.physical_density(e),
                result.stiffness_factors(e),
                static_cast<Scalar>(static_cast<int>(domain.tags()[static_cast<std::size_t>(e)])),
                result.element_strain_energy.size() > e
                    ? result.element_strain_energy(e)
                    : 0.0});
  }
  csv.close();
}

void ResultWriter::write_density_history(const TopologyOptimizationResult& result,
                                         int max_frames) const {
  if (result.snapshots.empty()) return;
  const int count = static_cast<int>(result.snapshots.size());
  const int stride = std::max(1, (count + max_frames - 1) / std::max(max_frames, 1));

  std::vector<int> frames;
  for (int i = 0; i < count; i += stride) frames.push_back(i);
  if (frames.back() != count - 1) frames.push_back(count - 1);

  const Eigen::Index ne = result.snapshots.front().size();
  std::vector<std::string> header{"iteration"};
  for (Eigen::Index e = 0; e < ne; ++e) {
    std::ostringstream os;
    os << "rho_" << e;
    header.push_back(os.str());
  }
  // Six significant digits keep the animation file small; the final density is
  // written at full precision in density_final.csv.
  CsvWriter csv(file("density_history.csv"), header, 6);
  for (int f : frames) {
    std::vector<Scalar> row;
    row.reserve(static_cast<std::size_t>(ne));
    const Vector& snapshot = result.snapshots[static_cast<std::size_t>(f)];
    for (Eigen::Index e = 0; e < ne; ++e) row.push_back(snapshot(e));
    csv.row(result.snapshot_iterations[static_cast<std::size_t>(f)], row);
  }
  csv.close();
  log::debug("wrote ", frames.size(), " density snapshots (of ", count,
             " recorded) to density_history.csv");
}

// ----------------------------------------------------------------------------
// Summaries
// ----------------------------------------------------------------------------

json::Value make_provenance(const Configuration& config) {
  json::Value out = json::Value::make_object();
  out.set("code", json::Value::make_string("SparLab"));
  out.set("version", json::Value::make_string(SPARLAB_VERSION_STRING));
  out.set("build_type", json::Value::make_string(SPARLAB_BUILD_TYPE));
  out.set("compiler", json::Value::make_string(SPARLAB_COMPILER));
  out.set("eigen_version", json::Value::make_string(SPARLAB_EIGEN_VERSION));
  out.set("config_source", json::Value::make_string(config.source_path));
  out.set("case", json::Value::make_string(config.name));
  out.set("description", json::Value::make_string(config.description));
  out.set("units", json::Value::make_string("SI: m, N, Pa, kg, kg/m^3, Hz, J"));

  const std::time_t now = std::time(nullptr);
  char buffer[64] = {0};
  std::tm utc {};
#if defined(_WIN32)
  gmtime_s(&utc, &now);
#else
  gmtime_r(&now, &utc);
#endif
  std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%SZ", &utc);
  out.set("timestamp_utc", json::Value::make_string(buffer));
  return out;
}

namespace {

json::Value tolerance_json(const Configuration& config) {
  json::Value out = json::Value::make_object();
  out.set("linear_solver", json::Value::make_string(
                               to_string(config.analysis.linear.type)));
  out.set("linear_residual_tolerance",
          json::Value::make_number(config.analysis.linear.residual_tolerance));
  out.set("linear_iterative_tolerance",
          json::Value::make_number(config.analysis.linear.iterative_tolerance));
  out.set("cholesky_pivot_tolerance",
          json::Value::make_number(config.analysis.linear.pivot_tolerance));
  out.set("equilibrium_tolerance",
          json::Value::make_number(config.analysis.equilibrium_tolerance));
  out.set("modal_tolerance", json::Value::make_number(config.modal.options.tolerance));
  out.set("modal_residual_tolerance",
          json::Value::make_number(config.modal.options.residual_tolerance));
  out.set("optimizer_change_tolerance",
          json::Value::make_number(config.topology.optimizer.change_tolerance));
  out.set("optimizer_volume_tolerance",
          json::Value::make_number(config.topology.optimizer.oc.volume_tolerance));
  return out;
}

}  // namespace

json::Value make_static_summary(const Configuration& config, const FemModel& model,
                                const ModelDiagnostics& diagnostics,
                                const std::vector<StaticSolution>& solutions,
                                const std::vector<StressField>& stresses,
                                const ModalResult* modal,
                                const TimingLedger& timings) {
  json::Value out = json::Value::make_object();
  out.set("provenance", make_provenance(config));
  out.set("mesh", mesh_stats_json(model));
  out.set("material", material_json(model.material(), model.stress_state()));
  out.set("tolerances", tolerance_json(config));
  out.set("model_diagnostics", diagnostics_json(diagnostics));

  json::Value cases = json::Value::make_array();
  for (std::size_t l = 0; l < solutions.size(); ++l) {
    const StaticSolution& sol = solutions[l];
    json::Value entry = json::Value::make_object();
    entry.set("name", json::Value::make_string(sol.load_case_name));
    entry.set("weight", json::Value::make_number(sol.weight));
    entry.set("compliance_J", json::Value::make_number(sol.compliance));
    entry.set("strain_energy_J", json::Value::make_number(sol.strain_energy));
    entry.set("compliance_over_twice_strain_energy",
              json::Value::make_number(
                  sol.strain_energy != 0.0
                      ? sol.compliance / (2.0 * sol.strain_energy)
                      : 0.0));
    entry.set("max_displacement_magnitude_m",
              json::Value::make_number(sol.max_displacement_magnitude));
    entry.set("max_displacement_node",
              json::Value::make_number(sol.max_displacement_node));
    entry.set("scaled_residual", json::Value::make_number(sol.scaled_residual));
    entry.set("solver_iterations", json::Value::make_number(sol.solver_iterations));
    entry.set("equilibrium", equilibrium_json(sol.equilibrium));
    if (l < stresses.size()) {
      entry.set("max_von_mises_Pa",
                json::Value::make_number(stresses[l].element_von_mises.maxCoeff()));
      entry.set("total_strain_energy_J",
                json::Value::make_number(stresses[l].element_strain_energy.sum()));
    }
    cases.push_back(entry);
  }
  out.set("load_cases", cases);

  if (modal != nullptr) out.set("modal", modal_json(*modal));
  out.set("timings_s", timings_json(timings));
  return out;
}

json::Value make_topology_summary(const Configuration& config, const FemModel& model,
                                  const DesignDomain& domain,
                                  const DensityFilter& filter,
                                  const TopologyOptimizationResult& result,
                                  const TopologyInterpretation* interpretation,
                                  const ModalResult* modal_initial,
                                  const ModalResult* modal_optimised,
                                  const ModalResult* modal_mass_matched,
                                  const TimingLedger& timings) {
  json::Value out = json::Value::make_object();
  out.set("provenance", make_provenance(config));
  out.set("mesh", mesh_stats_json(model));
  out.set("material", material_json(model.material(), model.stress_state()));
  out.set("tolerances", tolerance_json(config));

  {
    json::Value setup = json::Value::make_object();
    // The objective is a weighted mean over load cases, so the summary has to
    // name the cases and record their weights: `load_case_compliance_J` below
    // is in this order, and without the weights the reported total cannot be
    // reconstructed from the parts. Both the weights as written in the deck and
    // the normalised weights actually used are recorded, because only the
    // latter reproduce the objective:
    //   compliance_J = sum_l load_case_weights_normalised[l] *
    //                        load_case_compliance_J[l].
    {
      std::vector<std::string> names;
      std::vector<Scalar> weights;
      names.reserve(model.load_case_specs().size());
      weights.reserve(model.load_case_specs().size());
      for (const auto& spec : model.load_case_specs()) {
        names.push_back(spec.name);
        weights.push_back(spec.weight);
      }
      setup.set("load_case_names", json::array_of(names));
      setup.set("load_case_weights", json::array_of(weights));
      setup.set("load_case_weights_normalised",
                json::array_of(model.normalised_weights()));
    }
    setup.set("volume_fraction_target",
              json::Value::make_number(domain.volume_fraction()));
    setup.set("volume_target_m3", json::Value::make_number(domain.volume_target()));
    setup.set("domain_volume_m3", json::Value::make_number(domain.domain_volume()));
    setup.set("num_design_variables", json::Value::make_number(domain.num_elements()));
    setup.set("num_free_variables",
              json::Value::make_number(domain.num_free_variables()));
    setup.set("num_passive_solid", json::Value::make_number(domain.num_passive_solid()));
    setup.set("num_passive_void", json::Value::make_number(domain.num_passive_void()));
    setup.set("passive_solid_volume_m3",
              json::Value::make_number(domain.passive_solid_volume()));
    setup.set("passive_void_volume_m3",
              json::Value::make_number(domain.passive_void_volume()));
    setup.set("filter_type", json::Value::make_string(to_string(filter.type())));
    setup.set("filter_radius_m", json::Value::make_number(filter.radius()));
    setup.set("filter_average_support",
              json::Value::make_number(filter.average_support()));
    setup.set("simp_penalty",
              json::Value::make_number(config.topology.optimizer.simp.penalty));
    setup.set("simp_emin_ratio",
              json::Value::make_number(config.topology.optimizer.simp.emin_ratio));
    setup.set("mass_interpolation",
              json::Value::make_string(
                  to_string(config.topology.optimizer.simp.mass_law)));
    setup.set("move_limit",
              json::Value::make_number(config.topology.optimizer.oc.move_limit));
    setup.set("damping", json::Value::make_number(config.topology.optimizer.oc.damping));
    setup.set("continuation_steps",
              json::Value::make_number(config.topology.optimizer.continuation_steps));
    setup.set("penalty_start",
              json::Value::make_number(config.topology.optimizer.penalty_start));
    setup.set("max_iterations",
              json::Value::make_number(config.topology.optimizer.max_iterations));
    setup.set("change_tolerance",
              json::Value::make_number(config.topology.optimizer.change_tolerance));
    setup.set("objective_tolerance",
              json::Value::make_number(config.topology.optimizer.objective_tolerance));
    setup.set("objective_window",
              json::Value::make_number(config.topology.optimizer.objective_window));
    out.set("optimization_setup", setup);
  }

  {
    json::Value res = json::Value::make_object();
    res.set("converged", json::Value::make_bool(result.converged));
    res.set("stop_reason", json::Value::make_string(result.stop_reason));
    res.set("final_design_change",
            json::Value::make_number(result.final_design_change));
    res.set("final_relative_objective_change",
            json::Value::make_number(result.final_objective_change));
    res.set("iterations", json::Value::make_number(result.iterations));
    res.set("linear_solves", json::Value::make_number(result.linear_solves));
    res.set("compliance_J", json::Value::make_number(result.compliance));
    res.set("load_case_compliance_J", json::array_of(result.load_case_compliance));
    res.set("volume_m3", json::Value::make_number(result.volume));
    res.set("volume_fraction", json::Value::make_number(result.volume_fraction));
    res.set("volume_constraint_relative_violation",
            json::Value::make_number(result.volume_constraint_violation));
    res.set("volume_constraint_satisfied",
            json::Value::make_bool(result.volume_constraint_violation <= 1.0e-6));
    res.set("grey_level", json::Value::make_number(result.gray_level));
    res.set("mass_kg",
            json::Value::make_number(result.volume * model.material().density()));
    res.set("total_seconds", json::Value::make_number(result.total_seconds));
    res.set("seconds_per_iteration",
            json::Value::make_number(result.iterations > 0
                                         ? result.total_seconds / result.iterations
                                         : 0.0));
    res.set("warnings", json::array_of(result.warnings));
    if (!result.history.empty()) {
      res.set("initial_compliance_J",
              json::Value::make_number(result.history.front().compliance));
      res.set("compliance_reduction_factor",
              json::Value::make_number(result.compliance > 0.0
                                           ? result.history.front().compliance /
                                                 result.compliance
                                           : 0.0));
    }
    out.set("optimization_result", res);
  }

  if (interpretation != nullptr) {
    json::Value interp = json::Value::make_object();
    interp.set("threshold", json::Value::make_number(interpretation->threshold));
    interp.set("elements_above_threshold",
               json::Value::make_number(interpretation->elements_above_threshold));
    interp.set("elements_retained",
               json::Value::make_number(interpretation->elements_retained));
    interp.set("connected_groups_above_threshold",
               json::Value::make_number(interpretation->components_above_threshold));
    interp.set("volume_above_threshold_m3",
               json::Value::make_number(interpretation->volume_above_threshold));
    interp.set("volume_retained_m3",
               json::Value::make_number(interpretation->volume_retained));
    interp.set("volume_discarded_as_islands_m3",
               json::Value::make_number(interpretation->volume_discarded_as_islands));
    interp.set("note", json::Value::make_string(
                           "A SIMP density field is a relative material distribution. "
                           "This block records how it was turned into a solid body for "
                           "the modal comparison: threshold, largest edge-connected "
                           "group, and the material discarded in the process."));
    out.set("solid_interpretation", interp);
  }

  if (modal_initial != nullptr) out.set("modal_initial_solid", modal_json(*modal_initial));
  if (modal_mass_matched != nullptr) {
    out.set("modal_mass_matched_baseline", modal_json(*modal_mass_matched));
  }
  if (modal_optimised != nullptr) {
    out.set("modal_optimised_topology", modal_json(*modal_optimised));
  }

  out.set("timings_s", timings_json(timings));
  return out;
}

}  // namespace sparlab
