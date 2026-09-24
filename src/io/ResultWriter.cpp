#include "sparlab/io/ResultWriter.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/core/Version.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/io/StlWriter.hpp"
#include "sparlab/io/VtkWriter.hpp"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <map>
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

/// [x, y] on a 2-D model, [x, y, z] on a 3-D one.
json::Value point_json(const Vector3& v, int dim) {
  json::Value out = json::Value::make_array();
  for (int i = 0; i < dim; ++i) out.push_back(json::Value::make_number(v(i)));
  return out;
}

const char* component_name(int k) { return k == 0 ? "x" : k == 1 ? "y" : "z"; }

/// Vector magnitude over the model's components only.
Scalar magnitude(const Vector& full, Index n, int dim) {
  return dim == 2 ? std::hypot(full(n * dim + 0), full(n * dim + 1))
                  : std::hypot(full(n * dim + 0), full(n * dim + 1), full(n * dim + 2));
}

json::Value equilibrium_json(const EquilibriumCheck& eq, int dim) {
  json::Value out = json::Value::make_object();
  out.set("applied_force_N", point_json(eq.applied_force, dim));
  out.set("reaction_force_N", point_json(eq.reaction_force, dim));
  out.set("force_residual_N", point_json(eq.force_residual, dim));
  if (dim == 2) {
    // A plane model has a single moment component, about z.
    out.set("applied_moment_Nm", json::Value::make_number(eq.applied_moment.z()));
    out.set("reaction_moment_Nm", json::Value::make_number(eq.reaction_moment.z()));
    out.set("moment_residual_Nm", json::Value::make_number(eq.moment_residual.z()));
  } else {
    out.set("applied_moment_Nm", point_json(eq.applied_moment, 3));
    out.set("reaction_moment_Nm", point_json(eq.reaction_moment, 3));
    out.set("moment_residual_Nm", point_json(eq.moment_residual, 3));
  }
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
  if (!modal.linear_solver.empty()) {
    out.set("linear_solver", json::Value::make_string(modal.linear_solver));
    out.set("linear_iterations", json::Value::make_number(modal.linear_iterations));
  }
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

json::Value mesh_quality_json(const MeshQuality& q) {
  json::Value out = json::Value::make_object();
  out.set("metric", json::Value::make_string(q.metric));
  out.set("min", json::Value::make_number(q.min));
  out.set("mean", json::Value::make_number(q.mean));
  out.set("worst_element", json::Value::make_number(q.worst_element));
  out.set("poor_threshold", json::Value::make_number(q.poor_threshold));
  out.set("poor_elements", json::Value::make_number(q.poor_elements));
  return out;
}

json::Value count_map_json(const std::map<std::string, Index>& counts) {
  json::Value out = json::Value::make_object();
  for (const auto& entry : counts) out.set(entry.first, json::Value::make_number(entry.second));
  return out;
}

/// What the mesh reader found in a mesh file, so a summary records how the
/// file was interpreted (sets, repairs, warnings) and not only its size.
json::Value mesh_file_json(const Configuration& config) {
  const MeshReadReport& r = config.mesh_report;
  json::Value out = json::Value::make_object();
  out.set("path", json::Value::make_string(config.mesh_file.path));
  out.set("format", json::Value::make_string(r.format));
  if (!r.version.empty()) out.set("version", json::Value::make_string(r.version));
  out.set("scale", json::Value::make_number(r.scale));
  out.set("nodes_in_file", json::Value::make_number(r.nodes_in_file));
  out.set("nodes_used", json::Value::make_number(r.nodes_used));
  out.set("unreferenced_nodes_dropped", json::Value::make_number(r.unreferenced_nodes));
  out.set("elements_in_file", json::Value::make_number(r.elements_in_file));
  out.set("cells", json::Value::make_number(r.cells));
  out.set("boundary_elements", json::Value::make_number(r.boundary_elements));
  out.set("cells_reoriented", json::Value::make_number(r.reoriented));
  out.set("duplicate_nodes", json::Value::make_number(r.duplicate_nodes));
  out.set("duplicates_merged", json::Value::make_bool(r.duplicates_merged));
  out.set("node_sets", count_map_json(r.node_sets));
  out.set("element_sets", count_map_json(r.element_sets));
  out.set("ignored", count_map_json(r.ignored));
  json::Value warnings = json::Value::make_array();
  for (const std::string& w : r.warnings) warnings.push_back(json::Value::make_string(w));
  out.set("warnings", warnings);
  return out;
}

json::Value mesh_stats_json(const Configuration& config, const FemModel& model) {
  const Mesh& mesh = model.mesh();
  const int dim = mesh.dim();
  json::Value out = json::Value::make_object();
  out.set("source", json::Value::make_string(to_string(config.mesh_kind)));
  out.set("element_type", json::Value::make_string(to_string(mesh.element_type())));
  out.set("dim", json::Value::make_number(dim));
  out.set("num_nodes", json::Value::make_number(mesh.num_nodes()));
  out.set("num_elements", json::Value::make_number(mesh.num_elements()));
  out.set("num_dofs", json::Value::make_number(model.dofs().num_dofs()));
  out.set("num_free_dofs", json::Value::make_number(model.dofs().num_free()));
  out.set("num_prescribed_dofs",
          json::Value::make_number(model.dofs().num_constrained()));
  const BoundingBox bb = mesh.bounding_box();
  json::Value box = json::Value::make_array();
  for (int i = 0; i < dim; ++i) box.push_back(json::Value::make_number(bb.lower(i)));
  for (int i = 0; i < dim; ++i) box.push_back(json::Value::make_number(bb.upper(i)));
  out.set(dim == 2 ? "bounding_box_xmin_ymin_xmax_ymax_m"
                   : "bounding_box_xmin_ymin_zmin_xmax_ymax_zmax_m",
          box);
  out.set("thickness_m", json::Value::make_number(model.thickness()));
  out.set("domain_volume_m3", json::Value::make_number(model.domain_volume()));
  if (mesh.structured_info().has_value()) {
    const StructuredGridInfo& info = *mesh.structured_info();
    json::Value grid = json::Value::make_object();
    grid.set("nx", json::Value::make_number(info.nx));
    grid.set("ny", json::Value::make_number(info.ny));
    if (dim == 3) grid.set("nz", json::Value::make_number(info.nz));
    grid.set("lx_m", json::Value::make_number(info.lx));
    grid.set("ly_m", json::Value::make_number(info.ly));
    if (dim == 3) grid.set("lz_m", json::Value::make_number(info.lz));
    grid.set("uniform", json::Value::make_bool(info.uniform));
    out.set("structured_grid", grid);
  } else if (config.mesh_kind == MeshKind::StructuredTri ||
             config.mesh_kind == MeshKind::StructuredTet) {
    // The simplex box meshes split the cells of this grid.
    json::Value grid = json::Value::make_object();
    grid.set("nx", json::Value::make_number(config.mesh_spec.nx));
    grid.set("ny", json::Value::make_number(config.mesh_spec.ny));
    if (dim == 3) grid.set("nz", json::Value::make_number(config.mesh_spec.nz));
    grid.set("lx_m", json::Value::make_number(config.mesh_spec.lx));
    grid.set("ly_m", json::Value::make_number(config.mesh_spec.ly));
    if (dim == 3) grid.set("lz_m", json::Value::make_number(config.mesh_spec.lz));
    grid.set("elements_per_cell", json::Value::make_number(dim == 2 ? 2 : 6));
    out.set("split_grid", grid);
  }
  out.set("mean_element_size_m", json::Value::make_number(mesh.mean_element_size()));
  out.set("quality", mesh_quality_json(mesh.quality()));
  if (config.mesh_kind == MeshKind::File) out.set("file", mesh_file_json(config));
  return out;
}

json::Value material_json(const IsotropicMaterial& m, StressState state) {
  json::Value out = json::Value::make_object();
  out.set("name", json::Value::make_string(m.name()));
  out.set("youngs_modulus_Pa", json::Value::make_number(m.youngs_modulus()));
  out.set("poisson_ratio", json::Value::make_number(m.poisson_ratio()));
  out.set("density_kg_per_m3", json::Value::make_number(m.density()));
  out.set("shear_modulus_Pa", json::Value::make_number(m.shear_modulus()));
  out.set("stress_state", json::Value::make_string(to_string(state)));
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

/// Coordinate column headers for a node or centroid table.
std::vector<std::string> coordinate_headers(int dim, const char* prefix) {
  std::vector<std::string> out;
  for (int k = 0; k < dim; ++k) {
    out.push_back(std::string(prefix) + component_name(k) + "[m]");
  }
  return out;
}

/// Per-component headers such as ux[m], uy[m](, uz[m]).
std::vector<std::string> component_headers(int dim, const char* stem, const char* unit) {
  std::vector<std::string> out;
  for (int k = 0; k < dim; ++k) {
    out.push_back(std::string(stem) + component_name(k) + "[" + unit + "]");
  }
  return out;
}

std::vector<std::string> concat(std::vector<std::string> a,
                                const std::vector<std::string>& b) {
  a.insert(a.end(), b.begin(), b.end());
  return a;
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
  const int dim = mesh.dim();
  json::Value doc = json::Value::make_object();
  doc.set("case", json::Value::make_string(config_.name));
  doc.set("element_type", json::Value::make_string(to_string(mesh.element_type())));
  doc.set("dim", json::Value::make_number(dim));
  doc.set("nodes_per_element", json::Value::make_number(mesh.nodes_per_elem()));

  json::Value nodes = json::Value::make_array();
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    nodes.push_back(point_json(mesh.node(n), dim));
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
    entry.set("node", json::Value::make_number(d / dim));
    entry.set("component", json::Value::make_string(component_name(d % dim)));
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
      bool nonzero = false;
      for (int k = 0; k < dim; ++k) nonzero = nonzero || loads[l](n * dim + k) != 0.0;
      if (!nonzero) continue;
      json::Value f = json::Value::make_object();
      f.set("node", json::Value::make_number(n));
      for (int k = 0; k < dim; ++k) {
        f.set(std::string("f") + component_name(k) + "_N",
              json::Value::make_number(loads[l](n * dim + k)));
      }
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
  const int dim = mesh.dim();
  std::vector<std::string> header{"node"};
  header = concat(header, coordinate_headers(dim, ""));
  header = concat(header, component_headers(dim, "u", "m"));
  header.push_back("umag[m]");
  CsvWriter csv(file("displacement_" + sanitise(load_case) + ".csv"), header);
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    const Vector3 x = mesh.node(n);
    std::vector<Scalar> row;
    for (int k = 0; k < dim; ++k) row.push_back(x(k));
    for (int k = 0; k < dim; ++k) row.push_back(displacement(n * dim + k));
    row.push_back(magnitude(displacement, n, dim));
    csv.row(n, row);
  }
  csv.close();
}

void ResultWriter::write_stress(const Mesh& mesh, const std::string& load_case,
                                const StressField& field, const Vector* density) const {
  const int dim = mesh.dim();
  std::vector<std::string> header;
  if (dim == 2) {
    header = {"element", "cx[m]", "cy[m]", "area[m2]", "density[-]", "exx[-]",
              "eyy[-]", "gxy[-]", "sxx[Pa]", "syy[Pa]", "sxy[Pa]", "von_mises[Pa]",
              "solid_von_mises[Pa]", "principal_max[Pa]", "principal_min[Pa]",
              "strain_energy[J]"};
  } else {
    header = {"element", "cx[m]", "cy[m]", "cz[m]", "volume[m3]", "density[-]",
              "exx[-]", "eyy[-]", "ezz[-]", "gxy[-]", "gyz[-]", "gzx[-]",
              "sxx[Pa]", "syy[Pa]", "szz[Pa]", "sxy[Pa]", "syz[Pa]", "szx[Pa]",
              "von_mises[Pa]", "solid_von_mises[Pa]", "principal_max[Pa]",
              "principal_mid[Pa]", "principal_min[Pa]", "strain_energy[J]"};
  }
  CsvWriter csv(file("stress_" + sanitise(load_case) + ".csv"), header);
  const Eigen::Index nv = field.element_strain.rows();
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    const Vector3 c = mesh.element_centroid(e);
    std::vector<Scalar> row;
    for (int k = 0; k < dim; ++k) row.push_back(c(k));
    row.push_back(mesh.element_measure(e));
    row.push_back(density ? (*density)(e) : 1.0);
    for (Eigen::Index i = 0; i < nv; ++i) row.push_back(field.element_strain(i, e));
    for (Eigen::Index i = 0; i < nv; ++i) row.push_back(field.element_stress(i, e));
    row.push_back(field.element_von_mises(e));
    row.push_back(field.element_solid_von_mises(e));
    row.push_back(field.element_principal_max(e));
    if (dim == 3) row.push_back(field.element_principal_mid(e));
    row.push_back(field.element_principal_min(e));
    row.push_back(field.element_strain_energy(e));
    csv.row(e, row);
  }
  csv.close();
}

void ResultWriter::write_reactions(const Mesh& mesh, const DofManager& dofs,
                                   const std::string& load_case,
                                   const Vector& reactions) const {
  const int dim = mesh.dim();
  std::vector<std::string> header{"node"};
  header = concat(header, coordinate_headers(dim, ""));
  header = concat(header, component_headers(dim, "r", "N"));
  header.push_back("rmag[N]");
  CsvWriter csv(file("reactions_" + sanitise(load_case) + ".csv"), header);
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    bool constrained = false;
    for (int k = 0; k < dim; ++k) constrained = constrained || dofs.is_constrained(n * dim + k);
    if (!constrained) continue;
    const Vector3 x = mesh.node(n);
    std::vector<Scalar> row;
    for (int k = 0; k < dim; ++k) row.push_back(x(k));
    for (int k = 0; k < dim; ++k) row.push_back(reactions(n * dim + k));
    row.push_back(magnitude(reactions, n, dim));
    csv.row(n, row);
  }
  csv.close();
}

void ResultWriter::write_static_vtk(const Mesh& mesh, const std::string& load_case,
                                    const Vector& displacement,
                                    const StressField& field, const Vector* density,
                                    const Vector* stiffness_factor) const {
  const int dim = mesh.dim();
  VtkWriter writer(mesh, "SparLab static solution: " + config_.name + " / " + load_case);
  Vector mag(mesh.num_nodes());
  for (Index n = 0; n < mesh.num_nodes(); ++n) mag(n) = magnitude(displacement, n, dim);
  writer.add_point_vectors("displacement", displacement);
  writer.add_point_scalars("displacement_magnitude", mag);
  writer.add_point_scalars("nodal_von_mises", field.nodal_von_mises);
  if (dim == 2) {
    writer.add_cell_scalars("sigma_xx", field.element_stress.row(0).transpose());
    writer.add_cell_scalars("sigma_yy", field.element_stress.row(1).transpose());
    writer.add_cell_scalars("sigma_xy", field.element_stress.row(2).transpose());
  } else {
    writer.add_cell_scalars("sigma_xx", field.element_stress.row(0).transpose());
    writer.add_cell_scalars("sigma_yy", field.element_stress.row(1).transpose());
    writer.add_cell_scalars("sigma_zz", field.element_stress.row(2).transpose());
    writer.add_cell_scalars("sigma_xy", field.element_stress.row(3).transpose());
    writer.add_cell_scalars("sigma_yz", field.element_stress.row(4).transpose());
    writer.add_cell_scalars("sigma_zx", field.element_stress.row(5).transpose());
  }
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
  const int dim = mesh.dim();
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

  std::vector<std::string> header{"node"};
  header = concat(header, coordinate_headers(dim, ""));
  for (Eigen::Index i = 0; i < modal.mode_shapes.cols(); ++i) {
    for (int k = 0; k < dim; ++k) {
      std::ostringstream name;
      name << "u" << component_name(k) << "_mode" << i << "[m]";
      header.push_back(name.str());
    }
  }
  CsvWriter csv(file("mode_shapes" + suffix + ".csv"), header);
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    const Vector3 x = mesh.node(n);
    std::vector<Scalar> row;
    for (int k = 0; k < dim; ++k) row.push_back(x(k));
    for (Eigen::Index i = 0; i < modal.mode_shapes.cols(); ++i) {
      for (int k = 0; k < dim; ++k) row.push_back(modal.mode_shapes(n * dim + k, i));
    }
    csv.row(n, row);
  }
  csv.close();

  if (config_.output.write_vtk) {
    for (Eigen::Index i = 0; i < modal.mode_shapes.cols(); ++i) {
      std::ostringstream title;
      title << "SparLab mode " << i << " at " << modal.frequencies_hz(i) << " Hz";
      VtkWriter writer(mesh, title.str());
      const Vector shape = modal.mode_shapes.col(i);
      writer.add_point_vectors("mode_shape", shape);
      Vector mag(mesh.num_nodes());
      for (Index n = 0; n < mesh.num_nodes(); ++n) mag(n) = magnitude(shape, n, dim);
      writer.add_point_scalars("mode_shape_magnitude", mag);
      std::ostringstream name;
      name << "mode" << suffix << "_" << i << ".vtk";
      writer.write(file(name.str()));
    }
  }
}

void ResultWriter::write_history(const TopologyOptimizationResult& result) const {
  const bool mma = result.method == OptimizerMethod::MMA;
  std::vector<std::string> header{
      "iteration", "penalty[-]", "compliance[J]", "volume[m3]", "volume_fraction[-]",
      "max_design_change[-]", "lagrange_multiplier[-]", "grey_level[-]",
      mma ? "mma_iterations[-]" : "oc_bisections[-]", "volume_converged[-]",
      "seconds[s]"};
  if (mma) {
    header.insert(header.end(), {"max_stress_ratio[-]", "stress_constraint[-]",
                                 "constraint_violation[-]"});
  }
  header.push_back("linear_iterations[-]");
  if (result.projected) header.push_back("beta[-]");
  CsvWriter csv(file("history.csv"), header);
  for (const TopologyIteration& it : result.history) {
    std::vector<Scalar> row{it.penalty, it.compliance, it.volume, it.volume_fraction,
                            it.max_change, it.lambda, it.gray_level,
                            static_cast<Scalar>(it.bisections),
                            it.volume_converged ? 1.0 : 0.0, it.seconds};
    if (mma) {
      row.insert(row.end(), {it.max_stress_ratio, it.stress_constraint,
                             it.constraint_violation});
    }
    row.push_back(static_cast<Scalar>(it.linear_iterations));
    if (result.projected) row.push_back(it.beta);
    csv.row(it.iteration, row);
  }
  csv.close();
}

void ResultWriter::write_density(const Mesh& mesh, const DesignDomain& domain,
                                 const TopologyOptimizationResult& result) const {
  const int dim = mesh.dim();
  std::vector<std::string> header{"element"};
  header = concat(header, coordinate_headers(dim, "c"));
  if (dim == 2) header.push_back("area[m2]");
  header = concat(header, {"volume[m3]", "design_x[-]", "physical_density[-]",
                           "stiffness_factor[-]", "passive_tag[-]", "strain_energy[J]"});
  const bool projected =
      result.projected && result.filtered_density.size() == mesh.num_elements();
  if (projected) header.push_back("filtered_density[-]");
  CsvWriter csv(file("density_final.csv"), header);
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    const Vector3 c = mesh.element_centroid(e);
    std::vector<Scalar> row;
    for (int k = 0; k < dim; ++k) row.push_back(c(k));
    if (dim == 2) row.push_back(mesh.element_measure(e));
    row.push_back(domain.element_volumes()(e));
    row.push_back(result.design(e));
    row.push_back(result.physical_density(e));
    row.push_back(result.stiffness_factors(e));
    row.push_back(static_cast<Scalar>(
        static_cast<int>(domain.tags()[static_cast<std::size_t>(e)])));
    row.push_back(result.element_strain_energy.size() > e
                      ? result.element_strain_energy(e)
                      : 0.0);
    if (projected) row.push_back(result.filtered_density(e));
    csv.row(e, row);
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

json::Value ResultWriter::write_geometry(const Mesh& mesh, const Vector& density,
                                        Scalar thickness, const std::string& stem,
                                        const std::string& what) const {
  if (density.size() != mesh.num_elements()) {
    std::ostringstream os;
    os << "geometry export '" << stem << "': density has " << density.size()
       << " entries but the mesh has " << mesh.num_elements() << " cells";
    throw IoError(os.str());
  }
  const std::string vtk_name = stem + ".vtk";
  const std::string stl_name = stem + ".stl";

  VtkWriter vtk(mesh, "SparLab " + what + ": " + config_.name);
  vtk.add_cell_scalars("density", density);
  vtk.write(file(vtk_name));

  const Scalar extrusion = mesh.dim() == 2 ? thickness : 1.0;
  const TriangleSurface surface = boundary_surface(mesh, extrusion);
  write_stl(file(stl_name), surface, config_.name + " " + what);
  const SurfaceStats stats = surface_stats(surface);

  // Closure is checked exactly on the edge pairing. The enclosed volume is
  // compared with the cell volume as well: identical for planar faces, and
  // different by the flat-triangle approximation of a bilinear patch where a
  // distorted cell was cut, which is reported so it is never mistaken for
  // exactness.
  Scalar cell_volume = 0.0;
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    cell_volume += mesh.element_measure(e) * extrusion;
  }
  const Scalar mismatch =
      std::abs(stats.enclosed_volume - cell_volume) / std::max(cell_volume, 1.0e-300);
  if (!stats.closed) {
    log::warn("geometry export '", stem, "': the STL surface has ", stats.unmatched_edges,
              " unmatched directed edges, so it is not closed or not consistently "
              "oriented; do not use it as a solid");
  } else if (stats.enclosed_volume <= 0.0) {
    log::warn("geometry export '", stem, "': the STL surface encloses a non-positive "
              "volume (", stats.enclosed_volume, " m^3); its normals point inward");
  } else if (stats.non_manifold_edges > 0) {
    log::warn("geometry export '", stem, "': the STL surface is closed but has ",
              stats.non_manifold_edges, " non-manifold edge(s) where cells touch only "
              "along an edge; a slicer may split it into several shells, and the "
              "interpretation threshold or connectivity rule decides whether such "
              "cells are one part");
  }

  json::Value out = json::Value::make_object();
  out.set("vtk", json::Value::make_string(vtk_name));
  out.set("stl", json::Value::make_string(stl_name));
  out.set("stl_format", json::Value::make_string("binary, outward normals"));
  out.set("num_cells", json::Value::make_number(mesh.num_elements()));
  out.set("num_triangles", json::Value::make_number(stats.num_triangles));
  out.set("closed_surface", json::Value::make_bool(stats.closed));
  out.set("unmatched_edges", json::Value::make_number(stats.unmatched_edges));
  out.set("non_manifold_edges", json::Value::make_number(stats.non_manifold_edges));
  out.set("surface_area_m2", json::Value::make_number(stats.area));
  out.set("enclosed_volume_m3", json::Value::make_number(stats.enclosed_volume));
  out.set("cell_volume_m3", json::Value::make_number(cell_volume));
  out.set("volume_relative_mismatch", json::Value::make_number(mismatch));
  json::Value lower = json::Value::make_array();
  json::Value upper = json::Value::make_array();
  for (int i = 0; i < 3; ++i) {
    lower.push_back(json::Value::make_number(stats.bounds.lower(i)));
    upper.push_back(json::Value::make_number(stats.bounds.upper(i)));
  }
  out.set("bounds_lower_m", lower);
  out.set("bounds_upper_m", upper);
  if (mesh.dim() == 2) {
    out.set("extruded_thickness_m", json::Value::make_number(extrusion));
  }
  log::debug("wrote ", what, ": ", mesh.num_elements(), " cells, ", stats.num_triangles,
             " triangles, ", stats.enclosed_volume, " m^3");
  return out;
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

json::Value amg_stats_json(const AmgStats& st) {
  json::Value out = json::Value::make_object();
  out.set("levels", json::Value::make_number(static_cast<Scalar>(st.levels.size())));
  json::Value levels = json::Value::make_array();
  for (const AmgLevelStats& l : st.levels) {
    json::Value lv = json::Value::make_object();
    lv.set("unknowns", json::Value::make_number(l.unknowns));
    lv.set("nonzeros", json::Value::make_number(l.nonzeros));
    lv.set("lambda_max_Dinv_A", json::Value::make_number(l.lambda_max));
    levels.push_back(lv);
  }
  out.set("level_sizes", levels);
  out.set("operator_complexity", json::Value::make_number(st.operator_complexity));
  out.set("grid_complexity", json::Value::make_number(st.grid_complexity));
  out.set("near_null_space_dimension",
          json::Value::make_number(st.near_null_space_dimension));
  out.set("coarse_pivot_ratio", json::Value::make_number(st.coarse_pivot_ratio));
  out.set("last_setup_seconds", json::Value::make_number(st.setup_seconds));
  out.set("aggregates_reused_in_last_setup", json::Value::make_bool(st.reused_aggregates));
  return out;
}

json::Value linear_solver_json(const LinearSolverOptions& o) {
  json::Value out = json::Value::make_object();
  out.set("type", json::Value::make_string(to_string(o.type)));
  out.set("warm_start", json::Value::make_bool(o.warm_start));
  if (o.type == LinearSolverType::Auto) {
    out.set("auto_direct_limit_plane", json::Value::make_number(o.auto_direct_limit_2d));
    out.set("auto_direct_limit_solid", json::Value::make_number(o.auto_direct_limit_3d));
  }
  if (o.type == LinearSolverType::AmgCg || o.type == LinearSolverType::Auto) {
    json::Value amg = json::Value::make_object();
    amg.set("strength_threshold", json::Value::make_number(o.amg.strength_threshold));
    amg.set("max_levels", json::Value::make_number(o.amg.max_levels));
    amg.set("coarse_size", json::Value::make_number(o.amg.coarse_size));
    amg.set("smoother", json::Value::make_string(to_string(o.amg.smoother)));
    amg.set("smoother_degree", json::Value::make_number(o.amg.smoother_degree));
    amg.set("chebyshev_ratio", json::Value::make_number(o.amg.chebyshev_ratio));
    amg.set("prolongator_damping", json::Value::make_number(o.amg.prolongator_damping));
    amg.set("lanczos_steps", json::Value::make_number(o.amg.lanczos_steps));
    amg.set("reuse_aggregates", json::Value::make_bool(o.amg.reuse_aggregates));
    amg.set("coarse_pivot_tolerance",
            json::Value::make_number(o.amg.coarse_pivot_tolerance));
    out.set("multigrid", amg);
  }
  return out;
}

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
  out.set("linear_solver_settings", linear_solver_json(config.analysis.linear));
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
  out.set("mesh", mesh_stats_json(config, model));
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
    entry.set("linear_solver", json::Value::make_string(sol.solver_name));
    entry.set("equilibrium", equilibrium_json(sol.equilibrium, model.dim()));
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
  out.set("mesh", mesh_stats_json(config, model));
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
    setup.set("initial_density", json::Value::make_number(domain.initial_density()));
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
    setup.set("method", json::Value::make_string(to_string(result.method)));
    if (result.method == OptimizerMethod::MMA) {
      json::Value mma = json::Value::make_object();
      const MmaOptions& mo = config.topology.optimizer.mma;
      mma.set("move_limit", json::Value::make_number(mo.move_limit));
      mma.set("asymptote_init", json::Value::make_number(mo.asymptote_init));
      mma.set("asymptote_increase", json::Value::make_number(mo.asymptote_increase));
      mma.set("asymptote_decrease", json::Value::make_number(mo.asymptote_decrease));
      mma.set("constraint_penalty", json::Value::make_number(mo.c));
      mma.set("subproblem_tolerance", json::Value::make_number(mo.epsimin));
      mma.set("constraint_tolerance",
              json::Value::make_number(config.topology.optimizer.constraint_tolerance));
      setup.set("mma", mma);
    }
    {
      const ProjectionOptions& po = config.topology.optimizer.projection;
      json::Value pj = json::Value::make_object();
      pj.set("enabled", json::Value::make_bool(po.enabled));
      if (po.enabled) {
        pj.set("eta", json::Value::make_number(po.eta));
        pj.set("beta_start", json::Value::make_number(po.beta_start));
        pj.set("beta_max", json::Value::make_number(po.beta_max));
        pj.set("beta_factor", json::Value::make_number(po.beta_factor));
        pj.set("beta_interval", json::Value::make_number(po.beta_interval));
        pj.set("advance_on_convergence", json::Value::make_bool(po.advance_on_convergence));
        pj.set("formulation",
               json::Value::make_string(
                   "rho_bar = (tanh(beta eta) + tanh(beta (rho_tilde - eta))) / "
                   "(tanh(beta eta) + tanh(beta (1 - eta))) on the filtered density; "
                   "SIMP, volume and stress act on rho_bar"));
      }
      setup.set("projection", pj);
    }
    if (result.stress_constrained) {
      json::Value sc = json::Value::make_object();
      const StressConstraintOptions& so = config.topology.optimizer.stress;
      sc.set("limit_Pa", json::Value::make_number(so.limit));
      sc.set("p_norm", json::Value::make_number(so.p_norm));
      sc.set("relaxation_q", json::Value::make_number(so.relaxation));
      sc.set("scaling_blend", json::Value::make_number(so.scaling_blend));
      sc.set("feasibility_tolerance", json::Value::make_number(so.feasibility_tolerance));
      sc.set("formulation", json::Value::make_string(
                                "rho^q sigma_vm(solid) at the element centre, p-norm over "
                                "elements with adaptive scale c, one constraint per load "
                                "case: c * g_PN - 1 <= 0"));
      setup.set("stress_constraint", sc);
    }
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
    {
      json::Value ls = json::Value::make_object();
      ls.set("solver", json::Value::make_string(result.linear_solver));
      ls.set("iterative_iterations_total",
             json::Value::make_number(result.linear_iterations));
      ls.set("iterative_iterations_per_solve",
             json::Value::make_number(result.linear_solves > 0
                                          ? static_cast<Scalar>(result.linear_iterations) /
                                                static_cast<Scalar>(result.linear_solves)
                                          : 0.0));
      if (result.has_amg_stats) ls.set("multigrid", amg_stats_json(result.amg_stats));
      res.set("linear_solver", ls);
    }
    if (result.projected) {
      json::Value pj = json::Value::make_object();
      pj.set("final_beta", json::Value::make_number(result.final_beta));
      pj.set("eta", json::Value::make_number(result.projection_eta));
      pj.set("filtered_grey_level",
             json::Value::make_number(gray_level(result.filtered_density)));
      pj.set("note", json::Value::make_string(
                         "grey_level is that of the projected (physical) density; "
                         "filtered_grey_level that of the density before projection"));
      res.set("projection", pj);
    }
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
    if (result.method == OptimizerMethod::MMA) {
      res.set("constraint_violation", json::Value::make_number(result.constraint_violation));
      res.set("feasible", json::Value::make_bool(result.feasible));
    }
    if (result.stress_constrained) {
      json::Value stress = json::Value::make_object();
      stress.set("max_relaxed_stress_ratio",
                 json::Value::make_number(result.max_stress_ratio));
      stress.set("max_relaxed_stress_Pa",
                 json::Value::make_number(result.max_stress_ratio *
                                          config.topology.optimizer.stress.limit));
      json::Value cases = json::Value::make_array();
      for (const StressConstraintRecord& rec : result.stress) {
        json::Value c = json::Value::make_object();
        c.set("load_case", json::Value::make_string(rec.load_case));
        c.set("max_relaxed_stress_Pa", json::Value::make_number(rec.max_relaxed_stress_Pa));
        c.set("max_relaxed_stress_ratio", json::Value::make_number(rec.max_relaxed_ratio));
        c.set("max_element", json::Value::make_number(rec.max_element));
        c.set("p_norm_ratio", json::Value::make_number(rec.p_norm_ratio));
        c.set("scale", json::Value::make_number(rec.scale));
        c.set("constraint_value", json::Value::make_number(rec.constraint));
        c.set("max_solid_stress_ratio_retained",
              json::Value::make_number(rec.max_solid_ratio_retained));
        cases.push_back(c);
      }
      stress.set("load_cases", cases);
      stress.set("note", json::Value::make_string(
                             "The constraint bounds the relaxed (rho^q) aggregated stress "
                             "of the SIMP model. max_solid_stress_ratio_retained is the "
                             "unrelaxed solid-material stress over the elements at or above "
                             "the interpretation threshold, and "
                             "interpreted_solid_analysis.max_von_mises_over_limit is the "
                             "re-solve of the thresholded structure: those are the numbers "
                             "that say whether the structure meets the limit."));
      res.set("stress", stress);
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
