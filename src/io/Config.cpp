#include "sparlab/io/Config.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <cmath>
#include <filesystem>
#include <sstream>

namespace sparlab {
namespace {

StressState parse_stress_state(const std::string& text) {
  if (text == "plane_stress") return StressState::PlaneStress;
  if (text == "plane_strain") return StressState::PlaneStrain;
  if (text == "three_dimensional" || text == "3d") return StressState::ThreeDimensional;
  throw ConfigError("unknown stress state '" + text +
                    "' (expected plane_stress|plane_strain|three_dimensional)");
}

MassType parse_mass_type(const std::string& text) {
  if (text == "consistent") return MassType::Consistent;
  if (text == "lumped") return MassType::Lumped;
  throw ConfigError("unknown mass type '" + text + "' (expected consistent|lumped)");
}

int parse_axis(const ConfigNode& node, const std::string& key) {
  const std::string text = node.string_or(key, "z");
  if (text == "x") return 0;
  if (text == "y") return 1;
  if (text == "z") return 2;
  throw ConfigError("'" + node.path() + "." + key + "' must be \"x\", \"y\" or \"z\", got \"" +
                    text + "\"");
}

Selector parse_selector_primitive(const ConfigNode& node, int dim) {
  Selector sel;
  sel.tolerance = node.number_or("tolerance", 0.0);

  int matched = 0;
  if (node.child("all").exists()) {
    sel.kind = SelectorKind::All;
    ++matched;
  }
  if (node.child("box").exists()) {
    const ConfigNode box = node.child("box");
    sel.kind = SelectorKind::Box;
    sel.xmin = box.number_or("xmin", -std::numeric_limits<Scalar>::infinity());
    sel.xmax = box.number_or("xmax", std::numeric_limits<Scalar>::infinity());
    sel.ymin = box.number_or("ymin", -std::numeric_limits<Scalar>::infinity());
    sel.ymax = box.number_or("ymax", std::numeric_limits<Scalar>::infinity());
    sel.zmin = box.number_or("zmin", -std::numeric_limits<Scalar>::infinity());
    sel.zmax = box.number_or("zmax", std::numeric_limits<Scalar>::infinity());
    if (dim == 2 && (box.child("zmin").exists() || box.child("zmax").exists())) {
      throw ConfigError("'" + box.path() +
                        "' sets zmin/zmax but the mesh is two-dimensional");
    }
    if (sel.xmin > sel.xmax || sel.ymin > sel.ymax || sel.zmin > sel.zmax) {
      std::ostringstream os;
      os << "'" << box.path() << "' has an empty interval: x in [" << sel.xmin << ", "
         << sel.xmax << "], y in [" << sel.ymin << ", " << sel.ymax << "]";
      if (dim == 3) os << ", z in [" << sel.zmin << ", " << sel.zmax << "]";
      throw ConfigError(os.str());
    }
    ++matched;
  }
  if (node.child("circle").exists()) {
    const ConfigNode c = node.child("circle");
    sel.kind = SelectorKind::Circle;
    sel.center = c.require("center").vector3(dim);
    sel.radius = c.positive_number("radius");
    sel.axis = parse_axis(c, "axis");
    if (dim == 2 && sel.axis != 2) {
      throw ConfigError("'" + c.path() + ".axis' must be \"z\" on a two-dimensional mesh");
    }
    ++matched;
  }
  if (node.child("annulus").exists()) {
    const ConfigNode c = node.child("annulus");
    sel.kind = SelectorKind::Annulus;
    sel.center = c.require("center").vector3(dim);
    sel.inner_radius = c.number_or("inner_radius", 0.0);
    sel.radius = c.positive_number("radius");
    sel.axis = parse_axis(c, "axis");
    if (dim == 2 && sel.axis != 2) {
      throw ConfigError("'" + c.path() + ".axis' must be \"z\" on a two-dimensional mesh");
    }
    if (sel.inner_radius >= sel.radius) {
      std::ostringstream os;
      os << "'" << c.path() << "' needs inner_radius (" << sel.inner_radius
         << " m) strictly below radius (" << sel.radius << " m)";
      throw ConfigError(os.str());
    }
    ++matched;
  }
  if (node.child("sphere").exists()) {
    const ConfigNode c = node.child("sphere");
    sel.kind = SelectorKind::Sphere;
    sel.center = c.require("center").vector3(dim);
    sel.radius = c.positive_number("radius");
    ++matched;
  }
  if (node.child("node_ids").exists()) {
    sel.kind = SelectorKind::NodeIds;
    sel.ids = node.child("node_ids").index_list();
    if (sel.ids.empty()) {
      throw ConfigError("'" + node.child("node_ids").path() + "' is empty");
    }
    ++matched;
  }
  if (node.child("element_ids").exists()) {
    sel.kind = SelectorKind::ElementIds;
    sel.ids = node.child("element_ids").index_list();
    if (sel.ids.empty()) {
      throw ConfigError("'" + node.child("element_ids").path() + "' is empty");
    }
    ++matched;
  }
  if (node.child("nearest_node").exists()) {
    sel.kind = SelectorKind::NearestNode;
    sel.point = node.child("nearest_node").vector3(dim);
    ++matched;
  }
  if (node.child("group").exists()) {
    sel.kind = SelectorKind::Group;
    sel.group = node.child("group").string();
    if (sel.group.empty()) {
      throw ConfigError("'" + node.child("group").path() + "' is an empty group name");
    }
    ++matched;
  }

  if (matched == 0) {
    throw ConfigError(
        "'" + node.path() +
        "' does not name a region primitive; expected one of all, box, circle, "
        "annulus, sphere, node_ids, element_ids, nearest_node, group (or an 'any_of' "
        "list)");
  }
  if (matched > 1) {
    throw ConfigError("'" + node.path() +
                      "' names more than one region primitive; use 'any_of' to combine "
                      "several regions");
  }
  return sel;
}

}  // namespace

std::string to_string(MeshKind kind) {
  switch (kind) {
    case MeshKind::StructuredQuad: return "structured_quad";
    case MeshKind::StructuredHex: return "structured_hex";
    case MeshKind::StructuredTri: return "structured_tri";
    case MeshKind::StructuredTet: return "structured_tet";
    case MeshKind::File: return "file";
  }
  return "unknown";
}

bool is_structured(MeshKind kind) { return kind != MeshKind::File; }

SelectorGroup parse_region(const ConfigNode& node, const std::string& default_name,
                           int dim) {
  if (!node.exists()) {
    throw ConfigError("required region '" + node.path() + "' is missing");
  }
  SelectorGroup group;
  group.name = node.string_or("name", default_name);
  group.invert = node.boolean_or("invert", false);

  const std::vector<ConfigNode> any_of = node.array("any_of");
  if (!any_of.empty()) {
    for (const ConfigNode& item : any_of) {
      group.members.push_back(parse_selector_primitive(item, dim));
    }
  } else {
    group.members.push_back(parse_selector_primitive(node, dim));
  }
  return group;
}

const IsotropicMaterial& Configuration::material() const {
  if (!material_.has_value()) {
    throw ConfigError("the configuration has no material section");
  }
  return *material_;
}

void Configuration::set_material(IsotropicMaterial material) {
  material_ = std::move(material);
}

int Configuration::dim() const {
  switch (mesh_kind) {
    case MeshKind::StructuredQuad:
    case MeshKind::StructuredTri: return 2;
    case MeshKind::StructuredHex:
    case MeshKind::StructuredTet: return 3;
    case MeshKind::File:
      if (file_mesh == nullptr) {
        throw ConfigError("the configuration names a mesh file that has not been read");
      }
      return file_mesh->dim();
  }
  throw ConfigError("unhandled mesh type");
}

std::string Configuration::describe_mesh() const {
  std::ostringstream os;
  if (mesh_kind == MeshKind::File) {
    os << "file '" << mesh_file.path << "'";
    if (file_mesh != nullptr) {
      os << " (" << (mesh_elevated ? "Tet4 elevated to " : "")
         << to_string(file_mesh->element_type()) << ", " << file_mesh->num_elements()
         << " elements)";
    }
    return os.str();
  }
  os << to_string(mesh_kind) << " " << mesh_spec.nx << " x " << mesh_spec.ny;
  if (dim() == 3) os << " x " << mesh_spec.nz;
  if (mesh_order == 2) os << " (Tet10)";
  return os.str();
}

Scalar Configuration::resolved_filter_radius(const Mesh& mesh) const {
  if (topology.filter_radius > 0.0) return topology.filter_radius;
  if (!(topology.filter_radius_elements > 0.0)) {
    throw ConfigError(
        "topology.filter needs either 'radius' (metres) or a positive "
        "'radius_elements'");
  }
  // The side of the square / cube with the mean cell measure for Q4 and
  // Hex8, the mean edge length for triangles and tetrahedra (Mesh.hpp).
  return topology.filter_radius_elements * mesh.mean_element_size();
}

std::vector<std::size_t> Configuration::buckling_load_cases() const {
  std::vector<std::size_t> out;
  if (buckling.load_cases.empty()) {
    for (std::size_t l = 0; l < load_cases.size(); ++l) out.push_back(l);
    return out;
  }
  for (const std::string& wanted : buckling.load_cases) {
    std::size_t l = 0;
    while (l < load_cases.size() && load_cases[l].name != wanted) ++l;
    if (l == load_cases.size()) {
      throw ConfigError("'buckling.load_cases' names '" + wanted +
                        "', which is not a load case of the deck");
    }
    out.push_back(l);
  }
  return out;
}

namespace {

/// The eigensolver keys shared by the buckling check and the constraint.
void parse_buckling_solver(const ConfigNode& node, BucklingOptions& options) {
  options.max_iterations = node.integer_or("max_iterations", options.max_iterations);
  options.tolerance = node.number_or("tolerance", options.tolerance);
  options.residual_tolerance = node.number_or("residual_tolerance", options.residual_tolerance);
  options.seed = static_cast<unsigned int>(node.integer_or("seed", static_cast<int>(options.seed)));
}

std::vector<std::string> string_list(const ConfigNode& parent, const std::string& key) {
  std::vector<std::string> out;
  for (const ConfigNode& item : parent.array(key)) out.push_back(item.string());
  return out;
}

}  // namespace

Configuration parse_configuration(const json::Value& document, const std::string& source,
                                  bool strict, const std::string& base_directory) {
  if (!document.is_object()) {
    throw ConfigError(source + ": the top-level configuration must be a JSON object");
  }
  const ConfigNode root(document, "");

  Configuration config;
  config.document = document;
  config.source_path = source;
  config.name = root.string_or("name", "case");
  config.description = root.string_or("description", "");
  // Informational only, but read so that it is not flagged as an unused key.
  (void)root.string_or("units", "SI");

  // --- mesh ---------------------------------------------------------------
  {
    const ConfigNode mesh = root.require("mesh");
    const std::string type = mesh.string_or("type", "structured_quad");
    if (type == "structured_quad") {
      config.mesh_kind = MeshKind::StructuredQuad;
    } else if (type == "structured_hex") {
      config.mesh_kind = MeshKind::StructuredHex;
    } else if (type == "structured_tri") {
      config.mesh_kind = MeshKind::StructuredTri;
    } else if (type == "structured_tet") {
      config.mesh_kind = MeshKind::StructuredTet;
    } else if (type == "file") {
      config.mesh_kind = MeshKind::File;
    } else {
      throw ConfigError("'" + mesh.path() + ".type' must be one of structured_quad, "
                        "structured_tri, structured_hex, structured_tet or file; got '" +
                        type + "'");
    }
    if (config.mesh_kind == MeshKind::File) {
      for (const char* key : {"nx", "ny", "nz", "lx", "ly", "lz", "x0", "y0", "z0"}) {
        if (mesh.child(key).exists()) {
          throw ConfigError("'" + mesh.path() + "." + key +
                            "' does not apply to a mesh read from a file: its resolution "
                            "and extent come from the mesher. Remove the key, or refine "
                            "the mesh in the generator");
        }
      }
      MeshFileConfig& file = config.mesh_file;
      file.path = mesh.require("path").string();
      std::filesystem::path resolved(file.path);
      if (resolved.is_relative() && !base_directory.empty()) {
        resolved = std::filesystem::path(base_directory) / resolved;
      }
      file.resolved_path = resolved.lexically_normal().string();
      file.read.format = mesh.string_or("format", "auto");
      file.read.scale = mesh.number_or("scale", 1.0);
      if (!(file.read.scale > 0.0) || !std::isfinite(file.read.scale)) {
        std::ostringstream os;
        os << "'" << mesh.path() << ".scale' must be a positive factor (e.g. 0.001 for "
           << "a mesh in millimetres), got " << file.read.scale;
        throw ConfigError(os.str());
      }
      file.read.merge_duplicate_nodes = mesh.boolean_or("merge_duplicate_nodes", false);
      const ConfigNode order_node = mesh.child("order");
      file.read.duplicate_tolerance = mesh.number_or("duplicate_tolerance", 0.0);
      // The reader's messages name the file and line; the deck key is added
      // in front, and the error keeps its category (I/O versus mesh).
      const std::string where = "'" + mesh.path() + ".path' = \"" + file.path + "\": ";
      try {
        config.file_mesh = std::make_shared<const Mesh>(
            read_mesh_file(file.resolved_path, file.read, &config.mesh_report));
      } catch (const IoError& e) {
        throw IoError(where + e.what());
      } catch (const MeshError& e) {
        throw MeshError(where + e.what());
      } catch (const ConfigError& e) {
        throw ConfigError(where + e.what());
      }
      const ElementType read_type = config.file_mesh->element_type();
      config.mesh_order = read_type == ElementType::Tet10 ? 2 : 1;
      if (order_node.exists()) {
        const long long order = order_node.integer();
        if (order != 1 && order != 2) {
          throw ConfigError("'" + order_node.path() + "' must be 1 or 2, got " +
                            std::to_string(order));
        }
        if (order == 1 && read_type == ElementType::Tet10) {
          throw ConfigError("'" + order_node.path() + "' is 1 but " + file.path +
                            " holds 10-node tetrahedra; SparLab does not drop edge "
                            "nodes. Remove the key, or re-export the mesh with linear "
                            "elements");
        }
        if (order == 2 && read_type != ElementType::Tet4 && read_type != ElementType::Tet10) {
          throw ConfigError("'" + order_node.path() + "' is 2 but " + file.path +
                            " holds " + to_string(read_type) +
                            " cells; the quadratic element is the 10-node tetrahedron, "
                            "so order 2 needs a tetrahedral mesh");
        }
        if (order == 2 && read_type == ElementType::Tet4) {
          config.file_mesh = std::make_shared<const Mesh>(elevate_to_tet10(*config.file_mesh));
          config.mesh_elevated = true;
        }
        config.mesh_order = static_cast<int>(order);
      }
    } else {
      const bool solid = config.mesh_kind == MeshKind::StructuredHex ||
                         config.mesh_kind == MeshKind::StructuredTet;
      config.mesh_spec.nx = mesh.require("nx").integer();
      config.mesh_spec.ny = mesh.require("ny").integer();
      config.mesh_spec.lx = mesh.positive_number("lx");
      config.mesh_spec.ly = mesh.positive_number("ly");
      config.mesh_spec.x0 = mesh.number_or("x0", 0.0);
      config.mesh_spec.y0 = mesh.number_or("y0", 0.0);
      if (solid) {
        config.mesh_spec.nz = mesh.require("nz").integer();
        config.mesh_spec.lz = mesh.positive_number("lz");
        config.mesh_spec.z0 = mesh.number_or("z0", 0.0);
      }
      const ConfigNode order_node = mesh.child("order");
      if (order_node.exists()) {
        const long long order = order_node.integer();
        if (order != 1 && order != 2) {
          throw ConfigError("'" + order_node.path() + "' must be 1 or 2, got " +
                            std::to_string(order));
        }
        if (order == 2 && config.mesh_kind != MeshKind::StructuredTet) {
          throw ConfigError("'" + order_node.path() + "' is 2 but the quadratic element "
                            "is the 10-node tetrahedron; use \"type\": "
                            "\"structured_tet\" (or a tetrahedral mesh file) for order 2");
        }
        config.mesh_order = static_cast<int>(order);
      }
      if (!solid && (mesh.child("nz").exists() || mesh.child("lz").exists() ||
                     mesh.child("z0").exists())) {
        throw ConfigError("'" + mesh.path() + "' gives nz/lz/z0 for a " + type +
                          " mesh; use \"type\": \"structured_hex\" or "
                          "\"structured_tet\" for a solid mesh");
      }
    }
  }
  const int dim = config.dim();

  // --- material -----------------------------------------------------------
  {
    const ConfigNode mat = root.require("material");
    config.set_material(IsotropicMaterial(mat.positive_number("youngs_modulus"),
                                          mat.require("poisson_ratio").number(),
                                          mat.number_or("density", 0.0),
                                          mat.string_or("name", "material")));
  }

  // --- model --------------------------------------------------------------
  {
    const ConfigNode model = root.child("model");
    // A unit out-of-plane thickness is the conventional default for a 2-D
    // plane problem; it applies whether or not a "model" section is present.
    // A solid mesh has no thickness and rejects any other value.
    config.thickness = model.number_or("thickness", 1.0);
    if (!(config.thickness > 0.0)) {
      std::ostringstream os;
      os << "'model.thickness' must be positive, got " << config.thickness << " m";
      throw ConfigError(os.str());
    }
    if (dim == 3 && config.thickness != 1.0) {
      std::ostringstream os;
      os << "'model.thickness' is " << config.thickness << " m but the "
         << config.describe_mesh() << " mesh is a solid with no thickness; remove the key";
      throw ConfigError(os.str());
    }
    config.stress_state = parse_stress_state(
        model.string_or("stress_state", dim == 3 ? "three_dimensional" : "plane_stress"));
    if (stress_state_dimension(config.stress_state) != dim) {
      std::ostringstream os;
      os << "'model.stress_state' = \"" << to_string(config.stress_state) << "\" is a "
         << stress_state_dimension(config.stress_state) << "-D idealisation but the "
         << config.describe_mesh() << " mesh is " << dim << "-D";
      throw ConfigError(os.str());
    }
    const ConfigNode integ = model.child("integration");
    config.integration.stiffness_points = integ.integer_or("stiffness_points", 2);
    config.integration.mass_points = integ.integer_or("mass_points", 3);
    // "face_points" is the dimension-neutral name; "edge_points" is kept for
    // the existing 2-D decks.
    config.integration.edge_points =
        integ.integer_or("face_points", integ.integer_or("edge_points", 2));
  }

  // --- boundary conditions ------------------------------------------------
  {
    const std::vector<ConfigNode> bcs = root.array("boundary_conditions");
    if (bcs.empty()) {
      throw ConfigError(
          "'boundary_conditions' is missing or empty; an unconstrained model has a "
          "singular stiffness matrix");
    }
    int index = 0;
    for (const ConfigNode& bc : bcs) {
      DisplacementConstraint constraint;
      std::ostringstream default_name;
      default_name << "bc" << index++;
      constraint.region = parse_region(bc.require("region"),
                                       bc.string_or("name", default_name.str()), dim);
      const std::vector<ConfigNode> fix = bc.array("fix");
      if (fix.empty()) {
        throw ConfigError("'" + bc.path() +
                          ".fix' must list the components to constrain, e.g. "
                          "[\"x\", \"y\"]");
      }
      for (const ConfigNode& component : fix) {
        const std::string c = component.string();
        if (c == "x") {
          constraint.fix_x = true;
        } else if (c == "y") {
          constraint.fix_y = true;
        } else if (c == "z" && dim == 3) {
          constraint.fix_z = true;
        } else if (c == "z") {
          throw ConfigError("'" + component.path() +
                            "' fixes \"z\" but the mesh is two-dimensional");
        } else {
          throw ConfigError("'" + component.path() + "' must be \"x\", \"y\"" +
                            (dim == 3 ? " or \"z\"" : "") + ", got \"" + c + "\"");
        }
      }
      const Vector3 values = bc.vector3_or("value", Vector3::Zero(), dim);
      constraint.value_x = values.x();
      constraint.value_y = values.y();
      constraint.value_z = values.z();
      config.constraints.push_back(std::move(constraint));
    }
  }

  // --- load cases ---------------------------------------------------------
  {
    const std::vector<ConfigNode> cases = root.array("load_cases");
    if (cases.empty()) {
      throw ConfigError("'load_cases' is missing or empty; define at least one");
    }
    int index = 0;
    for (const ConfigNode& lc : cases) {
      LoadCaseSpec spec;
      std::ostringstream default_name;
      default_name << "case" << index++;
      spec.name = lc.string_or("name", default_name.str());
      spec.weight = lc.number_or("weight", 1.0);
      if (!(spec.weight >= 0.0)) {
        std::ostringstream os;
        os << "'" << lc.path() << ".weight' must be non-negative, got " << spec.weight;
        throw ConfigError(os.str());
      }

      int load_index = 0;
      for (const ConfigNode& pl : lc.array("point_loads")) {
        PointLoadSpec load;
        std::ostringstream ln;
        ln << spec.name << "_point" << load_index++;
        load.region =
            parse_region(pl.require("region"), pl.string_or("name", ln.str()), dim);
        load.force = pl.require("force").vector3(dim);
        const std::string distribution = pl.string_or("distribution", "total");
        if (distribution == "total") {
          load.distribute_total = true;
        } else if (distribution == "per_node") {
          load.distribute_total = false;
        } else {
          throw ConfigError("'" + pl.path() + ".distribution' must be 'total' or "
                            "'per_node', got '" + distribution + "'");
        }
        spec.point_loads.push_back(std::move(load));
      }

      int traction_index = 0;
      for (const ConfigNode& tr : lc.array("tractions")) {
        TractionLoadSpec load;
        std::ostringstream ln;
        ln << spec.name << "_traction" << traction_index++;
        load.region =
            parse_region(tr.require("region"), tr.string_or("name", ln.str()), dim);
        load.traction = tr.require("traction").vector3(dim);
        spec.tractions.push_back(std::move(load));
      }

      spec.prescribed_displacement_only =
          lc.boolean_or("prescribed_displacement_only", false);
      if (spec.point_loads.empty() && spec.tractions.empty() &&
          !spec.prescribed_displacement_only) {
        throw ConfigError(
            "load case '" + spec.name +
            "' defines neither point_loads nor tractions. If it is meant to be "
            "driven by prescribed displacements alone, set "
            "\"prescribed_displacement_only\": true");
      }
      config.load_cases.push_back(std::move(spec));
    }
  }

  // --- solver -------------------------------------------------------------
  {
    const ConfigNode solver = root.child("solver");
    const ConfigNode lin = solver.child("linear");
    config.analysis.linear.type =
        parse_linear_solver_type(lin.string_or("type", "auto"));
    config.analysis.linear.iterative_tolerance =
        lin.number_or("iterative_tolerance", 1.0e-12);
    config.analysis.linear.max_iterations = lin.integer_or("max_iterations", 0);
    config.analysis.linear.residual_tolerance =
        lin.number_or("residual_tolerance", 1.0e-8);
    config.analysis.linear.pivot_tolerance = lin.number_or("pivot_tolerance", 1.0e-14);
    config.analysis.linear.warm_start = lin.boolean_or("warm_start", true);
    {
      const ConfigNode limits = lin.child("auto_direct_limit");
      LinearSolverOptions& l = config.analysis.linear;
      l.auto_direct_limit_2d = limits.integer_or("plane", l.auto_direct_limit_2d);
      l.auto_direct_limit_3d = limits.integer_or("solid", l.auto_direct_limit_3d);
      if (l.auto_direct_limit_2d < 0 || l.auto_direct_limit_3d < 0) {
        throw ConfigError("'solver.linear.auto_direct_limit' entries must be non-negative");
      }
    }
    {
      const ConfigNode amg = lin.child("amg");
      AmgOptions& a = config.analysis.linear.amg;
      a.strength_threshold = amg.number_or("strength_threshold", a.strength_threshold);
      a.max_levels = amg.integer_or("max_levels", a.max_levels);
      a.coarse_size = amg.integer_or("coarse_size", a.coarse_size);
      a.smoother = parse_amg_smoother(amg.string_or("smoother", to_string(a.smoother)));
      a.smoother_degree = amg.integer_or("smoother_degree", a.smoother_degree);
      a.chebyshev_ratio = amg.number_or("chebyshev_ratio", a.chebyshev_ratio);
      a.prolongator_damping = amg.number_or("prolongator_damping", a.prolongator_damping);
      a.lanczos_steps = amg.integer_or("lanczos_steps", a.lanczos_steps);
      a.reuse_aggregates = amg.boolean_or("reuse_aggregates", a.reuse_aggregates);
      a.coarse_pivot_tolerance =
          amg.number_or("coarse_pivot_tolerance", a.coarse_pivot_tolerance);
      AmgPreconditioner check(a);  // validates, with the key in the message
      (void)check;
    }
    config.analysis.equilibrium_tolerance =
        solver.number_or("equilibrium_tolerance", 1.0e-6);
    config.analysis.check_model = solver.boolean_or("check_model", true);
  }

  // --- modal --------------------------------------------------------------
  {
    const ConfigNode modal = root.child("modal");
    config.modal.enabled = modal.boolean_or("enabled", false);
    config.modal.options.num_modes = modal.integer_or("num_modes", 6);
    config.modal.options.mass_type =
        parse_mass_type(modal.string_or("mass_type", "consistent"));
    config.modal.options.max_iterations = modal.integer_or("max_iterations", 200);
    config.modal.options.tolerance = modal.number_or("tolerance", 1.0e-10);
    config.modal.options.shift_factor = modal.number_or("shift_factor", 0.0);
    config.modal.options.rigid_body_ratio = modal.number_or("rigid_body_ratio", 1.0e-10);
    config.modal.options.residual_tolerance =
        modal.number_or("residual_tolerance", 1.0e-6);
    config.modal.options.seed =
        static_cast<unsigned int>(modal.integer_or("seed", 20240917));
    config.modal.analyse_optimised_topology =
        modal.boolean_or("analyse_optimised_topology", true);
    config.modal.compare_mass_matched_baseline =
        modal.boolean_or("compare_mass_matched_baseline", true);
    // The eigen solver factorises (or preconditions) with the deck's solver.
    config.modal.options.linear = config.analysis.linear;
    if (config.modal.enabled && config.material().density() <= 0.0) {
      throw ConfigError(
          "modal analysis is enabled but material.density is zero; set a positive "
          "density [kg/m^3]");
    }
  }

  // --- buckling -----------------------------------------------------------
  {
    const ConfigNode buckling = root.child("buckling");
    config.buckling.enabled = buckling.boolean_or("enabled", false);
    config.buckling.options.num_modes = buckling.integer_or("num_modes", 4);
    parse_buckling_solver(buckling, config.buckling.options);
    config.buckling.load_cases = string_list(buckling, "load_cases");
    config.buckling.analyse_optimised_topology =
        buckling.boolean_or("analyse_optimised_topology", true);
    config.buckling.options.linear = config.analysis.linear;
    if (config.buckling.enabled) {
      if (config.buckling.options.num_modes < 1) {
        throw ConfigError("'buckling.num_modes' must be at least 1");
      }
      (void)config.buckling_load_cases();  // validates the names
    }
  }

  // --- topology -----------------------------------------------------------
  {
    const ConfigNode topo = root.child("topology");
    config.topology.enabled = topo.boolean_or("enabled", false);
    config.topology.volume_fraction = topo.number_or("volume_fraction", 0.4);
    config.topology.initial_density = topo.number_or("initial_density", -1.0);

    const ConfigNode simp = topo.child("simp");
    config.topology.optimizer.simp.penalty = simp.number_or("penalty", 3.0);
    config.topology.optimizer.simp.emin_ratio = simp.number_or("emin_ratio", 1.0e-9);
    config.topology.optimizer.simp.mass_floor = simp.number_or("mass_floor", 1.0e-9);
    config.topology.optimizer.simp.mass_law = parse_mass_interpolation(
        simp.string_or("mass_interpolation", "penalty_matched"));

    const ConfigNode filter = topo.child("filter");
    config.topology.filter_type = parse_filter_type(filter.string_or("type", "density"));
    config.topology.filter_radius = filter.number_or("radius", 0.0);
    config.topology.filter_radius_elements = filter.number_or("radius_elements", 1.5);

    const ConfigNode opt = topo.child("optimizer");
    TopologyOptimizerOptions& o = config.topology.optimizer;
    o.max_iterations = opt.integer_or("max_iterations", 200);
    o.change_tolerance = opt.number_or("change_tolerance", 1.0e-2);
    o.objective_tolerance = opt.number_or("objective_tolerance", 0.0);
    o.objective_window = opt.integer_or("objective_window", 5);
    o.continuation_steps = opt.integer_or("continuation_steps", 1);
    o.penalty_start = opt.number_or("penalty_start", 1.0);
    o.continuation_iterations = opt.integer_or("continuation_iterations", 25);
    o.history_stride = opt.integer_or("history_stride", 1);
    o.interpretation_threshold = opt.number_or("interpretation_threshold", 0.5);
    o.oc.move_limit = opt.number_or("move_limit", 0.2);
    o.oc.damping = opt.number_or("damping", 0.5);
    o.oc.volume_tolerance = opt.number_or("volume_tolerance", 1.0e-10);
    o.oc.max_bisections = opt.integer_or("max_bisections", 200);
    o.analysis = config.analysis;

    // Update method: optimality criteria (default, volume constraint only) or
    // MMA (any number of constraints). MMA shares the move limit unless its
    // own block overrides it.
    o.method = parse_optimizer_method(opt.string_or("method", "oc"));
    const ConfigNode mma = opt.child("mma");
    o.mma.move_limit = mma.number_or("move_limit", o.oc.move_limit);
    o.mma.asymptote_init = mma.number_or("asymptote_init", 0.5);
    o.mma.asymptote_increase = mma.number_or("asymptote_increase", 1.2);
    o.mma.asymptote_decrease = mma.number_or("asymptote_decrease", 0.7);
    o.mma.c = mma.number_or("constraint_penalty", 1000.0);
    o.mma.epsimin = mma.number_or("subproblem_tolerance", 1.0e-7);
    o.mma.max_inner_iterations =
        mma.integer_or("max_newton_iterations", o.mma.max_inner_iterations);
    o.mma.validate();
    o.constraint_tolerance = opt.number_or("constraint_tolerance", 1.0e-4);

    // Heaviside projection with beta continuation.
    {
      const ConfigNode proj = topo.child("projection");
      ProjectionOptions& pr = o.projection;
      pr.enabled = proj.boolean_or("enabled", false);
      pr.eta = proj.number_or("eta", pr.eta);
      pr.beta_start = proj.number_or("beta_start", pr.beta_start);
      pr.beta_max = proj.number_or("beta_max", pr.beta_max);
      pr.beta_factor = proj.number_or("beta_factor", pr.beta_factor);
      pr.beta_interval = proj.integer_or("beta_interval", pr.beta_interval);
      pr.advance_on_convergence =
          proj.boolean_or("advance_on_convergence", pr.advance_on_convergence);
      pr.validate();
      if (pr.enabled && config.topology.filter_type == FilterType::Sensitivity) {
        throw ConfigError(
            "'topology.projection' needs the density filter (or none): the sensitivity "
            "filter has no chain rule to extend through the projection");
      }
    }

    // Aggregated von Mises stress constraint (MMA only).
    const ConfigNode stress = topo.child("stress");
    o.stress.enabled = stress.boolean_or("enabled", false);
    o.stress.limit = stress.number_or("limit", 0.0);
    o.stress.p_norm = stress.number_or("p_norm", 8.0);
    o.stress.relaxation = stress.number_or("relaxation", 0.5);
    o.stress.scaling_blend = stress.number_or("scaling_blend", 0.5);
    o.stress.feasibility_tolerance = stress.number_or("feasibility_tolerance", 1.0e-3);
    if (o.stress.enabled && o.method != OptimizerMethod::MMA) {
      throw ConfigError(
          "'topology.stress.enabled' needs 'topology.optimizer.method' = \"mma\"; the "
          "optimality-criteria update cannot handle a second constraint");
    }
    if (o.stress.enabled && config.topology.filter_type == FilterType::Sensitivity) {
      throw ConfigError(
          "'topology.stress' needs the density filter (or none): the sensitivity filter "
          "has no exact chain rule for the stress gradient");
    }
    if (config.topology.enabled) o.stress.validate(o.simp);

    // Aggregated lower bound on the buckling load factors (MMA only).
    const ConfigNode buckle = topo.child("buckling_constraint");
    o.buckling.enabled = buckle.boolean_or("enabled", false);
    o.buckling.min_load_factor = buckle.number_or("min_load_factor", 1.0);
    o.buckling.num_modes = buckle.integer_or("num_modes", 6);
    o.buckling.ks_parameter = buckle.number_or("ks_parameter", 40.0);
    o.buckling.solid_threshold = buckle.number_or("solid_threshold", 0.5);
    o.buckling.load_cases = string_list(buckle, "load_cases");
    parse_buckling_solver(buckle, o.buckling.eigen);
    o.buckling.eigen.linear = config.analysis.linear;
    if (o.buckling.enabled && o.method != OptimizerMethod::MMA) {
      throw ConfigError(
          "'topology.buckling_constraint.enabled' needs 'topology.optimizer.method' = "
          "\"mma\"; the optimality-criteria update cannot handle a second constraint");
    }
    if (o.buckling.enabled && config.topology.filter_type == FilterType::Sensitivity) {
      throw ConfigError(
          "'topology.buckling_constraint' needs the density filter (or none): the "
          "sensitivity filter has no exact chain rule for the load-factor gradient");
    }
    if (config.topology.enabled) o.buckling.validate();
    for (const std::string& name : o.buckling.load_cases) {
      bool found = false;
      for (const LoadCaseSpec& spec : config.load_cases) found = found || spec.name == name;
      if (!found) {
        throw ConfigError("'topology.buckling_constraint.load_cases' names '" + name +
                          "', which is not a load case of the deck");
      }
    }

    int passive_index = 0;
    for (const ConfigNode& region : topo.array("passive_regions")) {
      PassiveRegionSpec spec;
      std::ostringstream default_name;
      default_name << "passive" << passive_index++;
      spec.region = parse_region(region.require("region"),
                                 region.string_or("name", default_name.str()), dim);
      const std::string kind = region.string_or("type", "solid");
      if (kind == "solid") {
        spec.solid = true;
      } else if (kind == "void") {
        spec.solid = false;
      } else {
        throw ConfigError("'" + region.path() + ".type' must be 'solid' or 'void', got '" +
                          kind + "'");
      }
      spec.void_density = region.number_or("density", 0.0);
      config.topology.passive_regions.push_back(std::move(spec));
    }
  }

  // --- output -------------------------------------------------------------
  {
    const ConfigNode out = root.child("output");
    config.output.write_csv = out.boolean_or("csv", true);
    config.output.write_vtk = out.boolean_or("vtk", true);
    config.output.write_density_history = out.boolean_or("density_history", true);
    config.output.write_mode_shapes = out.boolean_or("mode_shapes", true);
  }

  // --- unknown keys -------------------------------------------------------
  const std::vector<std::string> unused = root.unused_keys();
  if (!unused.empty()) {
    std::ostringstream os;
    os << source << " contains " << unused.size()
       << " key(s) that SparLab does not recognise (a misspelled key silently takes "
          "its default, so this is reported):";
    for (const std::string& key : unused) os << "\n  - " << key;
    if (strict) throw ConfigError(os.str());
    log::warn(os.str());
  }

  return config;
}

Configuration load_configuration(const std::string& path, bool strict) {
  const json::Value document = json::parse_file(path);
  const std::string base = std::filesystem::path(path).parent_path().string();
  return parse_configuration(document, path, strict, base.empty() ? "." : base);
}

Mesh build_mesh(const Configuration& config) {
  switch (config.mesh_kind) {
    case MeshKind::StructuredQuad: return make_structured_quad_mesh(config.mesh_spec);
    case MeshKind::StructuredHex: return make_structured_hex_mesh(config.mesh_spec);
    case MeshKind::StructuredTri: return make_structured_tri_mesh(config.mesh_spec);
    case MeshKind::StructuredTet:
      return config.mesh_order == 2 ? make_structured_tet10_mesh(config.mesh_spec)
                                    : make_structured_tet_mesh(config.mesh_spec);
    case MeshKind::File:
      if (config.file_mesh == nullptr) {
        throw ConfigError("the configuration names a mesh file that has not been read");
      }
      return *config.file_mesh;
  }
  throw ConfigError("unhandled mesh type");
}

FemModel build_model(const Configuration& config) {
  Mesh mesh = build_mesh(config);
  FemModel model(std::move(mesh), config.material(), config.thickness,
                 config.stress_state, config.integration);
  model.constraints() = config.constraints;
  model.load_case_specs() = config.load_cases;
  model.finalize();
  return model;
}

DesignDomain build_design_domain(const Configuration& config, const FemModel& model) {
  return DesignDomain(model, config.topology.volume_fraction,
                      config.topology.initial_density,
                      config.topology.passive_regions);
}

}  // namespace sparlab
