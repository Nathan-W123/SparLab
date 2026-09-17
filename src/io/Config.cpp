#include "sparlab/io/Config.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"

#include <cmath>
#include <sstream>

namespace sparlab {
namespace {

StressState parse_stress_state(const std::string& text) {
  if (text == "plane_stress") return StressState::PlaneStress;
  if (text == "plane_strain") return StressState::PlaneStrain;
  throw ConfigError("unknown stress state '" + text +
                    "' (expected plane_stress|plane_strain)");
}

MassType parse_mass_type(const std::string& text) {
  if (text == "consistent") return MassType::Consistent;
  if (text == "lumped") return MassType::Lumped;
  throw ConfigError("unknown mass type '" + text + "' (expected consistent|lumped)");
}

Selector parse_selector_primitive(const ConfigNode& node) {
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
    if (sel.xmin > sel.xmax || sel.ymin > sel.ymax) {
      std::ostringstream os;
      os << "'" << box.path() << "' has an empty interval: x in [" << sel.xmin << ", "
         << sel.xmax << "], y in [" << sel.ymin << ", " << sel.ymax << "]";
      throw ConfigError(os.str());
    }
    ++matched;
  }
  if (node.child("circle").exists()) {
    const ConfigNode c = node.child("circle");
    sel.kind = SelectorKind::Circle;
    sel.center = c.require("center").vector2();
    sel.radius = c.positive_number("radius");
    ++matched;
  }
  if (node.child("annulus").exists()) {
    const ConfigNode c = node.child("annulus");
    sel.kind = SelectorKind::Annulus;
    sel.center = c.require("center").vector2();
    sel.inner_radius = c.number_or("inner_radius", 0.0);
    sel.radius = c.positive_number("radius");
    if (sel.inner_radius >= sel.radius) {
      std::ostringstream os;
      os << "'" << c.path() << "' needs inner_radius (" << sel.inner_radius
         << " m) strictly below radius (" << sel.radius << " m)";
      throw ConfigError(os.str());
    }
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
    sel.point = node.child("nearest_node").vector2();
    ++matched;
  }

  if (matched == 0) {
    throw ConfigError(
        "'" + node.path() +
        "' does not name a region primitive; expected one of all, box, circle, "
        "annulus, node_ids, element_ids, nearest_node (or an 'any_of' list)");
  }
  if (matched > 1) {
    throw ConfigError("'" + node.path() +
                      "' names more than one region primitive; use 'any_of' to combine "
                      "several regions");
  }
  return sel;
}

}  // namespace

SelectorGroup parse_region(const ConfigNode& node, const std::string& default_name) {
  if (!node.exists()) {
    throw ConfigError("required region '" + node.path() + "' is missing");
  }
  SelectorGroup group;
  group.name = node.string_or("name", default_name);
  group.invert = node.boolean_or("invert", false);

  const std::vector<ConfigNode> any_of = node.array("any_of");
  if (!any_of.empty()) {
    for (const ConfigNode& item : any_of) {
      group.members.push_back(parse_selector_primitive(item));
    }
  } else {
    group.members.push_back(parse_selector_primitive(node));
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

Scalar Configuration::resolved_filter_radius(const Mesh& mesh) const {
  if (topology.filter_radius > 0.0) return topology.filter_radius;
  if (!(topology.filter_radius_elements > 0.0)) {
    throw ConfigError(
        "topology.filter needs either 'radius' (metres) or a positive "
        "'radius_elements'");
  }
  Scalar area_sum = 0.0;
  for (Index e = 0; e < mesh.num_elements(); ++e) area_sum += mesh.element_area(e);
  const Scalar mean_size =
      std::sqrt(area_sum / static_cast<Scalar>(mesh.num_elements()));
  return topology.filter_radius_elements * mean_size;
}

Configuration parse_configuration(const json::Value& document, const std::string& source,
                                  bool strict) {
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
    if (type != "structured_quad") {
      throw ConfigError("'" + mesh.path() + ".type' must be 'structured_quad'; '" +
                        type +
                        "' is not implemented (the mesh layer is designed to accept "
                        "further generators)");
    }
    config.mesh_spec.nx = mesh.require("nx").integer();
    config.mesh_spec.ny = mesh.require("ny").integer();
    config.mesh_spec.lx = mesh.positive_number("lx");
    config.mesh_spec.ly = mesh.positive_number("ly");
    config.mesh_spec.x0 = mesh.number_or("x0", 0.0);
    config.mesh_spec.y0 = mesh.number_or("y0", 0.0);
  }

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
    config.thickness = model.number_or("thickness", 1.0);
    if (!(config.thickness > 0.0)) {
      std::ostringstream os;
      os << "'model.thickness' must be positive, got " << config.thickness << " m";
      throw ConfigError(os.str());
    }
    config.stress_state =
        parse_stress_state(model.string_or("stress_state", "plane_stress"));
    const ConfigNode integ = model.child("integration");
    config.integration.stiffness_points = integ.integer_or("stiffness_points", 2);
    config.integration.mass_points = integ.integer_or("mass_points", 3);
    config.integration.edge_points = integ.integer_or("edge_points", 2);
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
                                       bc.string_or("name", default_name.str()));
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
        } else {
          throw ConfigError("'" + component.path() + "' must be \"x\" or \"y\", got \"" +
                            c + "\"");
        }
      }
      const Vector2 values = bc.vector2_or("value", Vector2::Zero());
      constraint.value_x = values.x();
      constraint.value_y = values.y();
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
        load.region = parse_region(pl.require("region"), pl.string_or("name", ln.str()));
        load.force = pl.require("force").vector2();
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
        load.region = parse_region(tr.require("region"), tr.string_or("name", ln.str()));
        load.traction = tr.require("traction").vector2();
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
        parse_linear_solver_type(lin.string_or("type", "simplicial_ldlt"));
    config.analysis.linear.iterative_tolerance =
        lin.number_or("iterative_tolerance", 1.0e-12);
    config.analysis.linear.max_iterations = lin.integer_or("max_iterations", 0);
    config.analysis.linear.residual_tolerance =
        lin.number_or("residual_tolerance", 1.0e-8);
    config.analysis.linear.pivot_tolerance = lin.number_or("pivot_tolerance", 1.0e-14);
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
    if (config.modal.enabled && config.material().density() <= 0.0) {
      throw ConfigError(
          "modal analysis is enabled but material.density is zero; set a positive "
          "density [kg/m^3]");
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

    int passive_index = 0;
    for (const ConfigNode& region : topo.array("passive_regions")) {
      PassiveRegionSpec spec;
      std::ostringstream default_name;
      default_name << "passive" << passive_index++;
      spec.region =
          parse_region(region.require("region"),
                       region.string_or("name", default_name.str()));
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
  return parse_configuration(document, path, strict);
}

FemModel build_model(const Configuration& config) {
  Mesh mesh = make_structured_quad_mesh(config.mesh_spec);
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
