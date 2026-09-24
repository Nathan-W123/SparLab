/// \file sparlab_topopt.cpp
/// \brief SIMP minimum-compliance topology optimisation of one configuration.
///
/// Pipeline:
///   1. build the FE model and check that the *solid* model is well posed;
///   2. build the density filter and the design domain (passive regions);
///   3. run the optimality-criteria loop;
///   4. recover stresses on the final design;
///   5. interpret the density field as a solid body (threshold + largest
///      connected group) and, where the boundary conditions still apply,
///      analyse that body directly;
///   6. compare natural frequencies of the full solid domain, an equal-mass
///      uniform plate, and the optimised topology;
///   7. write every artefact plus a summary.
///
/// Step 6 uses a property specific to the 2-D idealisation: uniformly scaling
/// the thickness scales K and M equally, so an equal-mass *uniform* plate has
/// exactly the same natural frequencies as the full solid domain while its
/// compliance rises by 1/nu. The run verifies that numerically rather than
/// asserting it, and the equal-mass compliance is the baseline the optimised
/// design is judged against. A 3-D solid has no thickness to scale, so a 3-D
/// run reports the full solid as its only reference and says so.

#include "AppSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Timer.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/ModalAnalysis.hpp"
#include "sparlab/fem/ModelDiagnostics.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/io/Config.hpp"
#include "sparlab/io/ResultWriter.hpp"
#include "sparlab/mesh/SubMesh.hpp"
#include "sparlab/topopt/TopologyOptimizer.hpp"

#include <iostream>
#include <memory>

using namespace sparlab;

namespace {

/// Build a solid FE model on an extracted sub-mesh, re-applying the geometric
/// boundary conditions and (optionally) the load cases.
struct SubModelBuild {
  std::unique_ptr<FemModel> model;
  bool loads_applied = false;
  std::string note;
};

SubModelBuild build_solid_submodel(const Configuration& config, Mesh sub_mesh,
                                   bool with_loads) {
  SubModelBuild out;
  auto model = std::make_unique<FemModel>(std::move(sub_mesh), config.material(),
                                          config.thickness, config.stress_state,
                                          config.integration);
  model->constraints() = config.constraints;
  if (with_loads) {
    model->load_case_specs() = config.load_cases;
    try {
      model->finalize(/*require_load_cases=*/true);
      out.loads_applied = true;
    } catch (const ConfigError& e) {
      // The load region no longer intersects the retained material. Report it
      // instead of quietly skipping the comparison.
      out.note = std::string("load cases could not be applied to the interpreted "
                             "solid structure: ") +
                 e.what();
      model = std::make_unique<FemModel>(Mesh(model->mesh()), config.material(),
                                         config.thickness, config.stress_state,
                                         config.integration);
      model->constraints() = config.constraints;
      model->finalize(/*require_load_cases=*/false);
    }
  } else {
    model->finalize(/*require_load_cases=*/false);
  }
  out.model = std::move(model);
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  return app::run_guarded([&]() -> int {
    const std::vector<std::string> known = {
        "config",        "output",      "verbosity",   "strict-config",
        "volume-fraction", "penalty",   "filter-radius", "filter-radius-elements",
        "filter-type",   "max-iterations", "nx",       "ny",           "nz",
        "youngs-modulus", "load-weights", "modes",     "no-vtk",
        "no-csv",        "tag",         "help"};
    app::CommandLine cli(argc, argv, known);
    if (cli.has("help") || argc == 1) {
      return app::print_usage(
          "sparlab_topopt", "--config <deck.json> [--output <dir>] [overrides]",
          {{"--config <file>", "input deck with a topology section"},
           {"--output <dir>", "output directory (default results/<case>)"},
           {"--volume-fraction <v>", "override topology.volume_fraction"},
           {"--penalty <p>", "override topology.simp.penalty"},
           {"--filter-radius <r>", "override the filter radius in metres"},
           {"--filter-radius-elements <k>", "filter radius as a multiple of the cell size"},
           {"--filter-type <t>", "none|density|sensitivity"},
           {"--max-iterations <n>", "override optimizer.max_iterations"},
           {"--nx <n> --ny <n> [--nz <n>]", "override the mesh resolution"},
           {"--youngs-modulus <E>", "override material.youngs_modulus [Pa]"},
           {"--load-weights <w1,w2,...>", "override the per-load-case weights"},
           {"--modes <n>", "enable modal analysis with n modes"},
           {"--tag <name>", "suffix appended to the case name in summary.json"},
           {"--no-vtk / --no-csv", "skip the corresponding output"},
           {"--strict-config", "treat unknown configuration keys as errors"},
           {"--verbosity <lvl>", "trace|debug|info|warn|error|silent"},
           {"--help", "show this message"}});
    }
    app::apply_verbosity(cli);

    Configuration config =
        load_configuration(cli.require("config"), cli.has("strict-config"));

    // ---- command-line overrides (used by the parametric study) -------------
    if (cli.has("volume-fraction")) {
      config.topology.volume_fraction = cli.number("volume-fraction", 0.4);
    }
    if (cli.has("penalty")) {
      config.topology.optimizer.simp.penalty = cli.number("penalty", 3.0);
    }
    if (cli.has("filter-radius")) {
      config.topology.filter_radius = cli.number("filter-radius", 0.0);
      config.topology.filter_radius_elements = 0.0;
    }
    if (cli.has("filter-radius-elements")) {
      config.topology.filter_radius_elements =
          cli.number("filter-radius-elements", 1.5);
      config.topology.filter_radius = 0.0;
    }
    if (cli.has("filter-type")) {
      config.topology.filter_type = parse_filter_type(cli.value("filter-type"));
    }
    if (cli.has("max-iterations")) {
      config.topology.optimizer.max_iterations = cli.integer("max-iterations", 200);
    }
    if (cli.has("nx")) config.mesh_spec.nx = cli.integer("nx", config.mesh_spec.nx);
    if (cli.has("ny")) config.mesh_spec.ny = cli.integer("ny", config.mesh_spec.ny);
    if (cli.has("nz")) {
      if (config.dim() != 3) {
        throw ConfigError("--nz applies to a structured_hex mesh only");
      }
      config.mesh_spec.nz = cli.integer("nz", config.mesh_spec.nz);
    }
    if (cli.has("youngs-modulus")) {
      config.set_material(config.material().with_youngs_modulus(
          cli.number("youngs-modulus", config.material().youngs_modulus())));
    }
    if (cli.has("load-weights")) {
      const std::vector<Scalar> weights = cli.number_list("load-weights");
      if (weights.size() != config.load_cases.size()) {
        std::ostringstream os;
        os << "--load-weights supplied " << weights.size() << " weight(s) but the deck "
           << "defines " << config.load_cases.size() << " load case(s)";
        throw ConfigError(os.str());
      }
      for (std::size_t l = 0; l < weights.size(); ++l) {
        config.load_cases[l].weight = weights[l];
      }
    }
    if (cli.has("modes")) {
      config.modal.enabled = true;
      config.modal.options.num_modes = cli.integer("modes", 6);
    }
    if (cli.has("no-vtk")) config.output.write_vtk = false;
    if (cli.has("no-csv")) config.output.write_csv = false;
    if (cli.has("tag")) config.name += "_" + cli.value("tag");

    if (!config.topology.enabled) {
      throw ConfigError(
          "the deck has topology.enabled = false; sparlab_topopt needs a topology "
          "section. Use sparlab_solve for a plain analysis");
    }

    const std::string out_dir =
        cli.value("output", app::default_output_directory(config.name));

    TimingLedger timings;
    Timer wall;

    FemModel model = [&] {
      ScopedTimer t(timings, "model_build");
      return build_model(config);
    }();

    // The solid model must be well posed; the SIMP floor would otherwise mask a
    // missing boundary condition behind a very soft but non-singular matrix.
    require_well_posed(model);

    Assembler assembler(model);
    const Scalar filter_radius = config.resolved_filter_radius(model.mesh());
    DensityFilter filter(model.mesh(), config.topology.filter_type, filter_radius);
    DesignDomain domain = build_design_domain(config, model);

    TopologyOptimizationResult result;
    {
      ScopedTimer t(timings, "optimization");
      TopologyOptimizer optimizer(model, assembler, filter, domain,
                                  config.topology.optimizer);
      result = optimizer.run();
    }

    // ---- reference analyses on the solid domain ---------------------------
    Scalar solid_compliance = 0.0;
    {
      ScopedTimer t(timings, "solid_reference");
      StaticAnalysisOptions opts = config.analysis;
      StaticAnalysis analysis(model, assembler, opts);
      const std::vector<StaticSolution> sols = analysis.solve_all();
      solid_compliance =
          StaticAnalysis::weighted_compliance(sols, model.normalised_weights());
    }

    // ---- stresses of the optimised design ---------------------------------
    std::vector<StressField> stresses;
    {
      ScopedTimer t(timings, "stress_recovery");
      for (const Vector& u : result.displacements) {
        stresses.push_back(
            recover_stresses(model, assembler, u, &result.stiffness_factors));
      }
    }

    // ---- interpret the density field as a solid body -----------------------
    std::unique_ptr<TopologyInterpretation> interpretation;
    {
      ScopedTimer t(timings, "interpretation");
      interpretation = std::make_unique<TopologyInterpretation>(
          interpret_density_as_solid(model.mesh(), result.physical_density,
                                     domain.element_volumes(),
                                     config.topology.optimizer.interpretation_threshold));
    }

    // ---- analysis of the interpreted solid --------------------------------
    // Thresholding the density field is the only way a SIMP result becomes a
    // structure, so what *that* structure delivers is reported on every run,
    // not only when modal analysis happens to be switched on. It is one extra
    // factorisation of a smaller model.
    //
    // The interpreted structure can legitimately be ill-posed - the largest
    // connected group need not reach a support, and on a coarse mesh with an
    // unresolved filter it often does not. That is a fact about the design,
    // not a reason to throw away a completed optimisation, so the failure is
    // recorded in the summary and warned about rather than made fatal.
    std::unique_ptr<SubModelBuild> sub_build;
    std::unique_ptr<Assembler> sub_assembler;
    json::Value interpreted = json::Value::make_object();
    {
      ScopedTimer t(timings, "interpreted_analysis");
      try {
        sub_build = std::make_unique<SubModelBuild>(build_solid_submodel(
            config, interpretation->sub.mesh, /*with_loads=*/true));
        require_well_posed(*sub_build->model);
        sub_assembler = std::make_unique<Assembler>(*sub_build->model);

        const SubModelBuild& sub = *sub_build;
        interpreted.set("loads_applied", json::Value::make_bool(sub.loads_applied));
        if (!sub.note.empty()) {
          interpreted.set("note", json::Value::make_string(sub.note));
        }
        interpreted.set("num_elements",
                        json::Value::make_number(sub.model->mesh().num_elements()));
        interpreted.set("num_nodes",
                        json::Value::make_number(sub.model->mesh().num_nodes()));
        interpreted.set(
            "mass_kg",
            json::Value::make_number(sub.model->domain_volume() *
                                     config.material().density()));
        if (sub.loads_applied) {
          StaticAnalysisOptions opts = config.analysis;
          StaticAnalysis sub_analysis(*sub.model, *sub_assembler, opts);
          const std::vector<StaticSolution> sols = sub_analysis.solve_all();
          const Scalar c = StaticAnalysis::weighted_compliance(
              sols, sub.model->normalised_weights());
          interpreted.set("weighted_compliance_J", json::Value::make_number(c));
          // The gap between the penalised relaxed model and the structure a
          // threshold actually produces. It can fall either side of 1:
          // thresholding promotes every element above the threshold to full
          // material, which the SIMP penalty had been crediting with only
          // rho^p of its stiffness (stiffer, ratio < 1), and deletes
          // everything below it (softer, ratio > 1). Which effect wins is a
          // property of the design, so it is measured rather than assumed.
          interpreted.set("compliance_vs_simp_ratio",
                          json::Value::make_number(
                              result.compliance > 0.0 ? c / result.compliance : 0.0));
          Scalar max_vm = 0.0;
          for (const StaticSolution& s : sols) {
            const StressField f =
                recover_stresses(*sub.model, *sub_assembler, s.displacement);
            max_vm = std::max(max_vm, f.element_von_mises.maxCoeff());
          }
          interpreted.set("max_von_mises_Pa", json::Value::make_number(max_vm));
        }
      } catch (const std::exception& error) {
        sub_build.reset();
        sub_assembler.reset();
        interpreted.set("analysis_failed", json::Value::make_bool(true));
        interpreted.set("reason", json::Value::make_string(error.what()));
        log::warn("the density field thresholded at ",
                  config.topology.optimizer.interpretation_threshold,
                  " does not form an analysable structure: ", error.what(),
                  ". The optimisation result itself is unaffected and is still "
                  "reported; treat the interpreted geometry as invalid");
      }
    }

    // ---- modal comparison --------------------------------------------------
    std::unique_ptr<ModalResult> modal_solid;
    std::unique_ptr<ModalResult> modal_thin;
    std::unique_ptr<ModalResult> modal_topology;

    if (config.modal.enabled) {
      ScopedTimer t(timings, "modal_analysis");
      modal_solid = std::make_unique<ModalResult>(
          solve_modal(model, assembler, config.modal.options));

      if (config.modal.compare_mass_matched_baseline && config.dim() == 2) {
        // Equal-mass uniform plate: thickness scaled by the achieved volume
        // fraction. K and M both scale with thickness, so the frequencies must
        // match the full solid; the run checks that instead of assuming it.
        // The argument is a 2-D one, so a solid mesh skips it (see below).
        Configuration thin = config;
        thin.thickness = config.thickness * result.volume_fraction;
        FemModel thin_model = build_model(thin);
        Assembler thin_assembler(thin_model);
        modal_thin = std::make_unique<ModalResult>(
            solve_modal(thin_model, thin_assembler, config.modal.options));
      }

      if (config.modal.analyse_optimised_topology) {
        if (sub_build != nullptr) {
          modal_topology = std::make_unique<ModalResult>(
              solve_modal(*sub_build->model, *sub_assembler, config.modal.options));
        } else {
          log::warn("skipping the modal analysis of the optimised topology: the "
                    "interpreted structure could not be analysed (see the "
                    "interpreted_solid_analysis block of summary.json)");
        }
      }
    }

    // ---- output ------------------------------------------------------------
    ResultWriter writer(out_dir, config);
    {
      ScopedTimer t(timings, "output");
      writer.write_config();
      writer.write_mesh(model);
      writer.write_history(result);
      if (config.output.write_csv) {
        writer.write_density(model.mesh(), domain, result);
        for (std::size_t l = 0; l < result.displacements.size(); ++l) {
          const std::string& name = config.load_cases[l].name;
          writer.write_displacement(model.mesh(), name, result.displacements[l]);
          writer.write_stress(model.mesh(), name, stresses[l],
                              &result.physical_density);
        }
      }
      if (config.output.write_density_history) {
        writer.write_density_history(result);
      }
      if (config.output.write_vtk) {
        for (std::size_t l = 0; l < result.displacements.size(); ++l) {
          writer.write_static_vtk(model.mesh(), config.load_cases[l].name,
                                  result.displacements[l], stresses[l],
                                  &result.physical_density,
                                  &result.stiffness_factors);
        }
      }
      if (modal_solid) writer.write_modal(model.mesh(), *modal_solid, "solid");
      if (modal_topology) {
        writer.write_modal(interpretation->sub.mesh, *modal_topology, "topology");
      }
    }

    timings.add("total", wall.elapsed_seconds());

    json::Value summary = make_topology_summary(
        config, model, domain, filter, result, interpretation.get(), modal_solid.get(),
        modal_topology.get(), modal_thin.get(), timings);

    // Mass-stiffness comparison at equal mass.
    {
      json::Value comparison = json::Value::make_object();
      comparison.set("full_solid_compliance_J",
                     json::Value::make_number(solid_compliance));
      const Scalar equal_mass = result.volume_fraction > 0.0
                                    ? solid_compliance / result.volume_fraction
                                    : 0.0;
      if (config.dim() == 2) {
        comparison.set("equal_mass_uniform_plate_compliance_J",
                       json::Value::make_number(equal_mass));
      }
      comparison.set("optimised_compliance_J",
                     json::Value::make_number(result.compliance));
      if (config.dim() == 2) {
        comparison.set("stiffness_gain_over_equal_mass_plate",
                       json::Value::make_number(
                           result.compliance > 0.0 ? equal_mass / result.compliance : 0.0));
      }
      comparison.set("compliance_penalty_vs_full_solid",
                     json::Value::make_number(solid_compliance > 0.0
                                                  ? result.compliance / solid_compliance
                                                  : 0.0));
      comparison.set(
          "note",
          json::Value::make_string(
              config.dim() == 2
                  ? "In this 2-D idealisation K and M both scale linearly with "
                    "thickness, so an equal-mass uniform plate has the same natural "
                    "frequencies as the full solid domain and a compliance of C_solid / "
                    "volume_fraction. That plate is the equal-mass baseline for the "
                    "optimised design."
                  : "A 3-D solid has no thickness to scale, so the equal-mass uniform "
                    "plate of the 2-D runs does not exist here. The only reference "
                    "reported is the full solid domain; the interpreted structure "
                    "(interpreted_solid_analysis) is the design's own re-solve."));
      if (modal_thin && modal_solid && modal_solid->frequencies_hz.size() > 0) {
        const Eigen::Index n =
            std::min(modal_thin->frequencies_hz.size(), modal_solid->frequencies_hz.size());
        Scalar max_rel = 0.0;
        for (Eigen::Index i = 0; i < n; ++i) {
          const Scalar ref = std::max(std::abs(modal_solid->frequencies_hz(i)), 1e-30);
          max_rel = std::max(max_rel,
                             std::abs(modal_thin->frequencies_hz(i) -
                                      modal_solid->frequencies_hz(i)) /
                                 ref);
        }
        comparison.set("thickness_scaling_frequency_max_relative_difference",
                       json::Value::make_number(max_rel));
      }
      summary.set("mass_stiffness_comparison", comparison);
    }
    if (!interpreted.members().empty()) {
      summary.set("interpreted_solid_analysis", interpreted);
    }
    writer.write_json("summary.json", summary);

    // ---- console report ----------------------------------------------------
    std::cout << "case: " << config.name << "\n";
    std::cout << "  design:      " << domain.describe() << "\n";
    std::cout << "  filter:      " << to_string(filter.type()) << ", radius "
              << app::format(filter.radius()) << " m ("
              << app::format(filter.average_support()) << " elements support)\n";
    std::cout << "  converged:   " << (result.converged ? "yes" : "NO") << " after "
              << result.iterations << " iterations (" << result.linear_solves
              << " linear solves)\n";
    std::cout << "  compliance:  " << app::format(result.compliance)
              << " J (full solid " << app::format(solid_compliance) << " J";
    if (config.dim() == 2) {
      std::cout << ", equal-mass uniform plate "
                << app::format(result.volume_fraction > 0.0
                                   ? solid_compliance / result.volume_fraction
                                   : 0.0)
                << " J";
    }
    std::cout << ")\n";
    std::cout << "  volume:      " << app::format(result.volume) << " m^3, fraction "
              << app::format(result.volume_fraction) << " (target "
              << app::format(domain.volume_fraction()) << ", relative violation "
              << app::format(result.volume_constraint_violation) << ")\n";
    std::cout << "  grey level:  " << app::format(result.gray_level) << "\n";
    std::cout << "  interpreted: threshold "
              << app::format(interpretation->threshold) << " keeps "
              << interpretation->elements_retained << " of "
              << model.mesh().num_elements() << " elements in "
              << interpretation->components_above_threshold << " group(s), discarding "
              << app::format(interpretation->volume_discarded_as_islands)
              << " m^3 as islands\n";
    if (modal_solid) {
      std::cout << "  f1 solid:    " << app::format(modal_solid->frequencies_hz(0))
                << " Hz (mass " << app::format(modal_solid->total_mass) << " kg)\n";
    }
    if (modal_topology) {
      std::cout << "  f1 topology: " << app::format(modal_topology->frequencies_hz(0))
                << " Hz (mass " << app::format(modal_topology->total_mass) << " kg)\n";
    }
    std::cout << "  runtime:     " << app::format(timings.get("total")) << " s\n";
    std::cout << "  results:     " << out_dir << "\n";
    for (const std::string& w : result.warnings) {
      std::cout << "  warning:     " << w << "\n";
    }
    return 0;
  });
}
