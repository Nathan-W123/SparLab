/// \file sparlab_solve.cpp
/// \brief Linear static (and optional modal) analysis of one configuration.
///
/// Solves every load case of the deck, recovers stresses and reactions, checks
/// global equilibrium, optionally runs a modal analysis, and writes the full
/// result set to the output directory.

#include "AppSupport.hpp"

#include "sparlab/core/Timer.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/ModalAnalysis.hpp"
#include "sparlab/fem/ModelDiagnostics.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/io/CalculixWriter.hpp"
#include "sparlab/io/Config.hpp"
#include "sparlab/io/ResultWriter.hpp"

#include <iostream>

using namespace sparlab;

int main(int argc, char** argv) {
  return app::run_guarded([&]() -> int {
    const std::vector<std::string> known = {"config",   "output",        "verbosity",
                                            "strict-config", "modes",   "no-vtk",
                                            "no-csv",   "export-calculix", "help"};
    app::CommandLine cli(argc, argv, known);
    if (cli.has("help") || argc == 1) {
      return app::print_usage(
          "sparlab_solve", "--config <deck.json> [--output <dir>]",
          {{"--config <file>", "input deck describing mesh, material, BCs, loads"},
           {"--output <dir>", "output directory (default results/<case>)"},
           {"--modes <n>", "override modal.num_modes and enable modal analysis"},
           {"--no-vtk", "skip VTK output"},
           {"--no-csv", "skip per-node/per-element CSV output"},
           {"--export-calculix", "also write one CalculiX .inp per load case into the "
                                 "output directory (cross-validation)"},
           {"--strict-config", "treat unknown configuration keys as errors"},
           {"--verbosity <lvl>", "trace|debug|info|warn|error|silent"},
           {"--help", "show this message"}});
    }
    app::apply_verbosity(cli);

    Configuration config =
        load_configuration(cli.require("config"), cli.has("strict-config"));
    if (cli.has("modes")) {
      config.modal.enabled = true;
      config.modal.options.num_modes = cli.integer("modes", 6);
    }
    if (cli.has("no-vtk")) config.output.write_vtk = false;
    if (cli.has("no-csv")) config.output.write_csv = false;

    const std::string out_dir =
        cli.value("output", app::default_output_directory(config.name));

    TimingLedger timings;
    Timer wall;

    FemModel model = [&] {
      ScopedTimer t(timings, "model_build");
      return build_model(config);
    }();

    const ModelDiagnostics diagnostics = diagnose_model(model);
    if (!diagnostics.well_posed()) {
      // Report every problem, then fail: an ill-posed static model cannot
      // produce a meaningful answer.
      require_well_posed(model);
    }

    Assembler assembler(model);
    log::info("assembler: element matrix cache ",
              (assembler.uses_element_cache() ? "enabled (uniform structured mesh)"
                                              : "disabled (general mesh)"));

    std::vector<StaticSolution> solutions;
    {
      ScopedTimer t(timings, "static_solve");
      StaticAnalysis analysis(model, assembler, config.analysis);
      solutions = analysis.solve_all();
    }

    std::vector<StressField> stresses;
    {
      ScopedTimer t(timings, "stress_recovery");
      stresses.reserve(solutions.size());
      for (const StaticSolution& sol : solutions) {
        stresses.push_back(recover_stresses(model, assembler, sol.displacement));
      }
    }

    std::unique_ptr<ModalResult> modal;
    if (config.modal.enabled) {
      ScopedTimer t(timings, "modal_analysis");
      modal = std::make_unique<ModalResult>(
          solve_modal(model, assembler, config.modal.options));
    }

    ResultWriter writer(out_dir, config);
    {
      ScopedTimer t(timings, "output");
      writer.write_config();
      writer.write_mesh(model);
      for (std::size_t l = 0; l < solutions.size(); ++l) {
        const std::string& name = solutions[l].load_case_name;
        if (config.output.write_csv) {
          writer.write_displacement(model.mesh(), name, solutions[l].displacement);
          writer.write_stress(model.mesh(), name, stresses[l]);
          writer.write_reactions(model.mesh(), model.dofs(), name,
                                 solutions[l].reactions);
        }
        if (config.output.write_vtk) {
          writer.write_static_vtk(model.mesh(), name, solutions[l].displacement,
                                  stresses[l], nullptr, nullptr);
        }
      }
      if (modal) writer.write_modal(model.mesh(), *modal);
      if (cli.has("export-calculix")) {
        const std::vector<std::string> decks =
            write_calculix_decks(model, writer.file("calculix"), config.name);
        for (const std::string& deck : decks) log::info("wrote CalculiX deck ", deck);
      }
    }

    timings.add("total", wall.elapsed_seconds());
    writer.write_json("summary.json",
                      make_static_summary(config, model, diagnostics, solutions,
                                          stresses, modal.get(), timings));

    std::cout << "case: " << config.name << "\n";
    std::cout << "  mesh:        " << model.mesh().num_elements() << " elements, "
              << model.dofs().num_dofs() << " DOFs (" << model.dofs().num_free()
              << " free)\n";
    for (std::size_t l = 0; l < solutions.size(); ++l) {
      const StaticSolution& s = solutions[l];
      std::cout << "  load case '" << s.load_case_name << "': compliance "
                << app::format(s.compliance) << " J, strain energy "
                << app::format(s.strain_energy) << " J, max |u| "
                << app::format(s.max_displacement_magnitude) << " m, max von Mises "
                << app::format(stresses[l].element_von_mises.maxCoeff())
                << " Pa\n";
      const auto vec = [&](const Vector3& v) {
        std::string text = "(" + app::format(v.x()) + ", " + app::format(v.y());
        if (model.dim() == 3) text += ", " + app::format(v.z());
        return text + ")";
      };
      std::cout << "      equilibrium: applied " << vec(s.equilibrium.applied_force)
                << " N, reactions " << vec(s.equilibrium.reaction_force)
                << " N, relative error "
                << app::format(s.equilibrium.relative_force_error) << "\n";
    }
    if (modal) {
      std::cout << "  modal: ";
      for (Eigen::Index i = 0; i < modal->frequencies_hz.size(); ++i) {
        std::cout << (i ? ", " : "") << "f" << i + 1 << " = "
                  << app::format(modal->frequencies_hz(i)) << " Hz";
      }
      std::cout << "\n";
      std::cout << "      total mass " << app::format(modal->total_mass) << " kg\n";
    }
    std::cout << "  runtime:     " << app::format(timings.get("total")) << " s\n";
    std::cout << "  results:     " << out_dir << "\n";
    return 0;
  });
}
