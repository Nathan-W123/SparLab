/// \file sparlab_bench.cpp
/// \brief Runtime scaling benchmark.
///
/// Measures wall-clock cost against problem size for the four phases that
/// dominate a topology-optimisation run:
///   * assembly of the global stiffness matrix,
///   * sparse Cholesky factorisation of \f$K_{ff}\f$,
///   * back-substitution for one right-hand side,
///   * one full objective + gradient evaluation (assembly, factorisation,
///     one solve per load case, and the sensitivity loop).
///
/// Each size is repeated `--repeats` times and the *minimum* is reported, which
/// is the standard estimator for wall-clock benchmarks: it is the sample least
/// polluted by scheduling noise. Observed scaling exponents are fitted by
/// least squares on log(time) vs log(DOFs) over the largest three sizes, where
/// the asymptotic behaviour dominates.
///
/// `--solver` selects the linear solver (sparse Cholesky by default, or
/// multigrid-preconditioned CG, Jacobi CG, or the automatic choice) and
/// `--element` the cell type, so the same table compares solvers and
/// elements. For an iterative solver "factorize" is the multigrid setup and
/// "solve" the CG iterations for one load case from a zero initial guess;
/// the objective is timed with warm starts off (each repeat would otherwise
/// start from the previous answer) but with the hierarchy's aggregates
/// reused, which is the steady state of an optimisation loop.

#include "AppSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Timer.hpp"
#include "sparlab/core/Version.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/LinearSolver.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/io/Json.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"
#include "sparlab/topopt/Sensitivity.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#if defined(SPARLAB_HAVE_OPENMP)
#include <omp.h>
#endif

using namespace sparlab;

namespace {

struct BenchRow {
  Index nx = 0;
  Index ny = 0;
  Index nz = 0;
  Index num_elements = 0;
  Index num_dofs = 0;
  Index nonzeros = 0;
  Scalar assemble = 0.0;
  Scalar factorize = 0.0;
  Scalar solve = 0.0;
  Scalar objective = 0.0;
  Scalar compliance = 0.0;
  int iterations = 0;
  int levels = 0;
  Scalar operator_complexity = 0.0;
  Index solver_nonzeros = 0;
};

Mesh build_bench_mesh(ElementType element, const StructuredMeshSpec& spec) {
  switch (element) {
    case ElementType::Quad4: return make_structured_quad_mesh(spec);
    case ElementType::Tri3: return make_structured_tri_mesh(spec);
    case ElementType::Hex8: return make_structured_hex_mesh(spec);
    case ElementType::Tet4: return make_structured_tet_mesh(spec);
  }
  throw ConfigError("unhandled element type");
}

FemModel build_bench_model(ElementType element, Index nx, Index ny, Index nz) {
  StructuredMeshSpec spec;
  spec.nx = nx;
  spec.ny = ny;
  spec.nz = nz;
  spec.lx = 2.0;
  spec.ly = 1.0;
  spec.lz = 0.5;

  const bool solid = element_dimension(element) == 3;
  IsotropicMaterial material(70.0e9, 0.3, 2700.0, "bench");
  FemModel model(build_bench_mesh(element, spec), material, solid ? 1.0 : 0.01,
                 solid ? StressState::ThreeDimensional : StressState::PlaneStress,
                 IntegrationOptions());

  DisplacementConstraint root;
  root.region.name = "left_edge";
  Selector box;
  box.kind = SelectorKind::Box;
  box.xmax = 0.0;
  root.region.members.push_back(box);
  root.fix_x = true;
  root.fix_y = true;
  root.fix_z = solid;
  model.constraints().push_back(root);

  LoadCaseSpec load;
  load.name = "tip";
  PointLoadSpec tip;
  tip.region.name = "tip_node";
  Selector nearest;
  nearest.kind = SelectorKind::NearestNode;
  nearest.point = Vector3(spec.lx, 0.0, 0.0);
  tip.region.members.push_back(nearest);
  tip.force = Vector3(0.0, -1000.0, 0.0);
  load.point_loads.push_back(tip);
  model.load_case_specs().push_back(load);

  model.finalize();
  return model;
}

/// Least-squares slope of log(y) against log(x).
Scalar log_log_slope(const std::vector<Scalar>& x, const std::vector<Scalar>& y) {
  if (x.size() < 2 || x.size() != y.size()) return 0.0;
  Scalar sx = 0.0;
  Scalar sy = 0.0;
  Scalar sxx = 0.0;
  Scalar sxy = 0.0;
  std::size_t n = 0;
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (!(x[i] > 0.0) || !(y[i] > 0.0)) continue;
    const Scalar lx = std::log(x[i]);
    const Scalar ly = std::log(y[i]);
    sx += lx;
    sy += ly;
    sxx += lx * lx;
    sxy += lx * ly;
    ++n;
  }
  if (n < 2) return 0.0;
  const Scalar denom = static_cast<Scalar>(n) * sxx - sx * sx;
  if (std::abs(denom) < 1.0e-300) return 0.0;
  return (static_cast<Scalar>(n) * sxy - sx * sy) / denom;
}

}  // namespace

int main(int argc, char** argv) {
  return app::run_guarded([&]() -> int {
    const std::vector<std::string> known = {"output", "sizes", "repeats", "aspect",
                                            "dim", "solver", "element", "tolerance",
                                            "smoother", "no-objective", "verbosity",
                                            "help"};
    app::CommandLine cli(argc, argv, known);
    if (cli.has("help")) {
      return app::print_usage(
          "sparlab_bench", "[--sizes 20,40,80,160] [--repeats 3] [--output <dir>]",
          {{"--sizes <list>", "comma-separated nx values (ny = nx / aspect)"},
           {"--repeats <n>", "repetitions per size; the minimum is reported"},
           {"--aspect <a>", "nx / ny ratio of the benchmark plate (default 2)"},
           {"--dim <2|3>", "2 for a plate (default), 3 for a block with nz = ny / 2"},
           {"--element <type>", "quad|tri (2-D) or hex|tet (3-D); default quad / hex"},
           {"--solver <type>", "simplicial_ldlt (default), amg_cg, conjugate_gradient "
                               "or auto"},
           {"--tolerance <t>", "relative residual of the iterative solvers (1e-10)"},
           {"--smoother <s>", "multigrid smoother: chebyshev (default) or gauss_seidel"},
           {"--no-objective", "skip the objective + gradient timing"},
           {"--output <dir>", "output directory (default results/benchmark)"},
           {"--verbosity <lvl>", "trace|debug|info|warn|error|silent"},
           {"--help", "show this message"}});
    }
    app::apply_verbosity(cli);

    const std::string out_dir = cli.value("output", "results/benchmark");
    ensure_directory(out_dir);
    const int repeats = std::max(1, cli.integer("repeats", 3));
    const Scalar aspect = cli.number("aspect", 2.0);
    if (!(aspect > 0.0)) throw ConfigError("--aspect must be positive");
    const int dim = cli.integer("dim", 2);
    if (dim != 2 && dim != 3) throw ConfigError("--dim must be 2 or 3");
    const std::string element_name = cli.value("element", dim == 2 ? "quad" : "hex");
    ElementType element = ElementType::Quad4;
    if (element_name == "quad") {
      element = ElementType::Quad4;
    } else if (element_name == "tri") {
      element = ElementType::Tri3;
    } else if (element_name == "hex") {
      element = ElementType::Hex8;
    } else if (element_name == "tet") {
      element = ElementType::Tet4;
    } else {
      throw ConfigError("--element must be quad, tri, hex or tet");
    }
    if (element_dimension(element) != dim) {
      throw ConfigError("--element " + element_name + " does not match --dim " +
                        std::to_string(dim));
    }
    LinearSolverOptions linear;
    linear.type = parse_linear_solver_type(cli.value("solver", "simplicial_ldlt"));
    linear.iterative_tolerance = cli.number("tolerance", 1.0e-10);
    if (cli.has("smoother")) linear.amg.smoother = parse_amg_smoother(cli.value("smoother"));
    const bool run_objective = !cli.has("no-objective");

    std::vector<Scalar> sizes =
        cli.has("sizes") ? cli.number_list("sizes")
        : dim == 2       ? std::vector<Scalar>{20, 40, 80, 160, 240, 320}
                         : std::vector<Scalar>{8, 16, 24, 32, 40};

    // The default configuration keeps the historical file names.
    std::string stem = dim == 2 ? "runtime_scaling" : "runtime_scaling_3d";
    const bool default_element = element == (dim == 2 ? ElementType::Quad4 : ElementType::Hex8);
    if (!default_element) stem += "_" + element_name;
    if (linear.type != LinearSolverType::SimplicialLdlt) stem += "_" + to_string(linear.type);
    CsvWriter csv(path_join(out_dir, stem + ".csv"),
                  {"nx", "ny", "nz", "num_elements", "num_dofs", "stiffness_nonzeros",
                   "assemble[s]", "factorize[s]", "solve[s]", "objective_gradient[s]",
                   "compliance[J]", "iterations", "levels", "operator_complexity",
                   "solver_nonzeros"});

    std::vector<BenchRow> rows;
    std::cout << "element " << to_string(element) << ", solver " << to_string(linear.type)
              << "\n";
    std::cout << std::left << std::setw(8) << "nx" << std::setw(8) << "ny"
              << std::setw(8) << "nz" << std::setw(12) << "elements" << std::setw(12)
              << "DOFs"
              << std::setw(14) << "assemble[s]" << std::setw(14) << "factorize[s]"
              << std::setw(14) << "solve[s]" << std::setw(16) << "obj+grad[s]"
              << std::setw(8) << "iters" << "\n";
    std::cout << std::string(106, '-') << "\n";

    for (Scalar size : sizes) {
      const Index nx = static_cast<Index>(std::llround(size));
      const Index ny = std::max<Index>(1, static_cast<Index>(std::llround(size / aspect)));
      const Index nz = dim == 3 ? std::max<Index>(1, ny / 2) : 0;
      if (nx < 1) throw ConfigError("--sizes must contain positive element counts");

      BenchRow row;
      row.nx = nx;
      row.ny = ny;
      row.nz = nz;

      FemModel model = build_bench_model(element, nx, ny, nz);
      Assembler assembler(model);
      row.num_elements = model.mesh().num_elements();
      row.num_dofs = model.dofs().num_dofs();

      row.assemble = std::numeric_limits<Scalar>::max();
      row.factorize = std::numeric_limits<Scalar>::max();
      row.solve = std::numeric_limits<Scalar>::max();
      row.objective = std::numeric_limits<Scalar>::max();

      for (int r = 0; r < repeats; ++r) {
        Timer timer;
        const SparseMatrix k = assembler.assemble_stiffness();
        row.assemble = std::min(row.assemble, timer.elapsed_seconds());
        const SparseMatrix kff = assembler.reduce_free_free(k);
        row.nonzeros = static_cast<Index>(kff.nonZeros());

        auto solver = make_linear_solver(linear);
        DofLayout layout;
        layout.dim = model.dim();
        layout.coordinates = &model.mesh().coordinates();
        layout.unknowns = &model.dofs().free_dofs();
        solver->set_layout(layout);
        timer.reset();
        solver->factorize(kff);
        row.factorize = std::min(row.factorize, timer.elapsed_seconds());

        const Vector rhs =
            model.dofs().restrict_to_free(model.load_vectors().front());
        timer.reset();
        const Vector uf = solver->solve(rhs);
        row.solve = std::min(row.solve, timer.elapsed_seconds());
        row.compliance = rhs.dot(uf);
        row.iterations = solver->last_iterations();
        row.solver_nonzeros = solver->storage_nonzeros();
        if (const AmgStats* stats = solver->amg_stats()) {
          row.levels = static_cast<int>(stats->levels.size());
          row.operator_complexity = stats->operator_complexity;
        }
      }

      // One full objective + gradient evaluation, which is the unit of cost of
      // a topology-optimisation iteration.
      if (run_objective) {
        const Scalar cell = 2.0 / static_cast<Scalar>(nx);
        DensityFilter filter(model.mesh(), FilterType::Density, 1.5 * cell);
        DesignDomain domain(model, 0.5, 0.5, {});
        SimpOptions simp;
        StaticAnalysisOptions options;
        options.linear = linear;
        options.linear.warm_start = false;
        ComplianceObjective objective(model, assembler, filter, domain, simp, options);
        const Vector x = domain.initial_design();
        for (int r = 0; r < repeats; ++r) {
          Timer timer;
          const ObjectiveEvaluation eval = objective.evaluate(x, true);
          row.objective = std::min(row.objective, timer.elapsed_seconds());
          (void)eval;
        }
      }

      if (!run_objective) row.objective = 0.0;
      csv.row({static_cast<Scalar>(row.nx), static_cast<Scalar>(row.ny),
               static_cast<Scalar>(row.nz), static_cast<Scalar>(row.num_elements),
               static_cast<Scalar>(row.num_dofs), static_cast<Scalar>(row.nonzeros),
               row.assemble, row.factorize, row.solve, row.objective,
               row.compliance, static_cast<Scalar>(row.iterations),
               static_cast<Scalar>(row.levels), row.operator_complexity,
               static_cast<Scalar>(row.solver_nonzeros)});
      rows.push_back(row);

      std::cout << std::left << std::setw(8) << row.nx << std::setw(8) << row.ny
                << std::setw(8) << row.nz << std::setw(12) << row.num_elements << std::setw(12) << row.num_dofs
                << std::setw(14) << app::format(row.assemble, 4) << std::setw(14)
                << app::format(row.factorize, 4) << std::setw(14)
                << app::format(row.solve, 4) << std::setw(16)
                << app::format(row.objective, 4) << std::setw(8) << row.iterations << "\n";
    }
    csv.close();

    // Scaling exponents over the largest three sizes.
    const std::size_t tail = std::min<std::size_t>(3, rows.size());
    std::vector<Scalar> dofs;
    std::vector<Scalar> t_assemble;
    std::vector<Scalar> t_factorize;
    std::vector<Scalar> t_solve;
    std::vector<Scalar> t_objective;
    for (std::size_t i = rows.size() - tail; i < rows.size(); ++i) {
      dofs.push_back(static_cast<Scalar>(rows[i].num_dofs));
      t_assemble.push_back(rows[i].assemble);
      t_factorize.push_back(rows[i].factorize);
      t_solve.push_back(rows[i].solve);
      t_objective.push_back(rows[i].objective);
    }

    json::Value summary = json::Value::make_object();
    summary.set("code", json::Value::make_string("SparLab"));
    summary.set("version", json::Value::make_string(SPARLAB_VERSION_STRING));
    summary.set("build_type", json::Value::make_string(SPARLAB_BUILD_TYPE));
    summary.set("compiler", json::Value::make_string(SPARLAB_COMPILER));
    summary.set("repeats", json::Value::make_number(repeats));
    summary.set("estimator",
                json::Value::make_string("minimum over repeats (least noise-polluted)"));
    summary.set("linear_solver", json::Value::make_string(to_string(linear.type)));
    if (linear.type != LinearSolverType::SimplicialLdlt) {
      summary.set("iterative_tolerance", json::Value::make_number(linear.iterative_tolerance));
      json::Value amg = json::Value::make_object();
      amg.set("smoother", json::Value::make_string(to_string(linear.amg.smoother)));
      amg.set("smoother_degree", json::Value::make_number(linear.amg.smoother_degree));
      amg.set("strength_threshold", json::Value::make_number(linear.amg.strength_threshold));
      amg.set("coarse_size", json::Value::make_number(linear.amg.coarse_size));
      summary.set("multigrid", amg);
    }
#if defined(SPARLAB_HAVE_OPENMP)
    summary.set("openmp_threads", json::Value::make_number(omp_get_max_threads()));
#else
    summary.set("openmp_threads", json::Value::make_number(1));
#endif
    summary.set("dim", json::Value::make_number(dim));
    summary.set("element", json::Value::make_string(to_string(element)));

    json::Value exponents = json::Value::make_object();
    exponents.set("assemble", json::Value::make_number(log_log_slope(dofs, t_assemble)));
    exponents.set("factorize",
                  json::Value::make_number(log_log_slope(dofs, t_factorize)));
    exponents.set("solve", json::Value::make_number(log_log_slope(dofs, t_solve)));
    exponents.set("objective_gradient",
                  json::Value::make_number(log_log_slope(dofs, t_objective)));
    summary.set("scaling_exponent_vs_dofs", exponents);
    summary.set("note",
                json::Value::make_string(
                    dim == 2
                        ? "Exponents are least-squares slopes of log(time) vs log(DOFs) "
                          "over the largest three sizes. Assembly is O(n). A 2-D sparse "
                          "Cholesky with a good fill-reducing ordering is close to "
                          "O(n^1.5) in theory; on these sizes the measured value also "
                          "carries cache effects."
                        : "Exponents are least-squares slopes of log(time) vs log(DOFs) "
                          "over the largest three sizes. Assembly is O(n). A 3-D sparse "
                          "Cholesky with nested-dissection-quality ordering is O(n^2) "
                          "in theory and its fill grows as O(n^(4/3)); the AMD ordering "
                          "used here is somewhat worse, which is the cost that makes "
                          "3-D optimisation loops expensive."));

    json::Value records = json::Value::make_array();
    for (const BenchRow& row : rows) {
      json::Value rec = json::Value::make_object();
      rec.set("nx", json::Value::make_number(row.nx));
      rec.set("ny", json::Value::make_number(row.ny));
      rec.set("nz", json::Value::make_number(row.nz));
      rec.set("num_elements", json::Value::make_number(row.num_elements));
      rec.set("num_dofs", json::Value::make_number(row.num_dofs));
      rec.set("stiffness_nonzeros", json::Value::make_number(row.nonzeros));
      rec.set("assemble_s", json::Value::make_number(row.assemble));
      rec.set("factorize_s", json::Value::make_number(row.factorize));
      rec.set("solve_s", json::Value::make_number(row.solve));
      rec.set("objective_gradient_s", json::Value::make_number(row.objective));
      rec.set("iterations", json::Value::make_number(row.iterations));
      rec.set("levels", json::Value::make_number(row.levels));
      rec.set("operator_complexity", json::Value::make_number(row.operator_complexity));
      rec.set("solver_nonzeros", json::Value::make_number(row.solver_nonzeros));
      records.push_back(rec);
    }
    summary.set("records", records);

    std::ofstream out(path_join(out_dir, stem + ".json"));
    if (!out) throw IoError("cannot write the benchmark summary");
    out << json::dump(summary, 2) << '\n';
    out.close();

    std::cout << std::string(106, '-') << "\n";
    std::cout << "scaling exponents vs DOFs (largest three sizes): assemble "
              << app::format(log_log_slope(dofs, t_assemble), 3) << ", factorize "
              << app::format(log_log_slope(dofs, t_factorize), 3) << ", solve "
              << app::format(log_log_slope(dofs, t_solve), 3) << ", obj+grad "
              << app::format(log_log_slope(dofs, t_objective), 3) << "\n";
    std::cout << "results: " << out_dir << "\n";
    return 0;
  });
}
