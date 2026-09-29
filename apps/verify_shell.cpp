/// \file verify_shell.cpp
/// \brief Verification of the MITC4 shell against exact solutions of the
///        continuum model it discretises, and its validation against the
///        shell benchmarks of MacNeal and Harder.
///
/// Studies:
///   * `shell-patch`              constant membrane and bending states on
///                                distorted meshes in a turned plane, and
///                                rigid motions of curved panels: exact;
///   * `shell-plate`              a simply supported square plate under a
///                                uniform pressure, t / a from 1e-1 to 1e-4,
///                                on regular and distorted meshes, against the
///                                exact Reissner-Mindlin deflection (no shear
///                                locking: the error does not depend on t / a);
///                                and a thin clamped plate against Kirchhoff;
///   * `shell-plate-modes`        the six lowest frequencies of the simply
///                                supported plate against the exact
///                                Reissner-Mindlin frequencies (shear and
///                                rotary inertia);
///   * `shell-plate-buckling`     the simply supported plate under uniaxial
///                                and equal biaxial compression against the
///                                exact buckling loads of the model (k = 4 and
///                                2 in the thin limit);
///   * `shell-cylinder-pressure`  a slice of a long cylinder under internal
///                                pressure, its end rings held against
///                                rotation: the exact thick-ring state
///                                u_r = p R / (E ln((R + t/2) / (R - t/2))),
///                                within (t / R)^2 / 12 of p R^2 / (E t);
///   * `shell-scordelis-lo`, `shell-pinched-cylinder`,
///     `shell-pinched-hemisphere` the MacNeal-Harder benchmarks, against the
///                                references of thin-shell theory (validation);
///   * `shell-box-beam`           a cantilever of square box section - four
///                                walls meeting at folds - in bending and in
///                                torsion against beam theory and Bredt, and
///                                its sensitivity to the drilling stiffness.
///
/// The exact references of the plate studies are the solutions of the
/// Reissner-Mindlin plate with the shear factor 5/6, the model a flat MITC4
/// mesh discretises; for buckling, with the geometric stiffness of the
/// degenerated solid, which adds the in-plane stress acting on the fibres'
/// rotations, (t^2 / 12) N (psi_a,b)^2, to the classical N w_,a w_,b. Each is
/// computed here from its trigonometric series or its 3 x 3 modal problem.
#include "VerifySupport.hpp"

#include "AppSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/elements/Shell4.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/Buckling.hpp"
#include "sparlab/fem/Dynamics.hpp"
#include "sparlab/fem/ModalAnalysis.hpp"
#include "sparlab/fem/StaticAnalysis.hpp"
#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <sstream>
#include <vector>

namespace sparlab {
namespace verify {
namespace {

constexpr Scalar kPi = 3.14159265358979323846;
constexpr Scalar kInf = std::numeric_limits<Scalar>::infinity();
constexpr Scalar kShear = Shell4Element::kShearFactor;

std::string fmt(Scalar v, int digits = 6) { return app::format(v, digits); }

SelectorGroup region(const std::string& name, std::vector<Selector> members) {
  SelectorGroup g;
  g.name = name;
  g.members = std::move(members);
  return g;
}

Selector box(Scalar xmin, Scalar xmax, Scalar ymin = -kInf, Scalar ymax = kInf,
             Scalar zmin = -kInf, Scalar zmax = kInf) {
  Selector s;
  s.kind = SelectorKind::Box;
  s.xmin = xmin;
  s.xmax = xmax;
  s.ymin = ymin;
  s.ymax = ymax;
  s.zmin = zmin;
  s.zmax = zmax;
  return s;
}

Selector everything() { return Selector(); }

DisplacementConstraint fix(const SelectorGroup& where, std::initializer_list<int> components) {
  DisplacementConstraint c;
  c.region = where;
  for (int k : components) c.set(k, true);
  return c;
}

StaticSolution solve_static(FemModel& model) {
  model.finalize();
  Assembler assembler(model);
  StaticAnalysisOptions options;
  options.linear.type = LinearSolverType::SimplicialLdlt;
  StaticAnalysis analysis(model, assembler, options);
  return analysis.solve_all().front();
}

Index node_at(const Mesh& mesh, const Vector3& x) {
  Index best = -1;
  Scalar distance = kInf;
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    const Scalar d = (mesh.node(n) - x).norm();
    if (d < distance) {
      distance = d;
      best = n;
    }
  }
  return best;
}

/// Observed order between consecutive meshes of n and 2n cells per side.
Scalar order(Scalar coarse, Scalar fine) { return observed_order(2.0, coarse, 1.0, fine); }

// ---------------------------------------------------------------------------
// Exact plate solutions
// ---------------------------------------------------------------------------

/// Navier's series for the centre deflection of a simply supported square
/// plate under a uniform load, w D / (q a^4).
Scalar navier_centre() {
  Scalar sum = 0.0;
  for (int m = 1; m <= 401; m += 2) {
    for (int n = 1; n <= 401; n += 2) {
      const Scalar sign = (((m + n) / 2 - 1) % 2 == 0) ? 1.0 : -1.0;
      const Scalar mn2 = static_cast<Scalar>(m * m + n * n);
      sum += sign / (static_cast<Scalar>(m) * n * mn2 * mn2);
    }
  }
  return 16.0 / std::pow(kPi, 6) * sum;
}

/// The Marcus moment at the centre of the square: the solution of
/// grad^2 phi = -q with phi = 0 on the edges, phi / (q a^2), from its single
/// series (the Prandtl function of a square bar).
Scalar marcus_centre() {
  Scalar sum = 0.0;
  for (int m = 1; m <= 41; m += 2) {
    const Scalar sign = ((m - 1) / 2) % 2 == 0 ? 1.0 : -1.0;
    sum += sign / (std::pow(static_cast<Scalar>(m), 3) * std::cosh(m * kPi / 2.0));
  }
  return 0.125 - 4.0 / std::pow(kPi, 3) * sum;
}

/// The Reissner-Mindlin plate in the simply supported mode
/// w = W sin(ax) sin(by), psi_1 = P1 cos(ax) sin(by), psi_2 = P2 sin(ax) cos(by)
/// (a = m pi / lx, b = n pi / ly): its stiffness, and its mass and the
/// geometric stiffness of in-plane forces (nx, ny) per unit load, as 3 x 3
/// matrices in (W, P1, P2), the common factor lx ly / 4 dropped.
struct TrigMode {
  Matrix3 k;
  Matrix3 m;
  Matrix3 g;
};

TrigMode trig_mode(Scalar alpha, Scalar beta, Scalar e, Scalar nu, Scalar t, Scalar rho,
                   Scalar nx, Scalar ny) {
  const Scalar d = e * t * t * t / (12.0 * (1.0 - nu * nu));
  const Scalar kgt = kShear * e / (2.0 * (1.0 + nu)) * t;
  TrigMode out;
  out.k.setZero();
  // Bending: D [k11^2 + k22^2 + 2 nu k11 k22 + (1 - nu) / 2 (2 k12)^2].
  out.k(1, 1) = d * (alpha * alpha + 0.5 * (1.0 - nu) * beta * beta);
  out.k(2, 2) = d * (beta * beta + 0.5 * (1.0 - nu) * alpha * alpha);
  out.k(1, 2) = out.k(2, 1) = d * (nu + 0.5 * (1.0 - nu)) * alpha * beta;
  // Transverse shear: k G t [(w_x + psi_1)^2 + (w_y + psi_2)^2].
  out.k(0, 0) += kgt * (alpha * alpha + beta * beta);
  out.k(0, 1) += kgt * alpha;
  out.k(1, 0) += kgt * alpha;
  out.k(0, 2) += kgt * beta;
  out.k(2, 0) += kgt * beta;
  out.k(1, 1) += kgt;
  out.k(2, 2) += kgt;
  out.m = rho * t * Vector3(1.0, t * t / 12.0, t * t / 12.0).asDiagonal();
  const Scalar s = nx * alpha * alpha + ny * beta * beta;
  out.g = s * Vector3(1.0, t * t / 12.0, t * t / 12.0).asDiagonal();
  return out;
}

/// The smallest eigenvalue of the 3 x 3 pencil (k, b), b positive definite.
Scalar smallest(const Matrix3& k, const Matrix3& b) {
  Eigen::GeneralizedSelfAdjointEigenSolver<Matrix3> ges(k, b);
  return ges.eigenvalues()(0);
}

/// A plate lx x ly of n x n cells in z = 0 (n even), optionally distorted -
/// every interior node but the centre one moved - of thickness t.
FemModel plate(Index n, Scalar lx, Scalar ly, Scalar t, const IsotropicMaterial& material,
               Scalar perturbation = 0.0) {
  ShellMeshSpec spec;
  spec.n1 = n;
  spec.n2 = n;
  spec.lx = lx;
  spec.ly = ly;
  spec.perturbation = perturbation;
  spec.seed = 7u;
  const Mesh generated = make_structured_shell_mesh(spec);
  Matrix coords = generated.coordinates();
  coords.col((n / 2) * (n + 1) + n / 2) = Vector3(0.5 * lx, 0.5 * ly, 0.0);
  Mesh mesh(std::move(coords), generated.connectivity(), ElementType::Shell4);
  mesh.set_node_normals(generated.node_normals());
  return FemModel(std::move(mesh), material, t, StressState::Shell, IntegrationOptions());
}

/// Hard simple supports on the edges of the plate [0, lx] x [0, ly]: w = 0
/// and the rotation along each edge held; `in_plane` also holds the edges'
/// in-plane translations.
void simple_supports(FemModel& model, Scalar lx, Scalar ly, bool in_plane) {
  const SelectorGroup x_edges = region("x_edges", {box(0.0, 0.0), box(lx, lx)});
  const SelectorGroup y_edges = region("y_edges", {box(-kInf, kInf, 0.0, 0.0),
                                                   box(-kInf, kInf, ly, ly)});
  model.constraints().push_back(fix(x_edges, {2, 3}));
  model.constraints().push_back(fix(y_edges, {2, 4}));
  if (in_plane) {
    model.constraints().push_back(fix(x_edges, {0, 1}));
    model.constraints().push_back(fix(y_edges, {0, 1}));
  }
}

LoadCaseSpec uniform_pressure(Scalar p) {
  LoadCaseSpec lc;
  lc.name = "pressure";
  PressureLoadSpec load;
  load.region = region("all", {everything()});
  load.pressure = p;
  lc.pressures.push_back(load);
  return lc;
}

// ---------------------------------------------------------------------------
// Patch tests
// ---------------------------------------------------------------------------

Matrix3 turned_frame(Scalar angle) {
  return (Eigen::AngleAxisd(angle, Vector3(1.0, 2.0, -0.5).normalized()) *
          Eigen::AngleAxisd(-0.4, Vector3::UnitY()))
      .toRotationMatrix();
}

/// Every DOF of `nodes` prescribed from the nodal field `u`.
void prescribe(FemModel& model, const std::vector<Index>& nodes, const Vector& u) {
  for (Index n : nodes) {
    DisplacementConstraint c;
    c.region.name = "node" + std::to_string(n);
    c.region.members.emplace_back();
    c.region.members.back().kind = SelectorKind::NodeIds;
    c.region.members.back().ids.push_back(n);
    for (int k = 0; k < 6; ++k) c.set(k, true, u(6 * n + k));
    model.constraints().push_back(c);
  }
}

Vector solve_prescribed(FemModel& model) {
  LoadCaseSpec lc;
  lc.name = "prescribed";
  lc.prescribed_displacement_only = true;
  model.load_case_specs() = {lc};
  return solve_static(model).displacement;
}

/// The in-plane tensor with components (11, 22, 12) on (a, b), as 3 x 3.
Matrix3 tensor(const Vector3& c, const Vector3& a, const Vector3& b) {
  return c(0) * a * a.transpose() + c(1) * b * b.transpose() +
         c(2) * (a * b.transpose() + b * a.transpose());
}

struct PatchError {
  Scalar displacement = 0.0;  ///< interior nodes, over the largest exact value
  Scalar resultant = 0.0;     ///< N, M and Q, over the largest exact N (or M)
};

/// A constant membrane and / or bending state on the distorted plate turned
/// by `rot`: u = A x in the plane, w quadratic.
PatchError flat_patch(bool membrane, bool bending, Scalar perturbation, unsigned int seed,
                      Scalar angle) {
  const Scalar a = 1.0;
  const Scalar b = 0.8;
  const Scalar t = 0.02;
  const Scalar e_mod = 70.0e9;
  const Scalar nu = 0.3;
  const Matrix3 rot = turned_frame(angle);
  ShellMeshSpec spec;
  spec.n1 = 5;
  spec.n2 = 4;
  spec.lx = a;
  spec.ly = b;
  spec.perturbation = perturbation;
  spec.seed = seed;
  const Mesh flat = make_structured_shell_mesh(spec);
  Mesh mesh(rot * flat.coordinates(), flat.connectivity(), ElementType::Shell4);
  mesh.set_node_normals(rot * flat.node_normals());
  Eigen::Matrix2d amat = Eigen::Matrix2d::Zero();
  if (membrane) amat << 2.0e-4, -3.0e-4, 1.0e-4, -1.5e-4;
  const Scalar kx = bending ? 2.0e-3 : 0.0;
  const Scalar ky = bending ? -1.0e-3 : 0.0;
  const Scalar kxy = bending ? 1.5e-3 : 0.0;
  const Scalar omega = 0.5 * (amat(1, 0) - amat(0, 1));
  const Index nn = mesh.num_nodes();
  Vector exact(6 * nn);
  std::vector<Index> boundary;
  std::vector<Index> interior;
  for (Index k = 0; k < nn; ++k) {
    const Vector3 local = rot.transpose() * mesh.node(k);
    const Scalar x = local(0);
    const Scalar y = local(1);
    const Eigen::Vector2d u = amat * Eigen::Vector2d(x, y) + Eigen::Vector2d(1.0e-4, -2.0e-4);
    const Scalar w = 0.5 * kx * x * x + 0.5 * ky * y * y + kxy * x * y + 3.0e-4 - 2.0e-4 * x;
    const Scalar wx = kx * x + kxy * y - 2.0e-4;
    const Scalar wy = ky * y + kxy * x;
    exact.segment<3>(6 * k) = rot * Vector3(u(0), u(1), w);
    exact.segment<3>(6 * k + 3) = rot * Vector3(wy, -wx, omega);
    const Index i = k % (spec.n1 + 1);
    const Index j = k / (spec.n1 + 1);
    (i == 0 || i == spec.n1 || j == 0 || j == spec.n2 ? boundary : interior).push_back(k);
  }
  FemModel model(std::move(mesh), IsotropicMaterial(e_mod, nu, 2700.0), t, StressState::Shell,
                 IntegrationOptions());
  prescribe(model, boundary, exact);
  const Vector u = solve_prescribed(model);
  PatchError err;
  const Scalar scale = exact.cwiseAbs().maxCoeff();
  for (Index k : interior) {
    err.displacement = std::max(
        err.displacement,
        (u.segment<6>(6 * k) - exact.segment<6>(6 * k)).cwiseAbs().maxCoeff() / scale);
  }
  const Scalar d_b = e_mod * t * t * t / (12.0 * (1.0 - nu * nu));
  const Scalar c = e_mod * t / (1.0 - nu * nu);
  const Vector3 eps(amat(0, 0), amat(1, 1), 0.5 * (amat(0, 1) + amat(1, 0)));
  const Matrix3 n_exact =
      tensor(Vector3(c * (eps(0) + nu * eps(1)), c * (eps(1) + nu * eps(0)),
                     c * (1.0 - nu) * eps(2)),
             rot.col(0), rot.col(1));
  const Matrix3 m_exact = tensor(
      Vector3(-d_b * (kx + nu * ky), -d_b * (ky + nu * kx), -d_b * (1.0 - nu) * kxy),
      rot.col(0), rot.col(1));
  const Scalar n_scale = std::max(n_exact.cwiseAbs().maxCoeff(),
                                  m_exact.cwiseAbs().maxCoeff() / t);
  Assembler assembler(model);
  const ShellField field = recover_shell_resultants(model, assembler, u);
  for (const ShellResultants& r : field.element) {
    const Matrix3 n_h = tensor(r.membrane, r.e1, r.e2);
    const Matrix3 m_h = tensor(r.moment, r.e1, r.e2);
    err.resultant = std::max({err.resultant, (n_h - n_exact).cwiseAbs().maxCoeff() / n_scale,
                              (m_h - m_exact).cwiseAbs().maxCoeff() / (n_scale * t),
                              r.shear.cwiseAbs().maxCoeff() / n_scale});
  }
  return err;
}

/// A rigid motion prescribed on the boundary of a curved panel: every node
/// follows it and no element is strained.
PatchError curved_rigid(ShellShape shape) {
  ShellMeshSpec spec;
  spec.shape = shape;
  spec.radius = 2.0;
  spec.n1 = 5;
  spec.n2 = 4;
  spec.origin = Vector3(0.3, -0.2, 0.1);
  if (shape == ShellShape::Cylinder) {
    spec.axis = 1;
    spec.length = 1.5;
    spec.angle_start = -20.0;
    spec.angle_end = 55.0;
  } else {
    spec.angle_start = 10.0;
    spec.angle_end = 80.0;
    spec.polar_start = 40.0;
    spec.polar_end = 110.0;
  }
  Mesh mesh = make_structured_shell_mesh(spec);
  const Vector3 omega(0.3e-3, -0.5e-3, 0.8e-3);
  const Vector3 shift(1.0e-3, 2.0e-3, -1.0e-3);
  const Index nn = mesh.num_nodes();
  Vector exact(6 * nn);
  std::vector<Index> boundary;
  std::vector<Index> interior;
  for (Index k = 0; k < nn; ++k) {
    exact.segment<3>(6 * k) = shift + omega.cross(mesh.node(k));
    exact.segment<3>(6 * k + 3) = omega;
    const Index i = k % (spec.n1 + 1);
    const Index j = k / (spec.n1 + 1);
    (i == 0 || i == spec.n1 || j == 0 || j == spec.n2 ? boundary : interior).push_back(k);
  }
  const Scalar t = 0.05;
  FemModel model(std::move(mesh), IsotropicMaterial(70.0e9, 0.3, 2700.0), t, StressState::Shell,
                 IntegrationOptions());
  prescribe(model, boundary, exact);
  const Vector u = solve_prescribed(model);
  PatchError err;
  const Scalar scale = exact.cwiseAbs().maxCoeff();
  for (Index k : interior) {
    err.displacement = std::max(
        err.displacement,
        (u.segment<6>(6 * k) - exact.segment<6>(6 * k)).cwiseAbs().maxCoeff() / scale);
  }
  Assembler assembler(model);
  const ShellField field = recover_shell_resultants(model, assembler, u);
  // Against the membrane force of a strain as large as the rotation.
  const Scalar unit = 70.0e9 * t * omega.norm();
  for (const ShellResultants& r : field.element) {
    err.resultant = std::max({err.resultant, r.membrane.cwiseAbs().maxCoeff() / unit,
                              r.moment.cwiseAbs().maxCoeff() / (unit * t),
                              r.shear.cwiseAbs().maxCoeff() / unit});
  }
  return err;
}

}  // namespace

StudyOutcome study_shell_patch(const std::string& out_dir, json::Value& summary) {
  CsvWriter csv(path_join(out_dir, "shell_patch.csv"),
                {"case", "perturbation", "seed", "displacement_error[-]", "resultant_error[-]"});
  json::Value cases = json::Value::make_array();
  Scalar worst = 0.0;
  const auto record = [&](const std::string& name, Scalar perturbation, unsigned int seed,
                          const PatchError& e) {
    worst = std::max({worst, e.displacement, e.resultant});
    csv.raw_row({name, fmt(perturbation, 3), std::to_string(seed), fmt(e.displacement, 4),
                 fmt(e.resultant, 4)});
    json::Value c = json::Value::make_object();
    c.set("case", json::Value::make_string(name));
    c.set("perturbation", json::Value::make_number(perturbation));
    c.set("seed", json::Value::make_number(seed));
    c.set("displacement_error", json::Value::make_number(e.displacement));
    c.set("resultant_error", json::Value::make_number(e.resultant));
    cases.push_back(c);
  };
  for (const Scalar perturbation : {0.0, 0.1, 0.2}) {
    for (const unsigned int seed : {3u, 11u}) {
      if (perturbation == 0.0 && seed != 3u) continue;
      record("membrane", perturbation, seed, flat_patch(true, false, perturbation, seed, 0.7));
      record("bending", perturbation, seed, flat_patch(false, true, perturbation, seed, 0.7));
      record("membrane and bending", perturbation, seed,
             flat_patch(true, true, perturbation, seed, -1.3));
    }
  }
  record("rigid motion, cylinder panel", 0.0, 0u, curved_rigid(ShellShape::Cylinder));
  record("rigid motion, sphere zone", 0.0, 0u, curved_rigid(ShellShape::Sphere));
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification (exact states)"));
  block.set("cases", cases);
  block.set("note",
            json::Value::make_string(
                "A 1 x 0.8 m plate of 5 x 4 cells, t = 0.02 m, E = 70 GPa, nu = 0.3, its "
                "interior nodes moved by up to the stated fraction of a cell, turned out of "
                "every coordinate plane; its boundary nodes carry the exact field (all six "
                "DOFs): u = A x + c in the plane (the drilling rotation the in-plane "
                "rotation of A), w quadratic with its rotations (constant curvatures), or "
                "both. Interior displacements and rotations over the largest exact value; "
                "N, M and Q at every element centre against the exact constant resultants "
                "(Q = 0), in each element's own frame. Rigid motion: a cylinder panel and a "
                "sphere zone of 5 x 4 cells with a small rotation and translation on their "
                "boundary nodes: every node follows it, and no element carries a resultant "
                "(against E t |omega|)."));
  summary.set("shell_patch", block);

  StudyOutcome outcome;
  outcome.name = "shell patch tests: membrane, bending, curved rigid motion";
  outcome.kind = "verification";
  outcome.metric = "largest relative error of displacements and resultants";
  outcome.value = worst;
  outcome.tolerance = 1.0e-10;
  outcome.passed = worst <= outcome.tolerance;
  std::ostringstream note;
  note << "largest error " << fmt(worst, 3) << " over 17 cases on distorted, turned meshes";
  outcome.note = note.str();
  return outcome;
}

StudyOutcome study_shell_plate(const std::string& out_dir, json::Value& summary) {
  const Scalar a = 1.0;
  const Scalar e_mod = 1.0e9;
  const Scalar nu = 0.3;
  const Scalar q = 1.0e3;
  const IsotropicMaterial material(e_mod, nu, 0.0);
  const Scalar navier = navier_centre();
  const Scalar marcus = marcus_centre();
  const Scalar g = e_mod / (2.0 * (1.0 + nu));
  // Taylor and Govindjee (2004): the clamped square plate under a uniform
  // load, w_max D / (q a^4) = 0.001265319.
  const Scalar clamped = 0.001265319;
  const std::vector<Index> ladder{4, 8, 16, 32, 64};
  CsvWriter csv(path_join(out_dir, "shell_plate.csv"),
                {"support", "t/a", "mesh", "n", "w_centre[m]", "w_reference[m]",
                 "relative_error[-]", "observed_order"});
  json::Value cases = json::Value::make_array();
  Scalar worst_fine = 0.0;
  Scalar worst_order = kInf;
  Scalar locking_ratio = 0.0;
  std::vector<Scalar> error_at_16;
  for (const std::string support : {std::string("simple"), std::string("clamped")}) {
    const bool clamped_plate = support == "clamped";
    const std::vector<Scalar> thicknesses =
        clamped_plate ? std::vector<Scalar>{1.0e-3} : std::vector<Scalar>{1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4};
    for (const Scalar t : thicknesses) {
      for (const Scalar perturbation : {0.0, 0.2}) {
        const Scalar d_b = e_mod * t * t * t / (12.0 * (1.0 - nu * nu));
        const Scalar reference =
            clamped_plate ? clamped * q * a * a * a * a / d_b
                          : navier * q * a * a * a * a / d_b + marcus * q * a * a / (kShear * g * t);
        json::Value c = json::Value::make_object();
        c.set("support", json::Value::make_string(support));
        c.set("t_over_a", json::Value::make_number(t / a));
        c.set("mesh", json::Value::make_string(perturbation > 0.0 ? "distorted" : "regular"));
        c.set("reference_m", json::Value::make_number(reference));
        json::Value errors = json::Value::make_array();
        Scalar previous = 0.0;
        for (const Index n : ladder) {
          FemModel model = plate(n, a, a, t, material, perturbation);
          if (clamped_plate) {
            model.constraints().push_back(
                fix(region("edges", {box(0.0, 0.0), box(a, a), box(-kInf, kInf, 0.0, 0.0),
                                     box(-kInf, kInf, a, a)}),
                    {0, 1, 2, 3, 4, 5}));
          } else {
            simple_supports(model, a, a, true);
          }
          model.load_case_specs() = {uniform_pressure(q)};
          const StaticSolution sol = solve_static(model);
          const Index centre = (n / 2) * (n + 1) + n / 2;
          const Scalar w = -sol.displacement(6 * centre + 2);
          const Scalar error = std::abs(w - reference) / reference;
          const Scalar p = previous > 0.0 ? order(previous, error) : 0.0;
          csv.raw_row({support, fmt(t, 3), perturbation > 0.0 ? "distorted" : "regular",
                       std::to_string(n), fmt(w, 10), fmt(reference, 10), fmt(error, 4),
                       previous > 0.0 ? fmt(p, 4) : ""});
          errors.push_back(json::Value::make_number(error));
          if (n == ladder.back()) {
            worst_fine = std::max(worst_fine, error);
            if (perturbation == 0.0) worst_order = std::min(worst_order, p);
          }
          if (n == 16 && !clamped_plate && perturbation == 0.0 && t <= 1.0e-2) {
            error_at_16.push_back(error);
          }
          previous = error;
        }
        c.set("n", json::Value::make_number(0.0));
        c.set("errors", errors);
        cases.push_back(c);
      }
    }
  }
  csv.close();
  // MITC4 is free of shear locking on parallelogram meshes; on distorted ones
  // a mesh with few interior nodes can lock as t / a falls - measured on the
  // clamped plate, whose edges hold every DOF.
  CsvWriter distortion(path_join(out_dir, "shell_plate_distortion.csv"),
                       {"n", "t/a", "w_over_reference[-]"});
  json::Value distorted_rows = json::Value::make_array();
  Scalar coarse_thin = 0.0;
  Scalar fine_spread = 0.0;
  for (const Index n : {4, 8, 16}) {
    Scalar lo = kInf;
    Scalar hi = 0.0;
    for (const Scalar t : {1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4}) {
      FemModel model = plate(n, a, a, t, material, 0.2);
      model.constraints().push_back(
          fix(region("edges", {box(0.0, 0.0), box(a, a), box(-kInf, kInf, 0.0, 0.0),
                               box(-kInf, kInf, a, a)}),
              {0, 1, 2, 3, 4, 5}));
      model.load_case_specs() = {uniform_pressure(q)};
      const StaticSolution sol = solve_static(model);
      const Scalar d_b = e_mod * t * t * t / (12.0 * (1.0 - nu * nu));
      const Scalar ratio =
          -sol.displacement(6 * ((n / 2) * (n + 1) + n / 2) + 2) / (clamped * q * a * a * a * a / d_b);
      distortion.raw_row({std::to_string(n), fmt(t, 3), fmt(ratio, 6)});
      json::Value row = json::Value::make_object();
      row.set("n", json::Value::make_number(static_cast<Scalar>(n)));
      row.set("t_over_a", json::Value::make_number(t));
      row.set("w_over_reference", json::Value::make_number(ratio));
      distorted_rows.push_back(row);
      if (n == 4 && t == 1.0e-3) coarse_thin = ratio;
      if (t <= 1.0e-2) {
        lo = std::min(lo, ratio);
        hi = std::max(hi, ratio);
      }
    }
    if (n == 16) fine_spread = hi / lo - 1.0;
  }
  distortion.close();
  // No shear locking: at 16 x 16 the error of the thin plates, whose
  // solutions all approach Kirchhoff's, is the same whatever t / a.
  const Scalar low = *std::min_element(error_at_16.begin(), error_at_16.end());
  const Scalar high = *std::max_element(error_at_16.begin(), error_at_16.end());
  locking_ratio = high / low;

  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification (exact Reissner-Mindlin and "
                                             "thin-plate deflections)"));
  block.set("cases", cases);
  block.set("navier_coefficient", json::Value::make_number(navier));
  block.set("marcus_coefficient", json::Value::make_number(marcus));
  block.set("error_ratio_over_thickness_at_16", json::Value::make_number(locking_ratio));
  block.set("distorted_clamped", distorted_rows);
  block.set("distorted_clamped_4x4_thin_ratio", json::Value::make_number(coarse_thin));
  block.set("distorted_clamped_16x16_thin_spread", json::Value::make_number(fine_spread));
  block.set("note",
            json::Value::make_string(
                "Square plate a = 1 m, E = 1 GPa, nu = 0.3, uniform pressure 1 kPa, n x n "
                "cells (n = 4 ... 64), regular or with interior nodes moved by up to 0.2 of "
                "a cell. Simple support (hard: w and the rotation along the edge held): "
                "the exact Reissner-Mindlin centre deflection w_K + phi / (k G t), w_K from "
                "Navier's series and phi the Marcus moment (grad^2 phi = -q, phi = 0 on "
                "the edges), for t / a = 1e-1 ... 1e-4. Clamped (every DOF of the edges "
                "held): t / a = 1e-3 against the thin-plate value 0.001265319 q a^4 / D "
                "(Taylor and Govindjee 2004); the shear deflection of that plate is of "
                "order (t / a)^2 ~ 1e-6 of it. The distorted meshes keep their centre node "
                "at the centre. Observed orders from successive meshes; "
                "error_ratio_over_thickness_at_16 compares the thin plates (t / a = 1e-2 "
                "... 1e-4) at 16 x 16 - the t / a = 0.1 plate has a different solution, its "
                "shear deflection 5 % of it. distorted_clamped: the clamped plate on the "
                "distorted meshes for t / a = 1e-1 ... 1e-4 (w over the thin-plate value): "
                "MITC4's freedom from shear locking holds on parallelogram meshes; on a "
                "distorted 4 x 4 mesh of the clamped plate, whose 24 interior edges' shear "
                "constraints nearly exhaust its 27 interior DOFs, it locks as t / a falls; "
                "from 8 x 8 on it no longer falls as t / a does."));
  summary.set("shell_plate", block);

  StudyOutcome outcome;
  outcome.name = "shell plates: simply supported (t/a 1e-1 to 1e-4) and clamped";
  outcome.kind = "verification";
  outcome.metric = "largest centre-deflection error at 64 x 64";
  outcome.value = worst_fine;
  outcome.tolerance = 5.0e-4;
  outcome.passed = worst_fine <= outcome.tolerance && worst_order >= 1.8 &&
                   locking_ratio < 1.01 && fine_spread < 1.0e-2;
  std::ostringstream note;
  note << "smallest observed order " << fmt(worst_order, 3)
       << " (regular meshes, 32 -> 64); at 16 x 16 the errors for t/a = 1e-2 ... 1e-4 differ "
       << "by a factor " << fmt(locking_ratio, 4) << " (no shear locking); a distorted 4 x 4 "
       << "clamped mesh locks (" << fmt(coarse_thin, 3) << " of w at t/a = 1e-3), "
       << "16 x 16 does not (spread " << fmt(fine_spread, 3) << " over t/a = 1e-2 ... 1e-4)";
  outcome.note = note.str();
  return outcome;
}

StudyOutcome study_shell_plate_modes(const std::string& out_dir, json::Value& summary) {
  const Scalar a = 1.0;
  const Scalar t = 0.01;
  const Scalar e_mod = 70.0e9;
  const Scalar nu = 0.3;
  const Scalar rho = 2700.0;
  const IsotropicMaterial material(e_mod, nu, rho);
  // The exact frequencies of the modes (m, n), m, n = 1 ... 4, ascending.
  struct Exact {
    int m;
    int n;
    Scalar omega;
  };
  std::vector<Exact> exact;
  for (int m = 1; m <= 4; ++m) {
    for (int n = 1; n <= 4; ++n) {
      const TrigMode mode = trig_mode(m * kPi / a, n * kPi / a, e_mod, nu, t, rho, 0.0, 0.0);
      exact.push_back({m, n, std::sqrt(smallest(mode.k, mode.m))});
    }
  }
  std::sort(exact.begin(), exact.end(),
            [](const Exact& x, const Exact& y) { return x.omega < y.omega; });
  const int modes = 6;
  const std::vector<Index> ladder{8, 16, 32, 64};
  CsvWriter csv(path_join(out_dir, "shell_plate_modes.csv"),
                {"mass", "n", "mode", "m", "n_half_waves", "frequency[Hz]", "exact[Hz]",
                 "relative_error[-]", "observed_order"});
  Scalar worst_fine = 0.0;
  Scalar worst_lumped = 0.0;
  Scalar worst_order = kInf;
  json::Value meshes = json::Value::make_array();
  for (const MassType mass : {MassType::Consistent, MassType::Lumped}) {
    const std::string mass_name = mass == MassType::Consistent ? "consistent" : "lumped";
    std::vector<Scalar> previous(modes, 0.0);
    for (const Index n : ladder) {
      FemModel model = plate(n, a, a, t, material);
      simple_supports(model, a, a, true);
      model.load_case_specs() = {uniform_pressure(1.0)};
      model.finalize();
      Assembler assembler(model);
      ModalAnalysisOptions options;
      options.num_modes = modes;
      options.mass_type = mass;
      const ModalResult r = solve_modal(model, assembler, options);
      json::Value entry = json::Value::make_object();
      entry.set("mass", json::Value::make_string(mass_name));
      entry.set("n", json::Value::make_number(static_cast<Scalar>(n)));
      json::Value errors = json::Value::make_array();
      for (int i = 0; i < modes; ++i) {
        const Scalar error =
            std::abs(r.angular_frequencies(i) - exact[i].omega) / exact[i].omega;
        const Scalar p = previous[i] > 0.0 ? order(previous[i], error) : 0.0;
        csv.raw_row({mass_name, std::to_string(n), std::to_string(i + 1),
                     std::to_string(exact[i].m), std::to_string(exact[i].n),
                     fmt(r.frequencies_hz(i), 10), fmt(exact[i].omega / (2.0 * kPi), 10),
                     fmt(error, 4), previous[i] > 0.0 ? fmt(p, 4) : ""});
        errors.push_back(json::Value::make_number(error));
        if (n == ladder.back()) {
          if (mass == MassType::Consistent) {
            worst_fine = std::max(worst_fine, error);
          } else {
            worst_lumped = std::max(worst_lumped, error);
          }
          worst_order = std::min(worst_order, p);
        }
        previous[i] = error;
      }
      entry.set("errors", errors);
      meshes.push_back(entry);
    }
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification (exact Reissner-Mindlin "
                                             "frequencies)"));
  block.set("meshes", meshes);
  block.set("note",
            json::Value::make_string(
                "Simply supported (hard) square plate a = 1 m, t = 0.01 m, E = 70 GPa, "
                "nu = 0.3, rho = 2700 kg/m^3, consistent mass; the six lowest frequencies "
                "- the modes (1,1), (1,2) and (2,1), (2,2), (1,3) and (3,1) - against the "
                "exact Reissner-Mindlin frequencies with the rotary inertia rho t^3 / 12 and "
                "the shear factor 5/6 (the smallest root of the 3 x 3 modal problem in "
                "W sin sin, psi cos sin, sin cos). The rotations about the normal carry no "
                "mass; the eigensolver treats that semi-definite mass. Consistent mass and "
                "the lumped one (Hinton-Rock-Zienkiewicz: per DOF component the diagonal "
                "scaled to the element's total, rotary inertia included)."));
  block.set("largest_error_at_64_consistent", json::Value::make_number(worst_fine));
  block.set("largest_error_at_64_lumped", json::Value::make_number(worst_lumped));
  summary.set("shell_plate_modes", block);
  StudyOutcome outcome;
  outcome.name = "shell plate: six lowest natural frequencies, consistent and lumped mass";
  outcome.kind = "verification";
  outcome.metric = "largest frequency error at 64 x 64, consistent mass";
  outcome.value = worst_fine;
  outcome.tolerance = 3.0e-3;
  const Scalar lumped_tolerance = 1.0e-3;
  outcome.passed = worst_fine <= outcome.tolerance && worst_lumped <= lumped_tolerance &&
                   worst_order >= 1.9;
  std::ostringstream note;
  note << "lumped mass " << fmt(worst_lumped, 3) << " (tolerance " << fmt(lumped_tolerance, 3)
       << "); smallest observed order " << fmt(worst_order, 3) << " (32 -> 64, both masses)";
  outcome.note = note.str();
  return outcome;
}

// ---------------------------------------------------------------------------
// Harmonic response
// ---------------------------------------------------------------------------

/// The exact steady harmonic centre deflection of the simply supported
/// (hard) Reissner-Mindlin plate a x a under a uniform pressure q cos(omega t)
/// with the loss factor eta on its stiffness: the sum over the odd modes
/// (m, n) that the uniform load excites, 16 q / (pi^2 m n) each, of the
/// deflection amplitude W solving the mode's 3 x 3 problem
/// (k (1 + i eta) - omega^2 m) (W, P1, P2) = (16 q / (pi^2 m n), 0, 0).
std::complex<Scalar> harmonic_centre(Scalar a, Scalar e, Scalar nu, Scalar t, Scalar rho,
                                     Scalar q, Scalar omega, Scalar eta) {
  using Complex = std::complex<Scalar>;
  Complex sum = 0.0;
  for (int m = 1; m <= 401; m += 2) {
    for (int n = 1; n <= 401; n += 2) {
      const TrigMode mode = trig_mode(m * kPi / a, n * kPi / a, e, nu, t, rho, 0.0, 0.0);
      const Eigen::Matrix3cd dynamic = mode.k.cast<Complex>() * Complex(1.0, eta) -
                                       (omega * omega) * mode.m.cast<Complex>();
      Eigen::Vector3cd load = Eigen::Vector3cd::Zero();
      load(0) = 16.0 * q / (kPi * kPi * static_cast<Scalar>(m) * static_cast<Scalar>(n));
      const Eigen::Vector3cd x = dynamic.fullPivLu().solve(load);
      const Scalar sign = (((m + n) / 2 - 1) % 2 == 0) ? 1.0 : -1.0;
      sum += sign * x(0);
    }
  }
  return sum;
}

StudyOutcome study_shell_plate_harmonic(const std::string& out_dir, json::Value& summary) {
  const Scalar a = 1.0;
  const Scalar t = 0.01;
  const Scalar e_mod = 70.0e9;
  const Scalar nu = 0.3;
  const Scalar rho = 2700.0;
  const Scalar q = 1.0e3;
  const IsotropicMaterial material(e_mod, nu, rho);
  // Static, below the lowest resonance the uniform load excites - the mode
  // (1, 1) at 48.4 Hz - and between it and the next, (1, 3) and (3, 1) at
  // 241.6 Hz; undamped and with a loss factor.
  const std::vector<Scalar> frequencies{0.0, 30.0, 100.0, 200.0};
  const std::vector<Scalar> losses{0.0, 0.05};
  const std::vector<Index> ladder{8, 16, 32, 64};
  // The series against the closed form of the static deflection (Navier's
  // series plus the Marcus moment over the shear stiffness).
  const Scalar d_b = e_mod * t * t * t / (12.0 * (1.0 - nu * nu));
  const Scalar g = e_mod / (2.0 * (1.0 + nu));
  const Scalar closed_form =
      navier_centre() * q * a * a * a * a / d_b + marcus_centre() * q * a * a / (kShear * g * t);
  const Scalar series_check =
      std::abs(harmonic_centre(a, e_mod, nu, t, rho, q, 0.0, 0.0) - closed_form) / closed_form;
  CsvWriter csv(path_join(out_dir, "shell_plate_harmonic.csv"),
                {"mass", "loss_factor", "frequency[Hz]", "n", "w_real[m]", "w_imag[m]",
                 "exact_real[m]", "exact_imag[m]", "relative_error[-]", "observed_order"});
  Scalar worst_fine = 0.0;
  Scalar worst_order = kInf;
  json::Value cases = json::Value::make_array();
  for (const MassType mass : {MassType::Consistent, MassType::Lumped}) {
    const std::string mass_name = mass == MassType::Consistent ? "consistent" : "lumped";
    for (const Scalar eta : losses) {
      std::vector<std::complex<Scalar>> exact;
      for (const Scalar f : frequencies) {
        exact.push_back(harmonic_centre(a, e_mod, nu, t, rho, q, 2.0 * kPi * f, eta));
      }
      std::vector<Scalar> previous(frequencies.size(), 0.0);
      json::Value c = json::Value::make_object();
      c.set("mass", json::Value::make_string(mass_name));
      c.set("loss_factor", json::Value::make_number(eta));
      json::Value meshes = json::Value::make_array();
      for (const Index n : ladder) {
        FemModel model = plate(n, a, a, t, material);
        simple_supports(model, a, a, true);
        model.load_case_specs() = {uniform_pressure(q)};
        model.finalize();
        Assembler assembler(model);
        FrequencyResponseOptions options;
        options.frequencies = frequencies;
        options.mass_type = mass;
        options.structural_damping = eta;
        options.snapshot_frequencies = frequencies;
        const FrequencyResponseResult r = solve_frequency_response(model, assembler, 0, options);
        if (r.snapshots.size() != frequencies.size()) {
          throw SolverError("shell-plate-harmonic: expected one snapshot per frequency");
        }
        const Index centre = (n / 2) * (n + 1) + n / 2;
        json::Value errors = json::Value::make_array();
        for (std::size_t j = 0; j < frequencies.size(); ++j) {
          // The pressure pushes against the normal (+z): the plate moves down.
          const std::complex<Scalar> w = -r.snapshots[j].displacement(model.dofs().dof(centre, 2));
          const Scalar error = std::abs(w - exact[j]) / std::abs(exact[j]);
          const Scalar p = previous[j] > 0.0 ? order(previous[j], error) : 0.0;
          csv.raw_row({mass_name, fmt(eta, 3), fmt(frequencies[j], 6), std::to_string(n),
                       fmt(w.real(), 10), fmt(w.imag(), 10), fmt(exact[j].real(), 10),
                       fmt(exact[j].imag(), 10), fmt(error, 4),
                       previous[j] > 0.0 ? fmt(p, 4) : ""});
          errors.push_back(json::Value::make_number(error));
          if (n == ladder.back()) {
            worst_fine = std::max(worst_fine, error);
            worst_order = std::min(worst_order, p);
          }
          previous[j] = error;
        }
        json::Value entry = json::Value::make_object();
        entry.set("n", json::Value::make_number(static_cast<Scalar>(n)));
        entry.set("errors", errors);
        meshes.push_back(entry);
      }
      c.set("meshes", meshes);
      cases.push_back(c);
    }
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification (exact Reissner-Mindlin harmonic "
                                             "response)"));
  block.set("cases", cases);
  block.set("series_vs_closed_form_static", json::Value::make_number(series_check));
  block.set("note",
            json::Value::make_string(
                "The simply supported (hard) square plate of the frequency study (a = 1 m, "
                "t = 0.01 m, E = 70 GPa, nu = 0.3, rho = 2700 kg/m^3) under a uniform "
                "pressure of 1 kPa cos(omega t), solved directly at 0, 30, 100 and 200 Hz - "
                "static, below the lowest resonance the load excites, (1, 1) at 48.4 Hz, and "
                "between it and (1, 3), (3, 1) at 241.6 Hz - undamped and with the loss "
                "factor 0.05, on consistent and lumped mass. The exact centre amplitude is "
                "the series over the odd modes (m, n <= 401) of the 3 x 3 problems "
                "(k (1 + i eta) - omega^2 m) x = (16 q / (pi^2 m n), 0, 0); "
                "series_vs_closed_form_static checks it at 0 Hz against Navier's series "
                "plus the Marcus moment. Errors: |w - w_exact| / |w_exact| of the complex "
                "centre amplitude."));
  summary.set("shell_plate_harmonic", block);
  StudyOutcome outcome;
  outcome.name = "shell plate: harmonic response vs the exact Reissner-Mindlin series";
  outcome.kind = "verification";
  outcome.metric = "largest centre-amplitude error at 64 x 64";
  outcome.value = worst_fine;
  // The largest error is at 200 Hz, where the modes (1, 1) and (1, 3)
  // nearly cancel at the centre (the amplitude is a seventh of the static
  // one), so the small error of the (1, 3) resonance is a large part of it.
  outcome.tolerance = 1.0e-2;
  outcome.passed =
      worst_fine <= outcome.tolerance && worst_order >= 1.9 && series_check <= 1.0e-9;
  std::ostringstream note;
  note << "smallest observed order " << fmt(worst_order, 3)
       << " (32 -> 64; consistent and lumped mass, eta = 0 and 0.05, 0 ... 200 Hz); series "
          "vs closed form at 0 Hz "
       << fmt(series_check, 2);
  outcome.note = note.str();
  return outcome;
}

StudyOutcome study_shell_plate_buckling(const std::string& out_dir, json::Value& summary) {
  const Scalar a = 1.0;
  const Scalar t = 0.01;
  const Scalar e_mod = 70.0e9;
  const Scalar nu = 0.3;
  const IsotropicMaterial material(e_mod, nu, 2700.0);
  const Scalar d_b = e_mod * t * t * t / (12.0 * (1.0 - nu * nu));
  const Scalar sigma0 = 1.0e6;  // applied compressive stress [Pa]
  const std::vector<Index> ladder{8, 16, 32, 64};
  CsvWriter csv(path_join(out_dir, "shell_plate_buckling.csv"),
                {"loading", "n", "critical_force[N/m]", "exact[N/m]", "k_kirchhoff[-]",
                 "relative_error[-]", "observed_order"});
  json::Value cases = json::Value::make_array();
  Scalar worst_fine = 0.0;
  Scalar worst_order = kInf;
  for (const bool biaxial : {false, true}) {
    // The exact critical force of the model: the smallest load factor over
    // the modes (m, n) of the 3 x 3 problem.
    Scalar exact = kInf;
    for (int m = 1; m <= 4; ++m) {
      for (int n = 1; n <= 4; ++n) {
        const TrigMode mode =
            trig_mode(m * kPi / a, n * kPi / a, e_mod, nu, t, 2700.0, 1.0, biaxial ? 1.0 : 0.0);
        exact = std::min(exact, smallest(mode.k, mode.g));
      }
    }
    const Scalar kirchhoff = (biaxial ? 2.0 : 4.0) * kPi * kPi * d_b / (a * a);
    json::Value c = json::Value::make_object();
    c.set("loading", json::Value::make_string(biaxial ? "equal biaxial" : "uniaxial"));
    c.set("exact_N_per_m", json::Value::make_number(exact));
    c.set("kirchhoff_N_per_m", json::Value::make_number(kirchhoff));
    json::Value errors = json::Value::make_array();
    Scalar previous = 0.0;
    for (const Index n : ladder) {
      FemModel model = plate(n, a, a, t, material);
      simple_supports(model, a, a, false);
      // In-plane: u_x held on x = 0 and u_y on y = 0 (biaxial) or at the
      // corner (uniaxial), so the compression is uniform.
      model.constraints().push_back(fix(region("x0", {box(0.0, 0.0)}), {0}));
      model.constraints().push_back(
          fix(region("y0", {biaxial ? box(-kInf, kInf, 0.0, 0.0) : box(0.0, 0.0, 0.0, 0.0)}),
              {1}));
      LoadCaseSpec lc;
      lc.name = "compression";
      TractionLoadSpec tx;
      tx.region = region("x_a", {box(a, a)});
      tx.traction = Vector3(-sigma0, 0.0, 0.0);
      lc.tractions.push_back(tx);
      if (biaxial) {
        TractionLoadSpec ty;
        ty.region = region("y_a", {box(-kInf, kInf, a, a)});
        ty.traction = Vector3(0.0, -sigma0, 0.0);
        lc.tractions.push_back(ty);
      }
      model.load_case_specs() = {lc};
      model.finalize();
      Assembler assembler(model);
      BucklingOptions options;
      options.num_modes = 2;
      const BucklingResult r = analyse_buckling(model, assembler, 0, options);
      const Scalar critical = r.load_factors(0) * sigma0 * t;
      const Scalar error = std::abs(critical - exact) / exact;
      const Scalar p = previous > 0.0 ? order(previous, error) : 0.0;
      csv.raw_row({biaxial ? "equal biaxial" : "uniaxial", std::to_string(n), fmt(critical, 10),
                   fmt(exact, 10), fmt(critical * a * a / (kPi * kPi * d_b), 8), fmt(error, 4),
                   previous > 0.0 ? fmt(p, 4) : ""});
      errors.push_back(json::Value::make_number(error));
      if (n == ladder.back()) {
        worst_fine = std::max(worst_fine, error);
        worst_order = std::min(worst_order, p);
      }
      previous = error;
    }
    c.set("errors", errors);
    cases.push_back(c);
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification (exact buckling loads of the "
                                             "model)"));
  block.set("cases", cases);
  block.set("note",
            json::Value::make_string(
                "Simply supported (hard) square plate a = 1 m, t = 0.01 m, E = 70 GPa, "
                "nu = 0.3, compressed by a uniform edge stress (u_x held on x = 0; u_y on "
                "y = 0 for the biaxial case, at one corner otherwise). Exact: the smallest "
                "load factor over the modes (m, n) of the 3 x 3 problem of the "
                "Reissner-Mindlin plate with the degenerated solid's geometric stiffness, "
                "N [w_,x^2 + (t^2 / 12)(psi_1,x^2 + psi_2,x^2)] (and the same along y). It "
                "lies below Kirchhoff's k = 4 (uniaxial) and k = 2 (biaxial) by the shear "
                "deformation, N / (2 k G t) and N / (k G t) = 5.6e-4, and by the rotations' "
                "geometric term, 1.6e-4; k_kirchhoff is the computed load in units of "
                "pi^2 D / a^2."));
  summary.set("shell_plate_buckling", block);
  StudyOutcome outcome;
  outcome.name = "shell plate buckling: uniaxial (k = 4) and biaxial (k = 2) compression";
  outcome.kind = "verification";
  outcome.metric = "largest critical-load error at 64 x 64";
  outcome.value = worst_fine;
  outcome.tolerance = 1.0e-3;
  outcome.passed = worst_fine <= outcome.tolerance && worst_order >= 1.9;
  std::ostringstream note;
  note << "smallest observed order " << fmt(worst_order, 3) << " (32 -> 64)";
  outcome.note = note.str();
  return outcome;
}

StudyOutcome study_shell_cylinder_pressure(const std::string& out_dir, json::Value& summary) {
  const Scalar radius = 1.0;
  const Scalar length = 2.0;
  const Scalar t = 0.01;
  const Scalar e_mod = 70.0e9;
  const Scalar nu = 0.3;
  const Scalar p = 1.0e5;
  // The exact state of the shell's continuum model: each fibre moves out by
  // u_r and stays radial, so the hoop strain is u_r / (R + z) through the
  // thickness; the axial force and the pressure's balance, N = p R, then give
  // u_r = p R / (E ln((R + t/2) / (R - t/2))) - within (t/R)^2 / 12 of the
  // membrane formula p R^2 / (E t) - and the axial strain -nu u_r ln(..) / t.
  const Scalar log_ratio = std::log((radius + 0.5 * t) / (radius - 0.5 * t));
  const Scalar u_r = p * radius / (e_mod * log_ratio);
  const Scalar eps_z = -nu * u_r * log_ratio / t;
  const Scalar membrane = p * radius * radius / (e_mod * t);
  const std::vector<Index> ladder{8, 16, 32, 64, 128, 256};
  CsvWriter csv(path_join(out_dir, "shell_cylinder_pressure.csv"),
                {"n_around", "radial_error[-]", "axial_error[-]", "hoop_force_error[-]",
                 "observed_order"});
  json::Value errors = json::Value::make_array();
  Scalar previous = 0.0;
  Scalar fine = 0.0;
  Scalar fine_order = 0.0;
  Scalar worst_hoop = 0.0;
  for (const Index n : ladder) {
    ShellMeshSpec spec;
    spec.shape = ShellShape::Cylinder;
    spec.radius = radius;
    spec.length = length;
    spec.n1 = n;
    spec.n2 = 4;
    FemModel model(make_structured_shell_mesh(spec), IsotropicMaterial(e_mod, nu, 2700.0), t,
                   StressState::Shell, IntegrationOptions());
    const Scalar tol = 1.0e-9;
    // A slice of a long cylinder: u_z held on the ring z = 0, and the end
    // rings' rotations about x and y held - symmetry planes, which the exact
    // state satisfies (its rotations are zero); u_y at (+-R, 0, 0) and u_x at
    // (0, R, 0), zero in the exact state, leave no rigid motion.
    const SelectorGroup base = region("base", {box(-kInf, kInf, -kInf, kInf, 0.0, 0.0)});
    model.constraints().push_back(fix(base, {2}));
    model.constraints().push_back(
        fix(region("ends", {box(-kInf, kInf, -kInf, kInf, 0.0, 0.0),
                            box(-kInf, kInf, -kInf, kInf, length, length)}),
            {3, 4}));
    model.constraints().push_back(
        fix(region("y_points", {box(radius - tol, radius + tol, -tol, tol, 0.0, 0.0),
                                box(-radius - tol, -radius + tol, -tol, tol, 0.0, 0.0)}),
            {1}));
    model.constraints().push_back(
        fix(region("x_point", {box(-tol, tol, radius - tol, radius + tol, 0.0, 0.0)}), {0}));
    // Internal pressure pushes against the outward normal's side: -p.
    model.load_case_specs() = {uniform_pressure(-p)};
    const StaticSolution sol = solve_static(model);
    Scalar radial = 0.0;
    Scalar axial = 0.0;
    for (Index k = 0; k < model.mesh().num_nodes(); ++k) {
      const Vector3 x = model.mesh().node(k);
      const Vector3 u = sol.displacement.segment<3>(6 * k);
      const Vector3 e_r = Vector3(x.x(), x.y(), 0.0).normalized();
      radial = std::max(radial, std::abs(u.dot(e_r) - u_r) / u_r);
      axial = std::max(axial, std::abs(u.z() - eps_z * x.z()) / u_r);
    }
    Assembler assembler(model);
    const ShellField field = recover_shell_resultants(model, assembler, sol.displacement);
    Scalar hoop = 0.0;
    for (const ShellResultants& r : field.element) {
      // The hoop force: the normal force along e_3 x e_z.
      const Vector3 hoop_dir = r.e3.cross(Vector3::UnitZ()).normalized();
      const Matrix3 n_tensor = tensor(r.membrane, r.e1, r.e2);
      hoop = std::max(hoop, std::abs(hoop_dir.dot(n_tensor * hoop_dir) - p * radius) / (p * radius));
    }
    const Scalar ord = previous > 0.0 ? order(previous, radial) : 0.0;
    csv.raw_row({std::to_string(n), fmt(radial, 4), fmt(axial, 4), fmt(hoop, 4),
                 previous > 0.0 ? fmt(ord, 4) : ""});
    errors.push_back(json::Value::make_number(radial));
    previous = radial;
    fine = radial;
    fine_order = ord;
    worst_hoop = hoop;
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification (exact state of the model)"));
  block.set("radial_errors", errors);
  block.set("u_r_m", json::Value::make_number(u_r));
  block.set("membrane_formula_m", json::Value::make_number(membrane));
  block.set("note",
            json::Value::make_string(
                "A slice of a long cylinder, R = 1 m, 2 m long, t = 0.01 m, E = 70 GPa, "
                "nu = 0.3, internal pressure 0.1 MPa on the mid-surface, n cells round it "
                "and 4 along, its end rings held against rotation (symmetry planes). The "
                "shell's fibres diverge through its thickness, so its exact state is that "
                "of a thick ring in plane stress: u_r = p R / (E ln((R + t/2)/(R - t/2))), "
                "(t/R)^2 / 12 = 8.3e-6 below the membrane formula p R^2 / (E t), the axial "
                "strain -nu u_r ln(..) / t and the hoop force p R. Errors: the largest "
                "nodal ones over u_r, and the hoop force at the element centres over p R; "
                "the flat facets of the polygon converge to the circle at second order. "
                "With free ends instead the facets' bending under the pressure leaves a "
                "boundary layer about sqrt(R t) wide at each end, which a mesh this coarse "
                "along the axis does not resolve."));
  summary.set("shell_cylinder_pressure", block);
  StudyOutcome outcome;
  outcome.name = "shell cylinder under internal pressure: the thick-ring membrane state";
  outcome.kind = "verification";
  outcome.metric = "largest radial-displacement error at 256 cells round";
  outcome.value = fine;
  outcome.tolerance = 2.0e-4;
  outcome.passed = fine <= outcome.tolerance && fine_order >= 1.9;
  std::ostringstream note;
  note << "observed order " << fmt(fine_order, 3) << " (128 -> 256); hoop force error "
       << fmt(worst_hoop, 3);
  outcome.note = note.str();
  return outcome;
}

namespace {

/// One MacNeal-Harder benchmark on a mesh of n x n cells: the displacement
/// the benchmark reports.
using Benchmark = std::function<Scalar(Index)>;

Scalar scordelis_lo(Index n) {
  // Quarter of the roof: R = 25, L / 2 = 25, 40 degrees each side of the
  // crown, t = 0.25, E = 4.32e8, nu = 0, self-weight 90 per unit area.
  ShellMeshSpec spec;
  spec.shape = ShellShape::Cylinder;
  spec.axis = 0;
  spec.radius = 25.0;
  spec.length = 25.0;
  spec.angle_start = 50.0;
  spec.angle_end = 90.0;
  spec.n1 = n;
  spec.n2 = n;
  FemModel model(make_structured_shell_mesh(spec), IsotropicMaterial(4.32e8, 0.0, 1.0), 0.25,
                 StressState::Shell, IntegrationOptions());
  model.constraints().push_back(fix(region("symmetry_x0", {box(0.0, 0.0)}), {0, 4, 5}));
  model.constraints().push_back(fix(region("diaphragm", {box(25.0, 25.0)}), {1, 2}));
  model.constraints().push_back(fix(region("crown", {box(-kInf, kInf, 0.0, 0.0)}), {1, 3, 5}));
  LoadCaseSpec lc;
  lc.name = "self_weight";
  lc.gravity = Vector3(0.0, 0.0, -90.0 / 0.25);
  model.load_case_specs() = {lc};
  const StaticSolution sol = solve_static(model);
  const Scalar theta = 50.0 * kPi / 180.0;
  const Index k = node_at(model.mesh(), Vector3(0.0, 25.0 * std::cos(theta), 25.0 * std::sin(theta)));
  return -sol.displacement(6 * k + 2);
}

Scalar pinched_cylinder(Index n) {
  // An octant: R = 300, L / 2 = 300, t = 3, E = 3e6, nu = 0.3, end
  // diaphragms, a pair of unit loads pinching the middle section.
  ShellMeshSpec spec;
  spec.shape = ShellShape::Cylinder;
  spec.radius = 300.0;
  spec.length = 300.0;
  spec.angle_start = 0.0;
  spec.angle_end = 90.0;
  spec.n1 = n;
  spec.n2 = n;
  FemModel model(make_structured_shell_mesh(spec), IsotropicMaterial(3.0e6, 0.3, 1.0), 3.0,
                 StressState::Shell, IntegrationOptions());
  model.constraints().push_back(fix(region("symmetry_z0", {box(-kInf, kInf, -kInf, kInf, 0.0, 0.0)}), {2, 3, 4}));
  model.constraints().push_back(fix(region("symmetry_y0", {box(-kInf, kInf, 0.0, 0.0)}), {1, 3, 5}));
  model.constraints().push_back(fix(region("symmetry_x0", {box(0.0, 0.0)}), {0, 4, 5}));
  model.constraints().push_back(fix(region("diaphragm", {box(-kInf, kInf, -kInf, kInf, 300.0, 300.0)}), {0, 1}));
  LoadCaseSpec lc;
  lc.name = "pinch";
  PointLoadSpec p;
  p.region = region("load", {box(300.0, 300.0, 0.0, 0.0, 0.0, 0.0)});
  p.force = Vector3(-0.25, 0.0, 0.0);
  lc.point_loads.push_back(p);
  model.load_case_specs() = {lc};
  const StaticSolution sol = solve_static(model);
  return -sol.displacement(6 * node_at(model.mesh(), Vector3(300.0, 0.0, 0.0)));
}

Scalar pinched_hemisphere(Index n) {
  // A quarter: R = 10, t = 0.04, E = 6.825e7, nu = 0.3, an 18 degree hole
  // at the top, loads of 2 in and out at 90 degree intervals on the equator.
  ShellMeshSpec spec;
  spec.shape = ShellShape::Sphere;
  spec.radius = 10.0;
  spec.angle_start = 0.0;
  spec.angle_end = 90.0;
  spec.polar_start = 18.0;
  spec.polar_end = 90.0;
  spec.n1 = n;
  spec.n2 = n;
  FemModel model(make_structured_shell_mesh(spec), IsotropicMaterial(6.825e7, 0.3, 1.0), 0.04,
                 StressState::Shell, IntegrationOptions());
  model.constraints().push_back(fix(region("symmetry_y0", {box(-kInf, kInf, 0.0, 0.0)}), {1, 3, 5}));
  model.constraints().push_back(fix(region("symmetry_x0", {box(0.0, 0.0)}), {0, 4, 5}));
  model.constraints().push_back(fix(region("one_point", {box(10.0, 10.0, 0.0, 0.0, 0.0, 0.0)}), {2}));
  LoadCaseSpec lc;
  lc.name = "pinch";
  PointLoadSpec out;
  out.region = region("load_x", {box(10.0, 10.0, 0.0, 0.0, 0.0, 0.0)});
  out.force = Vector3(1.0, 0.0, 0.0);
  PointLoadSpec in;
  in.region = region("load_y", {box(0.0, 0.0, 10.0, 10.0, 0.0, 0.0)});
  in.force = Vector3(0.0, -1.0, 0.0);
  lc.point_loads = {out, in};
  model.load_case_specs() = {lc};
  const StaticSolution sol = solve_static(model);
  return sol.displacement(6 * node_at(model.mesh(), Vector3(10.0, 0.0, 0.0)));
}

StudyOutcome benchmark_study(const std::string& key, const std::string& title,
                             const Benchmark& run, Scalar reference, const std::string& unit,
                             const std::string& note_text, const std::string& out_dir,
                             json::Value& summary) {
  const std::vector<Index> ladder{4, 8, 16, 32, 64};
  CsvWriter csv(path_join(out_dir, key + ".csv"),
                {"n", "value[" + unit + "]", "reference[" + unit + "]", "normalised[-]"});
  json::Value values = json::Value::make_array();
  Scalar last = 0.0;
  for (const Index n : ladder) {
    const Scalar v = run(n);
    csv.raw_row({std::to_string(n), fmt(v, 10), fmt(reference, 10), fmt(v / reference, 6)});
    json::Value e = json::Value::make_object();
    e.set("n", json::Value::make_number(static_cast<Scalar>(n)));
    e.set("value", json::Value::make_number(v));
    e.set("normalised", json::Value::make_number(v / reference));
    values.push_back(e);
    last = v / reference;
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("validation (thin-shell theory reference)"));
  block.set("reference", json::Value::make_number(reference));
  block.set("meshes", values);
  block.set("note", json::Value::make_string(note_text));
  summary.set(key, block);
  StudyOutcome outcome;
  outcome.name = title;
  outcome.kind = "validation";
  outcome.metric = "|value / reference - 1| at 64 x 64";
  outcome.value = std::abs(last - 1.0);
  outcome.tolerance = 1.0e-2;
  outcome.passed = outcome.value <= outcome.tolerance;
  std::ostringstream note;
  note << "normalised " << fmt(last, 5) << " at 64 x 64";
  outcome.note = note.str();
  return outcome;
}

}  // namespace

StudyOutcome study_shell_scordelis_lo(const std::string& out_dir, json::Value& summary) {
  return benchmark_study(
      "shell_scordelis_lo", "Scordelis-Lo roof: vertical displacement at the free edge",
      scordelis_lo, 0.3024, "m",
      "Quarter of the roof (R = 25, L = 50, 40 degrees each side of the crown, t = 0.25, "
      "E = 4.32e8, nu = 0, self-weight 90 per unit area) on rigid diaphragms, n x n cells; "
      "the vertical displacement at the middle of the free edge against 0.3024, the "
      "deep-shell value MacNeal and Harder (1985) quote (0.3086 in the shallow-shell "
      "theory of Scordelis and Lo). The refinement approaches it from below at first "
      "order at the finest meshes, the free edge's boundary layer resolving slowly.",
      out_dir, summary);
}

StudyOutcome study_shell_pinched_cylinder(const std::string& out_dir, json::Value& summary) {
  return benchmark_study(
      "shell_pinched_cylinder", "pinched cylinder with end diaphragms: displacement under the load",
      pinched_cylinder, 1.8248e-5, "m",
      "An octant of the cylinder (R = 300, L = 600, t = 3, E = 3e6, nu = 0.3) with rigid "
      "end diaphragms, pinched by unit loads at the middle, n x n cells; the displacement "
      "under the load against 1.8248e-5 (MacNeal and Harder 1985, from Flugge's series "
      "solution). MITC4 converges slowly here (membrane locking of the coarse meshes); "
      "the finest meshes pass the reference, as the transverse-shear deflection under a "
      "point load of a shear-deformable shell grows as log(1/h), which thin-shell theory "
      "does not have.",
      out_dir, summary);
}

StudyOutcome study_shell_pinched_hemisphere(const std::string& out_dir, json::Value& summary) {
  return benchmark_study(
      "shell_pinched_hemisphere", "pinched hemisphere with an 18 degree hole: radial displacement",
      pinched_hemisphere, 0.094, "m",
      "A quarter of the hemisphere (R = 10, t = 0.04, E = 6.825e7, nu = 0.3, an 18 degree "
      "hole at the top), loads of 2 alternately in and out at 90 degree intervals on the "
      "equator, n x n cells; the displacement under a load against 0.094 (MacNeal and "
      "Harder 1985).",
      out_dir, summary);
}

namespace {

/// A square tube of mid-surface side b and length L along x: 4 nb cells
/// round it (nb per wall), nl along, its corners folds of 90 degrees; normals
/// outwards.
Mesh box_tube(Scalar b, Scalar length, Index nb, Index nl) {
  const Index np = 4 * nb;
  const Scalar h = 0.5 * b;
  const Eigen::Vector2d corner[4] = {{-h, -h}, {h, -h}, {h, h}, {-h, h}};
  Matrix coords(3, np * (nl + 1));
  for (Index j = 0; j <= nl; ++j) {
    const Scalar x = length * static_cast<Scalar>(j) / static_cast<Scalar>(nl);
    for (Index i = 0; i < np; ++i) {
      const Index side = i / nb;
      const Scalar s = static_cast<Scalar>(i % nb) / static_cast<Scalar>(nb);
      const Eigen::Vector2d p = corner[side] + s * (corner[(side + 1) % 4] - corner[side]);
      coords.col(j * np + i) = Vector3(x, p.x(), p.y());
    }
  }
  std::vector<Index> connectivity;
  for (Index j = 0; j < nl; ++j) {
    for (Index i = 0; i < np; ++i) {
      const Index i1 = (i + 1) % np;
      connectivity.insert(connectivity.end(),
                          {j * np + i, j * np + i1, (j + 1) * np + i1, (j + 1) * np + i});
    }
  }
  Mesh mesh(std::move(coords), std::move(connectivity), ElementType::Shell4);
  mesh.validate();
  return mesh;
}

struct BoxResult {
  Scalar deflection = 0.0;  ///< mean u_z of the tip section [m]
  Scalar twist = 0.0;       ///< mean rotation of the tip section about x [rad]
};

/// The tube clamped at x = 0 and loaded at its tip by a shear force P along
/// -z (a uniform shear on the two side walls) or a torque T (Bredt's shear
/// flow T / (2 b^2) round the walls).
BoxResult box_beam(Index nb, Scalar drilling, bool torsion) {
  const Scalar b = 0.1;
  const Scalar length = 2.0;
  const Scalar t = 0.002;
  const Scalar h = 0.5 * b;
  FemModel model(box_tube(b, length, nb, 20 * nb), IsotropicMaterial(70.0e9, 0.3, 2700.0), t,
                 StressState::Shell, IntegrationOptions());
  ShellOptions options;
  options.drilling_factor = drilling;
  model.set_shell_options(options);
  model.constraints().push_back(fix(region("root", {box(0.0, 0.0)}), {0, 1, 2, 3, 4, 5}));
  LoadCaseSpec lc;
  lc.name = torsion ? "torque" : "shear";
  const auto wall_tip = [&](const std::string& name, Scalar y0, Scalar y1, Scalar z0, Scalar z1,
                            const Vector3& traction) {
    TractionLoadSpec tr;
    tr.region = region(name, {box(length, length, y0, y1, z0, z1)});
    tr.traction = traction;
    lc.tractions.push_back(tr);
  };
  if (torsion) {
    const Scalar torque = 100.0;
    const Scalar tau = torque / (2.0 * b * b * t);  // shear flow over the wall
    wall_tip("bottom", -kInf, kInf, -h, -h, Vector3(0.0, tau, 0.0));
    wall_tip("right", h, h, -kInf, kInf, Vector3(0.0, 0.0, tau));
    wall_tip("top", -kInf, kInf, h, h, Vector3(0.0, -tau, 0.0));
    wall_tip("left", -h, -h, -kInf, kInf, Vector3(0.0, 0.0, -tau));
  } else {
    const Scalar force = 1000.0;
    const Scalar tau = force / (2.0 * b * t);
    wall_tip("right", h, h, -kInf, kInf, Vector3(0.0, 0.0, -tau));
    wall_tip("left", -h, -h, -kInf, kInf, Vector3(0.0, 0.0, -tau));
  }
  model.load_case_specs() = {lc};
  const StaticSolution sol = solve_static(model);
  BoxResult out;
  Index count = 0;
  for (Index k = 0; k < model.mesh().num_nodes(); ++k) {
    const Vector3 x = model.mesh().node(k);
    if (std::abs(x.x() - length) > 1.0e-12) continue;
    const Vector3 u = sol.displacement.segment<3>(6 * k);
    const Vector3 r(0.0, x.y(), x.z());
    out.deflection += u.z();
    out.twist += r.cross(u).x() / r.squaredNorm();
    ++count;
  }
  out.deflection /= static_cast<Scalar>(count);
  out.twist /= static_cast<Scalar>(count);
  return out;
}

}  // namespace

StudyOutcome study_shell_box_beam(const std::string& out_dir, json::Value& summary) {
  const Scalar b = 0.1;
  const Scalar length = 2.0;
  const Scalar t = 0.002;
  const Scalar e_mod = 70.0e9;
  const Scalar g = e_mod / 2.6;
  // Mid-line section: I = 2 t b^3 / 3, the side walls carry the shear
  // (A_s = 2 b t); Bredt's J = 4 A^2 t / s = b^3 t.
  const Scalar second_moment = 2.0 * t * b * b * b / 3.0;
  const Scalar beam = -1000.0 * length * length * length / (3.0 * e_mod * second_moment) -
                      1000.0 * length / (g * 2.0 * b * t);
  const Scalar bredt = 100.0 * length / (g * b * b * b * t);
  CsvWriter csv(path_join(out_dir, "shell_box_beam.csv"),
                {"cells_per_wall", "drilling_factor", "tip_deflection[m]", "beam_theory[m]",
                 "deflection_ratio[-]", "tip_twist[rad]", "bredt[rad]", "twist_ratio[-]"});
  json::Value rows = json::Value::make_array();
  Scalar bending_gap = 0.0;
  Scalar torsion_gap = 0.0;
  Scalar sensitivity = 0.0;
  BoxResult reference;
  for (const Index nb : {4, 8}) {
    std::vector<Scalar> factors{Shell4Element::kDefaultDrillingFactor};
    if (nb == 8) factors = {1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1};
    for (const Scalar drilling : factors) {
      const BoxResult bend = box_beam(nb, drilling, false);
      const BoxResult twist = box_beam(nb, drilling, true);
      BoxResult r{bend.deflection, twist.twist};
      csv.raw_row({std::to_string(nb), fmt(drilling, 2), fmt(r.deflection, 10), fmt(beam, 10),
                   fmt(r.deflection / beam, 6), fmt(r.twist, 10), fmt(bredt, 10),
                   fmt(r.twist / bredt, 6)});
      json::Value e = json::Value::make_object();
      e.set("cells_per_wall", json::Value::make_number(static_cast<Scalar>(nb)));
      e.set("drilling_factor", json::Value::make_number(drilling));
      e.set("deflection_ratio", json::Value::make_number(r.deflection / beam));
      e.set("twist_ratio", json::Value::make_number(r.twist / bredt));
      rows.push_back(e);
      if (nb == 8 && drilling == Shell4Element::kDefaultDrillingFactor) {
        reference = r;
        bending_gap = std::abs(r.deflection / beam - 1.0);
        torsion_gap = std::abs(r.twist / bredt - 1.0);
      }
    }
  }
  // The drilling factor from 1e-5 to 1e-2 against the default.
  for (const Scalar drilling : {1.0e-5, 1.0e-4, 1.0e-2}) {
    const BoxResult bend = box_beam(8, drilling, false);
    const BoxResult twist = box_beam(8, drilling, true);
    sensitivity = std::max({sensitivity,
                            std::abs(bend.deflection / reference.deflection - 1.0),
                            std::abs(twist.twist / reference.twist - 1.0)});
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("validation (beam theory, Bredt) and the "
                                             "drilling stiffness at folds"));
  block.set("rows", rows);
  block.set("beam_theory_m", json::Value::make_number(beam));
  block.set("bredt_rad", json::Value::make_number(bredt));
  block.set("drilling_sensitivity", json::Value::make_number(sensitivity));
  block.set("note",
            json::Value::make_string(
                "A cantilever of square box section, mid-surface side 0.1 m, walls 2 mm, 2 m "
                "long, E = 70 GPa, nu = 0.3, clamped at x = 0; its four walls meet at 90 "
                "degree folds, where each cell keeps its own normal and the rotations are "
                "coupled by the drilling stiffness. Bending: a 1 kN shear on the side "
                "walls' tip edges against P L^3 / (3 E I) + P L / (G A_s), I = 2 t b^3 / 3 "
                "and A_s = 2 b t. Torsion: 100 N m as Bredt's shear flow round the tip "
                "against T L / (G J), J = b^3 t (a square box does not warp). Both measured "
                "as the tip section's mean; the root clamp and the load's spread differ "
                "from beam theory within about a wall width of the ends. The drilling "
                "factor is swept from 1e-6 to 1e-1 of G t on 8 cells per wall; "
                "drilling_sensitivity is the largest change from 1e-5 to 1e-2 against the "
                "default 1e-3."));
  summary.set("shell_box_beam", block);
  StudyOutcome outcome;
  outcome.name = "box-section cantilever: bending, torsion, drilling stiffness at folds";
  outcome.kind = "validation";
  outcome.metric = "largest gap to beam theory and Bredt at 8 cells per wall";
  outcome.value = std::max(bending_gap, torsion_gap);
  outcome.tolerance = 2.0e-2;
  outcome.passed = outcome.value <= outcome.tolerance && sensitivity <= 1.0e-2;
  std::ostringstream note;
  note << "bending " << fmt(1.0 + (reference.deflection / beam - 1.0), 5) << ", torsion "
       << fmt(reference.twist / bredt, 5) << " of the theory; the drilling factor over "
       << "1e-5 ... 1e-2 changes them by at most " << fmt(sensitivity, 3);
  outcome.note = note.str();
  return outcome;
}

}  // namespace verify
}  // namespace sparlab
