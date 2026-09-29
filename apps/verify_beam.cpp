/// \file verify_beam.cpp
/// \brief Verification of the two-node Timoshenko beam against exact
///        solutions of the beam model it discretises.
///
/// Studies:
///   * `beam-exact`     an inclined cantilever in space under end forces, end
///                      moments and a uniform load along all three axes, and
///                      an L-shaped frame whose tip load bends one arm and
///                      twists the other: the nodal displacements and the end
///                      resultants of every element are the Timoshenko
///                      beam's on every mesh, one element per member
///                      included (the interpolation solves the beam's
///                      homogeneous equations);
///   * `beam-modes`     the twelve lowest natural frequencies of a simply
///                      supported beam - bending in both planes, torsion and
///                      stretching - against the exact Timoshenko
///                      frequencies (shear deformation and rotary inertia),
///                      with consistent and lumped mass;
///   * `beam-harmonic`  the same beam under a uniform load cos(omega t) in
///                      both planes, undamped and with a loss factor, against
///                      the exact series of its modes;
///   * `beam-buckling`  pinned and cantilever columns against the exact
///                      buckling loads of the model (Euler's with the shear
///                      deformation and the rotations' geometric term), and a
///                      section of small torsion constant, which buckles in
///                      torsion at G J A / I_p on every mesh;
///   * `beam-curved`    a quarter-circle cantilever of straight elements under
///                      out-of-plane and in-plane tip loads, against
///                      Castigliano's solution with bending, torsion,
///                      stretching and shear.
///
/// The references are the solutions of the model the element discretises:
/// the Timoshenko beam with Saint-Venant torsion, the rotary inertia
/// rho I and the geometric stiffness
/// N [u'^2 + v'^2 + w'^2 + (I_p/A) theta_x'^2 + (I_y/A) theta_y'^2 +
/// (I_z/A) theta_z'^2] of a uniform axial force (Beam2.hpp). A simply
/// supported beam's modes are v = V sin(kx), theta = Theta cos(kx) in each
/// plane, so each exact frequency or buckling load is the smallest root of a
/// 2 x 2 problem; a cantilever column's are v = V (1 - cos(kx)),
/// theta = Theta sin(kx), the same problem with k = (2j - 1) pi / (2L).
#include "VerifySupport.hpp"

#include "AppSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/elements/Beam2.hpp"
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
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <sstream>
#include <vector>

namespace sparlab {
namespace verify {
namespace {

constexpr Scalar kPi = 3.14159265358979323846;
constexpr Scalar kInf = std::numeric_limits<Scalar>::infinity();

// Steel, and a rectangle 40 mm wide (along y') by 100 mm deep (along z').
constexpr Scalar kE = 210.0e9;
constexpr Scalar kNu = 0.3;
constexpr Scalar kRho = 7850.0;
constexpr Scalar kWidth = 0.04;
constexpr Scalar kHeight = 0.1;

using Vector6 = Eigen::Matrix<Scalar, 6, 1>;

std::string fmt(Scalar v, int digits = 6) { return app::format(v, digits); }

IsotropicMaterial steel() { return IsotropicMaterial(kE, kNu, kRho); }

BeamSection rectangle(const Vector3& orientation = Vector3::Zero()) {
  BeamSection s;
  s.name = "rectangle";
  s.shape = BeamSectionShape::Rectangle;
  s.width = kWidth;
  s.height = kHeight;
  s.orientation = orientation;
  return s;
}

SelectorGroup point_set(const std::string& name) {
  SelectorGroup g;
  g.name = name;
  Selector s;
  s.kind = SelectorKind::Group;
  s.group = name;
  g.members.push_back(s);
  return g;
}

SelectorGroup everything() {
  SelectorGroup g;
  g.name = "all";
  g.members.emplace_back();
  return g;
}

DisplacementConstraint fix(const SelectorGroup& where, std::initializer_list<int> components) {
  DisplacementConstraint c;
  c.region = where;
  for (int k : components) c.set(k, true);
  return c;
}

/// A frame model of one section throughout.
FemModel frame_model(const FrameMeshSpec& spec, const BeamSection& section) {
  Mesh mesh = make_frame_mesh(spec);
  std::vector<Index> all(static_cast<std::size_t>(mesh.num_elements()));
  std::iota(all.begin(), all.end(), Index{0});
  FemModel model(std::move(mesh), steel(), 1.0, StressState::Beam, IntegrationOptions());
  model.assign_section(section, all);
  return model;
}

/// A straight member of n elements from `a` to `b` (points "A" and "B").
FemModel straight_beam(Index n, const Vector3& a, const Vector3& b, const BeamSection& section) {
  FrameMeshSpec spec;
  spec.points = {{"A", a}, {"B", b}};
  FrameMember member;
  member.name = "beam";
  member.from = "A";
  member.to = "B";
  member.elements = n;
  spec.members = {member};
  return frame_model(spec, section);
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

StaticSolution solve_static(FemModel& model) {
  model.finalize();
  Assembler assembler(model);
  StaticAnalysisOptions options;
  options.linear.type = LinearSolverType::SimplicialLdlt;
  StaticAnalysis analysis(model, assembler, options);
  return analysis.solve_all().front();
}

Vector6 node_dofs(const FemModel& model, const Vector& u, Index node) {
  Vector6 out;
  for (int k = 0; k < 6; ++k) out(k) = u(model.dofs().dof(node, k));
  return out;
}

/// Observed order between consecutive meshes of n and 2n elements.
Scalar order(Scalar coarse, Scalar fine) { return observed_order(2.0, coarse, 1.0, fine); }

/// The resolved rectangle and its moduli.
struct Stiffness {
  BeamSection s;
  Scalar g = 0.0;
  Scalar ea = 0.0;
  Scalar gj = 0.0;
  Scalar eiy = 0.0;
  Scalar eiz = 0.0;
  Scalar gay = 0.0;  ///< k_y G A
  Scalar gaz = 0.0;  ///< k_z G A
};

Stiffness stiffness_of(const BeamSection& section) {
  Stiffness k;
  k.s = resolve(section, kNu);
  k.g = kE / (2.0 * (1.0 + kNu));
  k.ea = kE * k.s.area;
  k.gj = k.g * k.s.torsion;
  k.eiy = kE * k.s.iy;
  k.eiz = kE * k.s.iz;
  k.gay = k.s.shear_y * k.g * k.s.area;
  k.gaz = k.s.shear_z * k.g * k.s.area;
  return k;
}

// ---------------------------------------------------------------------------
// Exactness under end loads and uniform loads
// ---------------------------------------------------------------------------

/// The loads of the inclined cantilever, in the member's axes: an end force
/// (N, P_y, P_z), an end moment (T, M_y, M_z) and a uniform load q.
struct CantileverLoads {
  Vector3 force;
  Vector3 moment;
  Vector3 q;
};

/// The Timoshenko cantilever clamped at x = 0, free at x = L: the
/// displacement (u, v, w, theta_x, theta_y, theta_z) at x, local.
Vector6 cantilever_field(const Stiffness& k, const CantileverLoads& p, Scalar length,
                         Scalar x) {
  const Scalar l = length;
  const Scalar n = p.force(0);
  const Scalar py = p.force(1);
  const Scalar pz = p.force(2);
  const Scalar t = p.moment(0);
  const Scalar my = p.moment(1);
  const Scalar mz = p.moment(2);
  const Scalar qx = p.q(0);
  const Scalar qy = p.q(1);
  const Scalar qz = p.q(2);
  const Scalar x2 = x * x;
  const Scalar tip = x2 * (3.0 * l - x) / 6.0;              // int int (L - s)
  const Scalar slope = l * x - 0.5 * x2;                    // int (L - s)
  const Scalar uniform = x2 * (6.0 * l * l - 4.0 * l * x + x2) / 24.0;
  const Scalar uniform_slope = x * (3.0 * l * l - 3.0 * l * x + x2) / 6.0;
  Vector6 f;
  f(0) = (n * x + qx * slope) / k.ea;
  f(1) = py * tip / k.eiz + py * x / k.gay + mz * 0.5 * x2 / k.eiz + qy * uniform / k.eiz +
         qy * slope / k.gay;
  f(2) = pz * tip / k.eiy + pz * x / k.gaz - my * 0.5 * x2 / k.eiy + qz * uniform / k.eiy +
         qz * slope / k.gaz;
  f(3) = t * x / k.gj;
  f(4) = -pz * slope / k.eiy + my * x / k.eiy - qz * uniform_slope / k.eiy;
  f(5) = py * slope / k.eiz + mz * x / k.eiz + qy * uniform_slope / k.eiz;
  return f;
}

/// The section resultants at x of that cantilever (the part beyond x acting
/// on the part before it), local: N, Q_y, Q_z, T, M_y, M_z.
Vector6 cantilever_resultants(const CantileverLoads& p, Scalar length, Scalar x) {
  const Scalar r = length - x;
  const Vector3 axis = Vector3::UnitX();
  Vector6 out;
  out.head<3>() = p.force + r * p.q;
  out.tail<3>() = p.moment + r * axis.cross(p.force) + 0.5 * r * r * axis.cross(p.q);
  return out;
}

struct ExactErrors {
  Scalar translation = 0.0;
  Scalar rotation = 0.0;
  Scalar force = 0.0;
  Scalar moment = 0.0;
  Scalar worst() const { return std::max({translation, rotation, force, moment}); }
};

ExactErrors inclined_cantilever(Index n, const CantileverLoads& loads) {
  const Vector3 a(0.2, -0.1, 0.3);
  const Scalar length = 1.3;
  const Vector3 axis = Vector3(2.0, 1.0, 1.5).normalized();
  const Vector3 b = a + length * axis;
  const BeamSection section = rectangle(Vector3(1.0, -1.0, 0.5));
  const Stiffness k = stiffness_of(section);
  FemModel model = straight_beam(n, a, b, section);
  const BeamFrame frame = Beam2Element::frame(model.element_geometry(0));
  const Matrix3& r = frame.rotation;
  model.constraints().push_back(fix(point_set("A"), {0, 1, 2, 3, 4, 5}));
  LoadCaseSpec lc;
  lc.name = "end and uniform loads";
  PointLoadSpec tip;
  tip.region = point_set("B");
  tip.force = r.transpose() * loads.force;
  tip.moment = r.transpose() * loads.moment;
  lc.point_loads.push_back(tip);
  if (loads.q.norm() > 0.0) {
    LineLoadSpec line;
    line.region = everything();
    line.force_per_length = r.transpose() * loads.q;
    lc.line_loads.push_back(line);
  }
  model.load_case_specs() = {lc};
  const StaticSolution solution = solve_static(model);
  const Mesh& mesh = model.mesh();
  ExactErrors err;
  Scalar u_scale = 0.0;
  Scalar r_scale = 0.0;
  for (Index node = 0; node < mesh.num_nodes(); ++node) {
    const Scalar x = (mesh.node(node) - a).dot(axis);
    const Vector6 exact = cantilever_field(k, loads, length, x);
    u_scale = std::max(u_scale, exact.head<3>().cwiseAbs().maxCoeff());
    r_scale = std::max(r_scale, exact.tail<3>().cwiseAbs().maxCoeff());
  }
  for (Index node = 0; node < mesh.num_nodes(); ++node) {
    const Scalar x = (mesh.node(node) - a).dot(axis);
    const Vector6 exact = cantilever_field(k, loads, length, x);
    const Vector6 h = node_dofs(model, solution.displacement, node);
    err.translation = std::max(
        err.translation, (h.head<3>() - r.transpose() * exact.head<3>()).cwiseAbs().maxCoeff() /
                             u_scale);
    err.rotation = std::max(
        err.rotation, (h.tail<3>() - r.transpose() * exact.tail<3>()).cwiseAbs().maxCoeff() /
                          r_scale);
  }
  Assembler assembler(model);
  const BeamField field = recover_beam_forces(model, assembler, solution.displacement, lc);
  const Vector6 root = cantilever_resultants(loads, length, 0.0);
  const Scalar f_scale = root.head<3>().cwiseAbs().maxCoeff();
  const Scalar m_scale = root.tail<3>().cwiseAbs().maxCoeff();
  for (Index e = 0; e < mesh.num_elements(); ++e) {
    const Index* nodes = mesh.element_nodes(e);
    const BeamEndForces& f = field.element[static_cast<std::size_t>(e)];
    // The element's x' runs from its node 0 to node 1; against the member's
    // axis its resultants change sign.
    const Scalar sign = f.frame.x_axis().dot(axis) > 0.0 ? 1.0 : -1.0;
    const Matrix3 turn = f.frame.rotation * r.transpose();
    for (int end = 0; end < 2; ++end) {
      const Scalar x = (mesh.node(nodes[end]) - a).dot(axis);
      const Vector6 exact = cantilever_resultants(loads, length, x);
      const Vector6& h = end == 0 ? f.start : f.end;
      err.force = std::max(err.force, (h.head<3>() - sign * turn * exact.head<3>())
                                              .cwiseAbs()
                                              .maxCoeff() /
                                          f_scale);
      err.moment = std::max(err.moment, (h.tail<3>() - sign * turn * exact.tail<3>())
                                                .cwiseAbs()
                                                .maxCoeff() /
                                            m_scale);
    }
  }
  return err;
}

/// An L-shaped frame in the x-y plane - an arm a along x from the clamp A,
/// an arm b along y from the corner B to the tip C - under a downward load P
/// at C: the arm AB bends and twists, BC bends. At C,
/// w = -P [(a^3 + b^3) / (3 E I_y) + (a + b) / (k G A) + a b^2 / (G J)],
/// theta_x = -P a b / (G J) - P b^2 / (2 E I_y), theta_y = P a^2 / (2 E I_y);
/// at the clamp the resultants in global axes are (0, 0, -P) and
/// (-P b, P a, 0).
ExactErrors l_frame(Index n, Scalar* tip_deflection) {
  const Scalar a = 1.2;
  const Scalar b = 0.8;
  const Scalar p = 2.0e3;
  const BeamSection section = rectangle();
  const Stiffness k = stiffness_of(section);
  FrameMeshSpec spec;
  spec.points = {{"A", Vector3::Zero()}, {"B", Vector3(a, 0.0, 0.0)}, {"C", Vector3(a, b, 0.0)}};
  FrameMember ab;
  ab.name = "AB";
  ab.from = "A";
  ab.to = "B";
  ab.elements = n;
  FrameMember bc = ab;
  bc.name = "BC";
  bc.from = "B";
  bc.to = "C";
  spec.members = {ab, bc};
  FemModel model = frame_model(spec, section);
  model.constraints().push_back(fix(point_set("A"), {0, 1, 2, 3, 4, 5}));
  LoadCaseSpec lc;
  lc.name = "tip load";
  PointLoadSpec load;
  load.region = point_set("C");
  load.force = Vector3(0.0, 0.0, -p);
  lc.point_loads.push_back(load);
  model.load_case_specs() = {lc};
  const StaticSolution solution = solve_static(model);
  Vector6 exact = Vector6::Zero();
  exact(2) = -p * ((a * a * a + b * b * b) / (3.0 * k.eiy) + (a + b) / k.gaz +
                   a * b * b / k.gj);
  exact(3) = -p * a * b / k.gj - p * b * b / (2.0 * k.eiy);
  exact(4) = p * a * a / (2.0 * k.eiy);
  const Index c = node_at(model.mesh(), Vector3(a, b, 0.0));
  const Vector6 h = node_dofs(model, solution.displacement, c);
  if (tip_deflection != nullptr) *tip_deflection = h(2);
  ExactErrors err;
  err.translation = (h.head<3>() - exact.head<3>()).cwiseAbs().maxCoeff() / std::abs(exact(2));
  err.rotation = (h.tail<3>() - exact.tail<3>()).cwiseAbs().maxCoeff() /
                 exact.tail<3>().cwiseAbs().maxCoeff();
  // The clamp's resultants, in AB's axes (x' = x, y' = y, z' = z).
  Assembler assembler(model);
  const BeamField field = recover_beam_forces(model, assembler, solution.displacement, lc);
  const Index root = node_at(model.mesh(), Vector3::Zero());
  Vector6 expected = Vector6::Zero();
  expected(2) = -p;
  expected(3) = -p * b;
  expected(4) = p * a;
  for (Index e = 0; e < model.mesh().num_elements(); ++e) {
    const Index* nodes = model.mesh().element_nodes(e);
    if (nodes[0] != root && nodes[1] != root) continue;
    const BeamEndForces& f = field.element[static_cast<std::size_t>(e)];
    // The resultant at x' acts from the side of larger x': the frame's on
    // the clamp at a start, the clamp's on the frame at an end.
    const bool start = nodes[0] == root;
    const Vector6& at_root = start ? f.start : f.end;
    const Scalar sign = start ? 1.0 : -1.0;
    const Vector3 force = sign * f.frame.rotation.transpose() * at_root.head<3>();
    const Vector3 moment = sign * f.frame.rotation.transpose() * at_root.tail<3>();
    err.force = std::max(err.force, (force - expected.head<3>()).cwiseAbs().maxCoeff() / p);
    err.moment =
        std::max(err.moment, (moment - expected.tail<3>()).cwiseAbs().maxCoeff() / (p * a));
  }
  return err;
}

}  // namespace

StudyOutcome study_beam_exact(const std::string& out_dir, json::Value& summary) {
  const std::vector<Index> ladder{1, 2, 4, 8, 16};
  struct Case {
    std::string name;
    CantileverLoads loads;
  };
  std::vector<Case> cases;
  {
    Case c;
    c.name = "end forces and moments";
    c.loads.force = Vector3(2.0e5, 1.0e3, 3.0e3);
    c.loads.moment = Vector3(500.0, 800.0, -400.0);
    c.loads.q = Vector3::Zero();
    cases.push_back(c);
    c.name = "uniform load";
    c.loads.force = Vector3::Zero();
    c.loads.moment = Vector3::Zero();
    c.loads.q = Vector3(1.0e5, -2.0e3, 4.0e3);
    cases.push_back(c);
    c.name = "both";
    c.loads.force = Vector3(-1.5e5, -2.0e3, 1.0e3);
    c.loads.moment = Vector3(-300.0, 600.0, 900.0);
    c.loads.q = Vector3(5.0e4, 3.0e3, -1.5e3);
    cases.push_back(c);
  }
  CsvWriter csv(path_join(out_dir, "beam_exact.csv"),
                {"model", "loads", "elements_per_member", "translation_error[-]",
                 "rotation_error[-]", "force_error[-]", "moment_error[-]"});
  Scalar worst = 0.0;
  json::Value rows = json::Value::make_array();
  for (const Case& c : cases) {
    for (const Index n : ladder) {
      const ExactErrors e = inclined_cantilever(n, c.loads);
      csv.raw_row({"inclined cantilever", c.name, std::to_string(n), fmt(e.translation, 3),
                   fmt(e.rotation, 3), fmt(e.force, 3), fmt(e.moment, 3)});
      worst = std::max(worst, e.worst());
      json::Value row = json::Value::make_object();
      row.set("model", json::Value::make_string("inclined cantilever"));
      row.set("loads", json::Value::make_string(c.name));
      row.set("elements", json::Value::make_number(static_cast<Scalar>(n)));
      row.set("worst", json::Value::make_number(e.worst()));
      rows.push_back(row);
    }
  }
  Scalar deflection = 0.0;
  for (const Index n : ladder) {
    const ExactErrors e = l_frame(n, &deflection);
    csv.raw_row({"L-frame", "tip load", std::to_string(n), fmt(e.translation, 3),
                 fmt(e.rotation, 3), fmt(e.force, 3), fmt(e.moment, 3)});
    worst = std::max(worst, e.worst());
    json::Value row = json::Value::make_object();
    row.set("model", json::Value::make_string("L-frame"));
    row.set("loads", json::Value::make_string("tip load"));
    row.set("elements", json::Value::make_number(static_cast<Scalar>(n)));
    row.set("worst", json::Value::make_number(e.worst()));
    rows.push_back(row);
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification (exact Timoshenko beam)"));
  block.set("rows", rows);
  block.set("l_frame_tip_deflection_m", json::Value::make_number(deflection));
  block.set("note",
            json::Value::make_string(
                "Steel (E = 210 GPa, nu = 0.3) rectangle 40 x 100 mm (Cowper's k), 1 to 16 "
                "elements per member. An inclined cantilever 1.3 m long, its y' axis set by an "
                "orientation vector, under an end force (N, P_y, P_z), an end moment "
                "(T, M_y, M_z) and a uniform load along all three axes, against the closed-form "
                "Timoshenko field at every node (bending with shear in both planes, torsion, "
                "stretching) and the statics at both ends of every element; and an L-frame "
                "(arms 1.2 m along x and 0.8 m along y) under a downward tip load, against "
                "Castigliano's tip displacement and rotations (bending of both arms, their "
                "shear, the twist of the first) and the clamp's resultants. Errors: the "
                "largest nodal difference over the largest exact value, translations and "
                "rotations apart; the end forces and moments over the root's."));
  summary.set("beam_exact", block);
  StudyOutcome outcome;
  outcome.name = "beam exactness: inclined cantilever and L-frame, 1 to 16 elements";
  outcome.kind = "verification";
  outcome.metric = "largest relative error of displacements, rotations and end resultants";
  outcome.value = worst;
  outcome.tolerance = 1.0e-10;
  outcome.passed = worst <= outcome.tolerance;
  outcome.note = "exact on every mesh: nodal loads and uniform loads alike";
  return outcome;
}

namespace {

// ---------------------------------------------------------------------------
// The simply supported beam: exact modes, harmonic response, buckling
// ---------------------------------------------------------------------------

/// One bending plane of the Timoshenko beam in the mode
/// v = V sin(kx), theta = Theta cos(kx) (per unit length, the factor 1/2 of
/// the averages dropped): stiffness, mass, and the geometric stiffness of a
/// unit compressive force.
struct PlaneMode {
  Eigen::Matrix2d k;
  Eigen::Matrix2d m;
  Eigen::Matrix2d g;
};

PlaneMode plane_mode(Scalar wavenumber, Scalar ei, Scalar kga, Scalar area, Scalar inertia) {
  const Scalar q = wavenumber;
  PlaneMode out;
  out.k << kga * q * q, -kga * q, -kga * q, ei * q * q + kga;
  out.m << kRho * area, 0.0, 0.0, kRho * inertia;
  out.g << q * q, 0.0, 0.0, inertia / area * q * q;
  return out;
}

Scalar smallest(const Eigen::Matrix2d& k, const Eigen::Matrix2d& b) {
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::Matrix2d> ges(k, b);
  return ges.eigenvalues()(0);
}

/// An exact frequency or load with the kind of mode it belongs to.
struct Exact {
  std::string kind;
  Scalar value;
};

/// A simply supported beam along x, length L, of n elements: v = w = 0 at
/// both ends, u and the twist held at x = 0.
FemModel simply_supported(Index n, Scalar length) {
  FemModel model = straight_beam(n, Vector3::Zero(), Vector3(length, 0.0, 0.0), rectangle());
  model.constraints().push_back(fix(point_set("A"), {0, 1, 2, 3}));
  model.constraints().push_back(fix(point_set("B"), {1, 2}));
  return model;
}

/// The exact natural frequencies [rad/s] of that beam below `limit`: bending
/// in the x-y plane (I_z, k_y) and the x-z plane (I_y, k_z) with n
/// half-waves, and the fixed-free torsion and stretching.
std::vector<Exact> exact_frequencies(Scalar length, Scalar limit) {
  const Stiffness k = stiffness_of(rectangle());
  std::vector<Exact> out;
  for (int n = 1; n <= 40; ++n) {
    const Scalar q = n * kPi / length;
    const PlaneMode y = plane_mode(q, k.eiz, k.gay, k.s.area, k.s.iz);
    const PlaneMode z = plane_mode(q, k.eiy, k.gaz, k.s.area, k.s.iy);
    out.push_back({"bending y' " + std::to_string(n), std::sqrt(smallest(y.k, y.m))});
    out.push_back({"bending z' " + std::to_string(n), std::sqrt(smallest(z.k, z.m))});
    const Scalar fixed_free = (2 * n - 1) * kPi / (2.0 * length);
    out.push_back({"torsion " + std::to_string(n),
                   fixed_free * std::sqrt(k.gj / (kRho * k.s.polar()))});
    out.push_back({"axial " + std::to_string(n), fixed_free * std::sqrt(kE / kRho)});
  }
  std::sort(out.begin(), out.end(), [](const Exact& a, const Exact& b) { return a.value < b.value; });
  out.erase(std::remove_if(out.begin(), out.end(), [&](const Exact& e) { return e.value > limit; }),
            out.end());
  return out;
}

}  // namespace

StudyOutcome study_beam_modes(const std::string& out_dir, json::Value& summary) {
  const Scalar length = 1.0;
  const int modes = 12;
  const std::vector<Exact> exact = exact_frequencies(length, kInf);
  const std::vector<Index> ladder{16, 32, 64, 128};
  CsvWriter csv(path_join(out_dir, "beam_modes.csv"),
                {"mass", "n", "mode", "kind", "frequency[Hz]", "exact[Hz]", "relative_error[-]",
                 "observed_order"});
  Scalar worst_consistent = 0.0;
  Scalar worst_lumped = 0.0;
  Scalar worst_order = kInf;
  // The smallest relative gap between neighbouring exact frequencies: the
  // modes pair up by order while the errors stay below half of it.
  Scalar gap = kInf;
  for (int i = 0; i + 1 <= modes; ++i) {
    gap = std::min(gap, (exact[static_cast<std::size_t>(i) + 1].value -
                         exact[static_cast<std::size_t>(i)].value) /
                            exact[static_cast<std::size_t>(i)].value);
  }
  Scalar largest = 0.0;
  json::Value meshes = json::Value::make_array();
  for (const MassType mass : {MassType::Consistent, MassType::Lumped}) {
    const std::string mass_name = mass == MassType::Consistent ? "consistent" : "lumped";
    std::vector<Scalar> previous(modes, 0.0);
    for (const Index n : ladder) {
      FemModel model = simply_supported(n, length);
      LoadCaseSpec lc;
      lc.name = "none";
      lc.prescribed_displacement_only = true;
      model.load_case_specs() = {lc};
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
        const Exact& x = exact[static_cast<std::size_t>(i)];
        const Scalar error = std::abs(r.angular_frequencies(i) - x.value) / x.value;
        const Scalar p = previous[static_cast<std::size_t>(i)] > 0.0
                             ? order(previous[static_cast<std::size_t>(i)], error)
                             : 0.0;
        csv.raw_row({mass_name, std::to_string(n), std::to_string(i + 1), x.kind,
                     fmt(r.frequencies_hz(i), 10), fmt(x.value / (2.0 * kPi), 10),
                     fmt(error, 4), previous[static_cast<std::size_t>(i)] > 0.0 ? fmt(p, 4) : ""});
        errors.push_back(json::Value::make_number(error));
        largest = std::max(largest, error);
        if (n == ladder.back()) {
          (mass == MassType::Consistent ? worst_consistent : worst_lumped) =
              std::max(mass == MassType::Consistent ? worst_consistent : worst_lumped, error);
          worst_order = std::min(worst_order, p);
        }
        previous[static_cast<std::size_t>(i)] = error;
      }
      entry.set("errors", errors);
      meshes.push_back(entry);
    }
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification (exact Timoshenko frequencies)"));
  block.set("meshes", meshes);
  json::Value reference = json::Value::make_array();
  for (int i = 0; i < modes; ++i) {
    json::Value e = json::Value::make_object();
    e.set("kind", json::Value::make_string(exact[static_cast<std::size_t>(i)].kind));
    e.set("hz", json::Value::make_number(exact[static_cast<std::size_t>(i)].value / (2.0 * kPi)));
    reference.push_back(e);
  }
  block.set("exact", reference);
  block.set("smallest_relative_gap", json::Value::make_number(gap));
  block.set("note",
            json::Value::make_string(
                "Steel beam 1 m long, rectangle 40 x 100 mm, simply supported (v = w = 0 at "
                "both ends; u and the twist held at x = 0): the twelve lowest frequencies - "
                "bending in the x-y plane (I_z) and the x-z plane (I_y), torsion and "
                "stretching - against the exact frequencies of the Timoshenko beam with "
                "Cowper's shear coefficient and the rotary inertia rho I (the smallest root of "
                "the 2 x 2 problem of each mode v = V sin(n pi x / L), theta = Theta cos) and "
                "of Saint-Venant torsion with the polar inertia rho I_p. Consistent mass, and "
                "the lumped one: rho A L / 2 on each node's translations and half the "
                "sections' inertia tensor, rho L / 2 diag(I_p, I_y, I_z) turned to global "
                "axes, on its rotations. Both converge at O(h^2): the twist and the stretch "
                "are linear, and with elements shorter than the depth (Phi = 12 E I / "
                "(k G A L^2) from 1.3 to 500 here) the interpolated rotation is nearly linear "
                "too - the consistent mass above the exact frequencies, the lumped one "
                "below."));
  block.set("largest_error_consistent", json::Value::make_number(worst_consistent));
  block.set("largest_error_lumped", json::Value::make_number(worst_lumped));
  summary.set("beam_modes", block);
  StudyOutcome outcome;
  outcome.name = "beam: twelve lowest natural frequencies, consistent and lumped mass";
  outcome.kind = "verification";
  outcome.metric = "largest frequency error at 128 elements, consistent mass";
  outcome.value = worst_consistent;
  // The largest error is the third torsional mode's, of the linear twist.
  outcome.tolerance = 2.0e-4;
  const Scalar lumped_tolerance = 2.0e-4;
  outcome.passed = worst_consistent <= outcome.tolerance && worst_lumped <= lumped_tolerance &&
                   worst_order >= 1.9 && largest < 0.5 * gap;
  std::ostringstream note;
  note << "lumped mass " << fmt(worst_lumped, 3) << " (tolerance " << fmt(lumped_tolerance, 3)
       << "); smallest observed order " << fmt(worst_order, 3)
       << " (64 -> 128, both masses); every error below half the smallest gap between exact "
          "frequencies, "
       << fmt(gap, 3);
  outcome.note = note.str();
  return outcome;
}

namespace {

/// The exact steady midspan amplitude of the simply supported beam under a
/// uniform load q cos(omega t) in one plane, loss factor eta: the series over
/// the odd modes, 4 q / (n pi) each, of (k (1 + i eta) - omega^2 m) x =
/// (4 q / (n pi), 0), times sin(n pi / 2).
std::complex<Scalar> harmonic_midspan(Scalar length, Scalar ei, Scalar kga, Scalar area,
                                      Scalar inertia, Scalar q, Scalar omega, Scalar eta) {
  using Complex = std::complex<Scalar>;
  Complex sum = 0.0;
  for (int n = 1; n <= 40001; n += 2) {
    const PlaneMode mode = plane_mode(n * kPi / length, ei, kga, area, inertia);
    const Eigen::Matrix2cd dynamic =
        mode.k.cast<Complex>() * Complex(1.0, eta) - (omega * omega) * mode.m.cast<Complex>();
    Eigen::Vector2cd load(4.0 * q / (n * kPi), 0.0);
    const Eigen::Vector2cd x = dynamic.fullPivLu().solve(load);
    sum += (((n - 1) / 2) % 2 == 0 ? 1.0 : -1.0) * x(0);
  }
  return sum;
}

}  // namespace

StudyOutcome study_beam_harmonic(const std::string& out_dir, json::Value& summary) {
  const Scalar length = 1.0;
  const Vector3 q(0.0, 1.0e3, 2.0e3);  // [N/m] along y' = y and z' = z
  const Stiffness k = stiffness_of(rectangle());
  // Static; below the lowest resonance the load excites (bending y' 1 at
  // 93.6 Hz); between bending z' 1 (230.7 Hz) and y' 3 (824.9 Hz); between
  // y' 3 and z' 3 (1860 Hz).
  const std::vector<Scalar> frequencies{0.0, 60.0, 400.0, 1000.0};
  const std::vector<Scalar> losses{0.0, 0.05};
  const std::vector<Index> ladder{16, 32, 64, 128};
  // The series against the closed-form static midspan deflection,
  // 5 q L^4 / (384 E I) + q L^2 / (8 k G A).
  Scalar series_check = 0.0;
  for (int plane = 1; plane <= 2; ++plane) {
    const Scalar ei = plane == 1 ? k.eiz : k.eiy;
    const Scalar kga = plane == 1 ? k.gay : k.gaz;
    const Scalar inertia = plane == 1 ? k.s.iz : k.s.iy;
    const Scalar closed = 5.0 * q(plane) * std::pow(length, 4) / (384.0 * ei) +
                          q(plane) * length * length / (8.0 * kga);
    const Scalar series =
        harmonic_midspan(length, ei, kga, k.s.area, inertia, q(plane), 0.0, 0.0).real();
    series_check = std::max(series_check, std::abs(series - closed) / closed);
  }
  CsvWriter csv(path_join(out_dir, "beam_harmonic.csv"),
                {"mass", "loss_factor", "frequency[Hz]", "n", "direction", "real[m]", "imag[m]",
                 "exact_real[m]", "exact_imag[m]", "relative_error[-]", "observed_order"});
  Scalar worst_fine = 0.0;
  Scalar worst_order = kInf;
  json::Value cases = json::Value::make_array();
  for (const MassType mass : {MassType::Consistent, MassType::Lumped}) {
    const std::string mass_name = mass == MassType::Consistent ? "consistent" : "lumped";
    for (const Scalar eta : losses) {
      // exact[j][plane - 1]
      std::vector<std::array<std::complex<Scalar>, 2>> exact;
      for (const Scalar f : frequencies) {
        std::array<std::complex<Scalar>, 2> x;
        for (int plane = 1; plane <= 2; ++plane) {
          x[static_cast<std::size_t>(plane - 1)] = harmonic_midspan(
              length, plane == 1 ? k.eiz : k.eiy, plane == 1 ? k.gay : k.gaz, k.s.area,
              plane == 1 ? k.s.iz : k.s.iy, q(plane), 2.0 * kPi * f, eta);
        }
        exact.push_back(x);
      }
      std::vector<Scalar> previous(2 * frequencies.size(), 0.0);
      json::Value c = json::Value::make_object();
      c.set("mass", json::Value::make_string(mass_name));
      c.set("loss_factor", json::Value::make_number(eta));
      json::Value meshes = json::Value::make_array();
      for (const Index n : ladder) {
        FemModel model = simply_supported(n, length);
        LoadCaseSpec lc;
        lc.name = "uniform";
        LineLoadSpec line;
        line.region = everything();
        line.force_per_length = q;
        lc.line_loads.push_back(line);
        model.load_case_specs() = {lc};
        model.finalize();
        Assembler assembler(model);
        FrequencyResponseOptions options;
        options.frequencies = frequencies;
        options.mass_type = mass;
        options.structural_damping = eta;
        options.snapshot_frequencies = frequencies;
        const FrequencyResponseResult r = solve_frequency_response(model, assembler, 0, options);
        if (r.snapshots.size() != frequencies.size()) {
          throw SolverError("beam-harmonic: expected one snapshot per frequency");
        }
        const Index mid = node_at(model.mesh(), Vector3(0.5 * length, 0.0, 0.0));
        json::Value errors = json::Value::make_array();
        for (std::size_t j = 0; j < frequencies.size(); ++j) {
          for (int plane = 1; plane <= 2; ++plane) {
            const std::complex<Scalar> h =
                r.snapshots[j].displacement(model.dofs().dof(mid, plane));
            const std::complex<Scalar> x = exact[j][static_cast<std::size_t>(plane - 1)];
            const Scalar error = std::abs(h - x) / std::abs(x);
            const std::size_t slot = 2 * j + static_cast<std::size_t>(plane - 1);
            const Scalar p = previous[slot] > 0.0 ? order(previous[slot], error) : 0.0;
            csv.raw_row({mass_name, fmt(eta, 3), fmt(frequencies[j], 6), std::to_string(n),
                         plane == 1 ? "y" : "z", fmt(h.real(), 10), fmt(h.imag(), 10),
                         fmt(x.real(), 10), fmt(x.imag(), 10), fmt(error, 4),
                         previous[slot] > 0.0 ? fmt(p, 4) : ""});
            errors.push_back(json::Value::make_number(error));
            if (n == ladder.back()) {
              worst_fine = std::max(worst_fine, error);
              // At 0 Hz the nodal values are exact: no order to observe.
              if (frequencies[j] > 0.0) worst_order = std::min(worst_order, p);
            }
            previous[slot] = error;
          }
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
  block.set("kind", json::Value::make_string("verification (exact Timoshenko harmonic "
                                             "response)"));
  block.set("cases", cases);
  block.set("series_vs_closed_form_static", json::Value::make_number(series_check));
  block.set("note",
            json::Value::make_string(
                "The simply supported beam of the frequency study under a uniform load "
                "(0, 1, 2) kN/m cos(omega t), solved directly at 0, 60, 400 and 1000 Hz - "
                "static, below the lowest resonance the load excites (y' 1 at 93.6 Hz), "
                "between z' 1 (230.7 Hz) and y' 3 (824.9 Hz), and between y' 3 and z' 3 "
                "(1860 Hz) - undamped and with the loss factor 0.05, on consistent and lumped "
                "mass. The exact midspan amplitude in each plane is the series over the odd "
                "modes (n <= 40001) of the 2 x 2 problems (k (1 + i eta) - omega^2 m) x = "
                "(4 q / (n pi), 0); series_vs_closed_form_static checks it at 0 Hz against "
                "5 q L^4 / (384 E I) + q L^2 / (8 k G A). Errors: |v - v_exact| / |v_exact| "
                "of the complex amplitude."));
  summary.set("beam_harmonic", block);
  StudyOutcome outcome;
  outcome.name = "beam: harmonic response vs the exact Timoshenko series";
  outcome.kind = "verification";
  outcome.metric = "largest midspan-amplitude error at 128 elements";
  outcome.value = worst_fine;
  // The largest error is at 1000 Hz in the x-y plane, where the modes
  // nearly cancel at midspan: the amplitude is 1/1700 of the static one.
  outcome.tolerance = 1.0e-3;
  outcome.passed =
      worst_fine <= outcome.tolerance && worst_order >= 1.9 && series_check <= 1.0e-9;
  std::ostringstream note;
  note << "smallest observed order " << fmt(worst_order, 3)
       << " (64 -> 128; both masses, eta = 0 and 0.05, 60 ... 1000 Hz); series vs closed form "
          "at 0 Hz "
       << fmt(series_check, 2);
  outcome.note = note.str();
  return outcome;
}

StudyOutcome study_beam_buckling(const std::string& out_dir, json::Value& summary) {
  const Scalar length = 1.0;
  const Scalar p0 = 1.0e3;  // applied compressive force [N]
  const int modes = 4;
  const Stiffness k = stiffness_of(rectangle());
  const std::vector<Index> ladder{4, 8, 16, 32};
  CsvWriter csv(path_join(out_dir, "beam_buckling.csv"),
                {"support", "n", "mode", "kind", "critical_force[N]", "exact[N]",
                 "relative_error[-]", "observed_order"});
  json::Value cases = json::Value::make_array();
  Scalar worst_fine = 0.0;
  Scalar worst_order = kInf;
  for (const bool cantilever : {false, true}) {
    const std::string support = cantilever ? "cantilever" : "pinned";
    // The exact loads: k = n pi / L pinned, (2n - 1) pi / (2L) fixed-free.
    std::vector<Exact> exact;
    for (int n = 1; n <= 10; ++n) {
      const Scalar q = (cantilever ? (2 * n - 1) * kPi / 2.0 : n * kPi) / length;
      const PlaneMode y = plane_mode(q, k.eiz, k.gay, k.s.area, k.s.iz);
      const PlaneMode z = plane_mode(q, k.eiy, k.gaz, k.s.area, k.s.iy);
      exact.push_back({"bending y' " + std::to_string(n), smallest(y.k, y.g)});
      exact.push_back({"bending z' " + std::to_string(n), smallest(z.k, z.g)});
    }
    std::sort(exact.begin(), exact.end(),
              [](const Exact& a, const Exact& b) { return a.value < b.value; });
    json::Value c = json::Value::make_object();
    c.set("support", json::Value::make_string(support));
    const Scalar euler = kPi * kPi * k.eiz / (cantilever ? 4.0 : 1.0) / (length * length);
    c.set("euler_N", json::Value::make_number(euler));
    c.set("exact_N", json::Value::make_number(exact.front().value));
    std::vector<Scalar> previous(modes, 0.0);
    json::Value meshes = json::Value::make_array();
    for (const Index n : ladder) {
      FemModel model = straight_beam(n, Vector3::Zero(), Vector3(length, 0.0, 0.0), rectangle());
      if (cantilever) {
        model.constraints().push_back(fix(point_set("A"), {0, 1, 2, 3, 4, 5}));
      } else {
        model.constraints().push_back(fix(point_set("A"), {0, 1, 2, 3}));
        model.constraints().push_back(fix(point_set("B"), {1, 2}));
      }
      LoadCaseSpec lc;
      lc.name = "compression";
      PointLoadSpec load;
      load.region = point_set("B");
      load.force = Vector3(-p0, 0.0, 0.0);
      lc.point_loads.push_back(load);
      model.load_case_specs() = {lc};
      model.finalize();
      Assembler assembler(model);
      BucklingOptions options;
      options.num_modes = modes;
      const BucklingResult r = analyse_buckling(model, assembler, 0, options);
      json::Value errors = json::Value::make_array();
      for (int i = 0; i < modes; ++i) {
        const Exact& x = exact[static_cast<std::size_t>(i)];
        const Scalar critical = r.load_factors(i) * p0;
        const Scalar error = std::abs(critical - x.value) / x.value;
        const Scalar p = previous[static_cast<std::size_t>(i)] > 0.0
                             ? order(previous[static_cast<std::size_t>(i)], error)
                             : 0.0;
        csv.raw_row({support, std::to_string(n), std::to_string(i + 1), x.kind,
                     fmt(critical, 10), fmt(x.value, 10), fmt(error, 4),
                     previous[static_cast<std::size_t>(i)] > 0.0 ? fmt(p, 4) : ""});
        errors.push_back(json::Value::make_number(error));
        if (n == ladder.back()) {
          worst_fine = std::max(worst_fine, error);
          worst_order = std::min(worst_order, p);
        }
        previous[static_cast<std::size_t>(i)] = error;
      }
      json::Value entry = json::Value::make_object();
      entry.set("n", json::Value::make_number(static_cast<Scalar>(n)));
      entry.set("errors", errors);
      meshes.push_back(entry);
    }
    c.set("meshes", meshes);
    cases.push_back(c);
  }
  // Torsional buckling: with Saint-Venant torsion alone the twist's
  // geometric stiffness N (I_p / A) theta_x'^2 is proportional to its
  // stiffness G J theta_x'^2, so every twist buckles at P = G J A / I_p - on
  // every mesh. A section of small torsion constant buckles so first.
  BeamSection open;
  open.name = "open";
  open.area = k.s.area;
  open.iy = k.s.iy;
  open.iz = k.s.iz;
  open.torsion = 1.0e-9;
  open.shear_y = k.s.shear_y;
  open.shear_z = k.s.shear_z;
  const Scalar torsional = k.g * open.torsion * open.area / open.polar();
  Scalar torsion_error = 0.0;
  for (const Index n : ladder) {
    FemModel model = straight_beam(n, Vector3::Zero(), Vector3(length, 0.0, 0.0), open);
    model.constraints().push_back(fix(point_set("A"), {0, 1, 2, 3}));
    model.constraints().push_back(fix(point_set("B"), {1, 2}));
    LoadCaseSpec lc;
    lc.name = "compression";
    PointLoadSpec load;
    load.region = point_set("B");
    load.force = Vector3(-p0, 0.0, 0.0);
    lc.point_loads.push_back(load);
    model.load_case_specs() = {lc};
    model.finalize();
    Assembler assembler(model);
    BucklingOptions options;
    options.num_modes = 2;
    const BucklingResult r = analyse_buckling(model, assembler, 0, options);
    for (int i = 0; i < 2; ++i) {
      const Scalar error = std::abs(r.load_factors(i) * p0 - torsional) / torsional;
      torsion_error = std::max(torsion_error, error);
      csv.raw_row({"pinned, J = 1e-9 m^4", std::to_string(n), std::to_string(i + 1), "torsion",
                   fmt(r.load_factors(i) * p0, 10), fmt(torsional, 10), fmt(error, 4), ""});
    }
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification (exact buckling loads of the "
                                             "model)"));
  block.set("cases", cases);
  block.set("torsional_N", json::Value::make_number(torsional));
  block.set("torsional_error", json::Value::make_number(torsion_error));
  block.set("note",
            json::Value::make_string(
                "Steel column 1 m long, rectangle 40 x 100 mm, compressed by an end force: "
                "pinned (v = w = 0 at both ends, u and the twist held at the foot) and a "
                "cantilever (clamped foot, free top). Exact: the smallest roots over the "
                "modes of each plane of the 2 x 2 problem of the Timoshenko beam with the "
                "geometric stiffness N [v'^2 + (I/A) theta'^2], the four lowest compared; "
                "they lie below Euler's by the shear deformation and the rotations' term. "
                "The observed order falls from 3 towards 2 as the elements grow shorter than "
                "the depth and the interpolated rotation nearly linear. And the same pinned "
                "column with the torsion constant 1e-9 m^4, which buckles in twist at "
                "G J A / I_p (Saint-Venant torsion, no warping stiffness) with every mode at "
                "once."));
  summary.set("beam_buckling", block);
  StudyOutcome outcome;
  outcome.name = "beam buckling: pinned and cantilever columns, torsional buckling";
  outcome.kind = "verification";
  outcome.metric = "largest critical-load error at 32 elements";
  outcome.value = worst_fine;
  outcome.tolerance = 3.0e-4;
  outcome.passed = worst_fine <= outcome.tolerance && worst_order >= 1.9 &&
                   torsion_error <= 1.0e-10;
  std::ostringstream note;
  note << "smallest observed order " << fmt(worst_order, 3)
       << " (16 -> 32); torsional buckling G J A / I_p within " << fmt(torsion_error, 2)
       << " on every mesh";
  outcome.note = note.str();
  return outcome;
}

StudyOutcome study_beam_curved(const std::string& out_dir, json::Value& summary) {
  const Scalar radius = 1.0;
  const Scalar p = 1.0e3;
  const Stiffness k = stiffness_of(rectangle());
  // Castigliano on the quarter circle theta in [0, pi / 2], clamped at
  // theta = 0, the tip at theta = pi / 2.
  const Scalar r3 = radius * radius * radius;
  const Scalar out_of_plane = p * r3 * kPi / (4.0 * k.eiy) +
                              p * r3 * (0.75 * kPi - 2.0) / k.gj +
                              p * radius * kPi / (2.0 * k.gaz);
  const Scalar in_plane = p * r3 * (0.75 * kPi - 2.0) / k.eiz + p * radius * kPi / (4.0 * k.ea) +
                          p * radius * kPi / (4.0 * k.gay);
  const Scalar in_plane_y =
      -p * r3 / (2.0 * k.eiz) + p * radius / (2.0 * k.ea) - p * radius / (2.0 * k.gay);
  const std::vector<Index> ladder{4, 8, 16, 32, 64, 128};
  CsvWriter csv(path_join(out_dir, "beam_curved.csv"),
                {"load", "n", "component", "displacement[m]", "exact[m]", "relative_error[-]",
                 "observed_order"});
  Scalar worst_fine = 0.0;
  Scalar worst_order = kInf;
  std::array<Scalar, 3> previous{0.0, 0.0, 0.0};
  json::Value rows = json::Value::make_array();
  for (const Index n : ladder) {
    FrameMeshSpec spec;
    spec.points = {{"A", Vector3(radius, 0.0, 0.0)}, {"B", Vector3(0.0, radius, 0.0)}};
    FrameMember arc;
    arc.name = "arc";
    arc.from = "A";
    arc.to = "B";
    arc.elements = n;
    arc.arc = true;
    arc.arc_centre = Vector3::Zero();
    arc.arc_axis = Vector3::UnitZ();
    spec.members = {arc};
    FemModel model = frame_model(spec, rectangle());
    model.constraints().push_back(fix(point_set("A"), {0, 1, 2, 3, 4, 5}));
    LoadCaseSpec down;
    down.name = "out of plane";
    PointLoadSpec fz;
    fz.region = point_set("B");
    fz.force = Vector3(0.0, 0.0, -p);
    down.point_loads.push_back(fz);
    LoadCaseSpec in;
    in.name = "in plane";
    PointLoadSpec fx;
    fx.region = point_set("B");
    fx.force = Vector3(-p, 0.0, 0.0);
    in.point_loads.push_back(fx);
    model.load_case_specs() = {down, in};
    model.finalize();
    Assembler assembler(model);
    StaticAnalysisOptions options;
    options.linear.type = LinearSolverType::SimplicialLdlt;
    StaticAnalysis analysis(model, assembler, options);
    const std::vector<StaticSolution> solutions = analysis.solve_all();
    const Index tip = node_at(model.mesh(), Vector3(0.0, radius, 0.0));
    const std::array<Scalar, 3> h{
        -solutions[0].displacement(model.dofs().dof(tip, 2)),
        -solutions[1].displacement(model.dofs().dof(tip, 0)),
        solutions[1].displacement(model.dofs().dof(tip, 1))};
    const std::array<Scalar, 3> x{out_of_plane, in_plane, in_plane_y};
    const std::array<const char*, 3> load{"out of plane (-z)", "in plane (-x)", "in plane (-x)"};
    const std::array<const char*, 3> component{"-u_z", "-u_x", "u_y"};
    json::Value row = json::Value::make_object();
    row.set("n", json::Value::make_number(static_cast<Scalar>(n)));
    json::Value errors = json::Value::make_array();
    for (std::size_t i = 0; i < 3; ++i) {
      const Scalar error = std::abs(h[i] - x[i]) / std::abs(x[i]);
      const Scalar q = previous[i] > 0.0 ? order(previous[i], error) : 0.0;
      csv.raw_row({load[i], std::to_string(n), component[i], fmt(h[i], 10), fmt(x[i], 10),
                   fmt(error, 4), previous[i] > 0.0 ? fmt(q, 4) : ""});
      errors.push_back(json::Value::make_number(error));
      if (n == ladder.back()) {
        worst_fine = std::max(worst_fine, error);
        worst_order = std::min(worst_order, q);
      }
      previous[i] = error;
    }
    row.set("errors", errors);
    rows.push_back(row);
  }
  csv.close();
  json::Value block = json::Value::make_object();
  block.set("kind", json::Value::make_string("verification (Castigliano, curved beam)"));
  block.set("rows", rows);
  block.set("exact_out_of_plane_m", json::Value::make_number(out_of_plane));
  block.set("exact_in_plane_m", json::Value::make_number(in_plane));
  block.set("exact_in_plane_y_m", json::Value::make_number(in_plane_y));
  block.set("note",
            json::Value::make_string(
                "A quarter circle of radius 1 m in the x-y plane, clamped at (1, 0, 0), of n "
                "straight elements with their nodes on the circle; steel rectangle 40 x 100 mm "
                "(z' out of the plane). A tip load of 1 kN out of the plane (-z) bends it "
                "about the radius and twists it: Castigliano gives "
                "P R^3 pi / (4 E I_y) + P R^3 (3 pi / 4 - 2) / (G J) + P R pi / (2 k G A). One "
                "in the plane (-x) bends, stretches and shears it: "
                "P R^3 (3 pi / 4 - 2) / (E I_z) + P R pi / (4 E A) + P R pi / (4 k G A) along "
                "the load, and -P R^3 / (2 E I_z) + P R / (2 E A) - P R / (2 k G A) along y. "
                "The polygon of chords converges to the circle at O(h^2)."));
  summary.set("beam_curved", block);
  StudyOutcome outcome;
  outcome.name = "beam: quarter-circle cantilever vs Castigliano";
  outcome.kind = "verification";
  outcome.metric = "largest tip-displacement error at 128 elements";
  outcome.value = worst_fine;
  outcome.tolerance = 1.0e-4;
  outcome.passed = worst_fine <= outcome.tolerance && worst_order >= 1.9;
  std::ostringstream note;
  note << "smallest observed order " << fmt(worst_order, 3) << " (64 -> 128)";
  outcome.note = note.str();
  return outcome;
}

}  // namespace verify
}  // namespace sparlab
