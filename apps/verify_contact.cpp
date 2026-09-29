/// \file verify_contact.cpp
/// \brief Verification of unilateral contact against exact solutions of the
///        discrete problem and against Hertz's theory.
///
/// Studies:
///   * `contact-patch`   homogeneous states the linear elements reproduce on
///                       distorted meshes, so the contact solution must too,
///                       to round-off: a block pressed onto a rigid plane
///                       (with and without an initial gap) and two blocks
///                       with non-matching meshes pressed together (the
///                       contact patch test of the dual mortar method), on
///                       Q4, Tri3, Hex8 and Tet4; and a block dragged along
///                       a plane in full slip, whose tangential force is mu
///                       times its normal force exactly in 2-D;
///   * `hertz-line`      Hertz's line contact in plane strain: an elastic
///                       cylinder on a rigid flat, a rigid cylinder pressed
///                       into an elastic block, and an elastic cylinder on a
///                       block of the same material (a mortar pair with
///                       non-matching meshes), on graded Q4 meshes; the
///                       pressure against p0 sqrt(1 - x^2 / a^2), with a and
///                       p0 from the total load, as the mesh is refined -
///                       and two series that identify the floor the
///                       refinement reaches as the finite body and, for the
///                       curved rigid obstacle, a / R;
///   * `hertz-point`     Hertz's point contact: an elastic sphere on a rigid
///                       flat and on an elastic block (mortar, non-matching),
///                       quarter models on graded Hex8 meshes, the pressure
///                       against p0 sqrt(1 - r^2 / a^2).
///
/// Hertz's theory takes each body as a half-space with a parabolic profile
/// and the load normal to the contact plane; the models are finite bodies,
/// so the comparison measures the discretisation error down to the gap
/// between the model and the theory, which the refinement exposes as a
/// floor.

#include "VerifySupport.hpp"

#include "AppSupport.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/core/Logging.hpp"
#include "sparlab/fem/Assembler.hpp"
#include "sparlab/fem/Contact.hpp"
#include "sparlab/fem/NonlinearStatic.hpp"
#include "sparlab/io/CsvWriter.hpp"
#include "sparlab/mesh/StructuredMesh.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <sstream>
#include <vector>

namespace sparlab {
namespace verify {
namespace {

constexpr Scalar kPi = 3.14159265358979323846;
constexpr Scalar kInf = std::numeric_limits<Scalar>::infinity();

std::string fmt(Scalar v, int digits = 6) {
  return std::isnan(v) ? std::string() : app::format(v, digits);
}

SelectorGroup box(Scalar xmin, Scalar xmax, Scalar ymin = -kInf, Scalar ymax = kInf,
                  Scalar zmin = -kInf, Scalar zmax = kInf) {
  SelectorGroup g;
  Selector s;
  s.kind = SelectorKind::Box;
  s.xmin = xmin;
  s.xmax = xmax;
  s.ymin = ymin;
  s.ymax = ymax;
  s.zmin = zmin;
  s.zmax = zmax;
  g.members.push_back(s);
  return g;
}

SelectorGroup nodes_of(std::vector<Index> ids) {
  SelectorGroup g;
  Selector s;
  s.kind = SelectorKind::NodeIds;
  s.ids = std::move(ids);
  g.members.push_back(s);
  return g;
}

DisplacementConstraint fix(const SelectorGroup& region, int comp, Scalar value = 0.0) {
  DisplacementConstraint c;
  c.region = region;
  c.set(comp, true, value);
  return c;
}

std::string element_name(ElementType type) { return to_string(type); }

/// Positions on [0, length]: uniform steps h0 up to x1, then steps growing
/// by `growth`, the last one landing on `length`.
std::vector<Scalar> graded(Scalar h0, Scalar x1, Scalar length, Scalar growth) {
  const int n1 = std::max(1, static_cast<int>(std::lround(x1 / h0)));
  const Scalar h = x1 / n1;
  std::vector<Scalar> x;
  for (int i = 0; i <= n1; ++i) x.push_back(i * h);
  Scalar step = h;
  while (x.back() < length * (1.0 - 1.0e-12)) {
    step *= growth;
    if (x.back() + step * (1.0 + growth) > length) {
      x.push_back(length);
      break;
    }
    x.push_back(x.back() + step);
  }
  return x;
}

/// A structured mesh with unit cells (node (i, j[, k]) at integer
/// coordinates) whose nodes are then moved by `map`: a monotone map keeps
/// every cell's orientation.
Mesh remapped(ElementType type, Index nx, Index ny, Index nz,
              const std::function<Vector3(Index, Index, Index)>& map) {
  StructuredMeshSpec spec;
  spec.nx = nx;
  spec.ny = ny;
  spec.nz = nz;
  spec.lx = static_cast<Scalar>(nx);
  spec.ly = static_cast<Scalar>(ny);
  spec.lz = static_cast<Scalar>(nz);
  const Mesh grid = type == ElementType::Quad4 ? make_structured_quad_mesh(spec)
                                               : make_structured_hex_mesh(spec);
  Matrix coords = grid.coordinates();
  for (Index n = 0; n < grid.num_nodes(); ++n) {
    const Vector3 p = grid.node(n);
    const Vector3 x = map(static_cast<Index>(std::lround(p.x())),
                          static_cast<Index>(std::lround(p.y())),
                          static_cast<Index>(std::lround(p.z())));
    coords.col(n) = x.head(grid.dim());
  }
  return Mesh(coords, grid.connectivity(), type);
}

Mesh merge(const Mesh& a, const Mesh& b) {
  Matrix coords(a.dim(), a.num_nodes() + b.num_nodes());
  coords << a.coordinates(), b.coordinates();
  std::vector<Index> connectivity = a.connectivity();
  for (Index n : b.connectivity()) connectivity.push_back(n + a.num_nodes());
  return Mesh(coords, connectivity, a.element_type());
}

/// The lower cap of a cylinder (2-D) or sphere (3-D) of radius R touching
/// y = 0 at the origin: x (and z) in [0, W] graded from the origin, y from
/// the circle y_b = R - sqrt(R^2 - r^2) up to H, graded from the bottom.
Mesh cap_mesh(int dim, Scalar radius, Scalar width, Scalar height, Scalar h0, Scalar fine,
              Scalar growth) {
  const std::vector<Scalar> xs = graded(h0, fine, width, growth);
  const std::vector<Scalar> ys = graded(h0, fine, height, growth);
  const Index n = static_cast<Index>(xs.size()) - 1;
  const Index m = static_cast<Index>(ys.size()) - 1;
  return remapped(dim == 2 ? ElementType::Quad4 : ElementType::Hex8, n, m, n,
                  [&](Index i, Index j, Index k) {
                    const Scalar x = xs[static_cast<std::size_t>(i)];
                    const Scalar z = dim == 3 ? xs[static_cast<std::size_t>(k)] : 0.0;
                    const Scalar yb = radius - std::sqrt(radius * radius - x * x - z * z);
                    const Scalar y =
                        yb + (height - yb) * ys[static_cast<std::size_t>(j)] / height;
                    return Vector3(x, y, z);
                  });
}

/// A block [0, W] x [-H, 0] (x [0, W] in z), graded towards the origin of
/// its top face.
Mesh block_mesh(int dim, Scalar width, Scalar height, Scalar h0, Scalar fine, Scalar growth) {
  const std::vector<Scalar> xs = graded(h0, fine, width, growth);
  const std::vector<Scalar> ys = graded(h0, fine, height, growth);
  const Index n = static_cast<Index>(xs.size()) - 1;
  const Index m = static_cast<Index>(ys.size()) - 1;
  return remapped(dim == 2 ? ElementType::Quad4 : ElementType::Hex8, n, m, n,
                  [&](Index i, Index j, Index k) {
                    const Scalar x = xs[static_cast<std::size_t>(i)];
                    const Scalar z = dim == 3 ? xs[static_cast<std::size_t>(k)] : 0.0;
                    // j = 0 at the bottom, m at the top face y = 0.
                    const Scalar y = -ys[static_cast<std::size_t>(m - j)];
                    return Vector3(x, y, z);
                  });
}

NonlinearOptions contact_run(const std::vector<ContactPairSpec>& pairs, int steps = 1) {
  NonlinearOptions o;
  o.kinematics = Kinematics::SmallStrain;
  o.steps = steps;
  o.residual_tolerance = 1.0e-11;
  o.displacement_tolerance = 1.0e-11;
  o.max_iterations = 60;
  o.contact.enabled = true;
  o.contact.pairs = pairs;
  return o;
}

// ---------------------------------------------------------------------------
// contact-patch
// ---------------------------------------------------------------------------
Mesh patch_block(ElementType type, Index nx, Index ny, Index nz, Scalar lx, Scalar ly, Scalar lz,
                 Scalar y0, unsigned int seed) {
  StructuredMeshSpec spec;
  spec.nx = nx;
  spec.ny = ny;
  spec.nz = nz;
  spec.lx = lx;
  spec.ly = ly;
  spec.lz = lz;
  spec.y0 = y0;
  switch (type) {
    case ElementType::Quad4: return make_perturbed_quad_mesh(spec, 0.25, seed);
    case ElementType::Tri3: return make_perturbed_tri_mesh(spec, 0.25, seed);
    case ElementType::Hex8: return make_perturbed_hex_mesh(spec, 0.2, seed);
    default: return make_perturbed_tet_mesh(spec, 0.2, seed);
  }
}

}  // namespace

StudyOutcome study_contact_patch(const std::string& out_dir, json::Value& summary) {
  const IsotropicMaterial m(70.0e9, 0.3, 2700.0, "aluminium");
  const Scalar e = m.youngs_modulus();
  const Scalar nu = m.poisson_ratio();
  const Scalar lx = 0.4;
  const Scalar lz = 0.3;
  const Scalar h1 = 0.2;
  const Scalar h2 = 0.15;
  const Scalar delta = 1.0e-4;
  const Scalar g0 = 0.25 * delta;
  CsvWriter csv(path_join(out_dir, "contact_patch.csv"),
                {"element", "case", "nodes_in_contact", "pressure_error[-]",
                 "displacement_error[-]", "max_gap[m]", "tangential_over_normal[-]"});
  json::Value block = json::Value::make_object();
  json::Value cases = json::Value::make_array();
  Scalar worst = 0.0;
  bool all_completed = true;
  for (const ElementType type :
       {ElementType::Quad4, ElementType::Tri3, ElementType::Hex8, ElementType::Tet4}) {
    for (const std::string& kind : {std::string("rigid plane"), std::string("rigid plane, gap"),
                                    std::string("mortar, non-matching"),
                                    std::string("rigid plane, full slip")}) {
      const bool mortar = kind.rfind("mortar", 0) == 0;
      const bool gap = kind == "rigid plane, gap";
      const bool slip = kind == "rigid plane, full slip";
      if (slip && (type == ElementType::Tri3 || type == ElementType::Tet4)) continue;
      const Mesh lower = patch_block(type, 4, 3, 3, lx, h1, lz, 0.0, 11u);
      Mesh mesh = mortar ? merge(lower, patch_block(type, 5, 2, 4, lx, h2, lz, h1 + g0, 17u))
                         : patch_block(type, 4, 3, 3, lx, h1, lz, 0.0, 11u);
      const int dim = mesh.dim();
      const Scalar thickness = dim == 2 ? 0.01 : 1.0;
      FemModel model(std::move(mesh), m, thickness,
                     dim == 2 ? StressState::PlaneStress : StressState::ThreeDimensional,
                     IntegrationOptions());
      const Scalar top = mortar ? h1 + g0 + h2 : h1;
      const SelectorGroup top_face = box(-kInf, kInf, top, kInf);
      model.constraints().push_back(fix(top_face, 1, -delta));
      std::vector<Scalar> corners{0.0};
      if (mortar) {
        model.constraints().push_back(fix(box(-kInf, kInf, -kInf, 0.0), 1));
        corners.push_back(h1 + g0);
      }
      if (slip) {
        model.constraints().push_back(fix(top_face, 0, 5.0 * delta));
        if (dim == 3) model.constraints().push_back(fix(top_face, 2));
      } else {
        for (const Scalar y : corners) {
          model.constraints().push_back(fix(box(-kInf, 0.0, y, y, -kInf, 0.0), 0));
          if (dim == 3) {
            model.constraints().push_back(fix(box(-kInf, 0.0, y, y, -kInf, 0.0), 2));
            model.constraints().push_back(fix(box(lx, kInf, y, y, -kInf, 0.0), 2));
          }
        }
      }
      LoadCaseSpec lc;
      lc.name = "press";
      lc.prescribed_displacement_only = true;
      model.load_case_specs().push_back(lc);
      model.finalize();
      Assembler assembler(model);
      ContactPairSpec pair;
      pair.name = "interface";
      if (mortar) {
        pair.rigid = false;
        pair.slave = box(-kInf, kInf, h1 + g0, h1 + g0);
        pair.master = box(-kInf, kInf, h1, h1);
      } else {
        pair.slave = box(-kInf, kInf, -kInf, 0.0);
        pair.obstacle.kind = RigidObstacle::Kind::Plane;
        pair.obstacle.point = Vector3(0.0, gap ? -g0 : 0.0, 0.0);
        pair.obstacle.direction = Vector3::UnitY();
      }
      if (slip) pair.friction = 0.3;
      NonlinearOptions options = contact_run({pair}, slip ? 4 : 1);
      NonlinearMonitor fx;
      fx.name = "fx";
      fx.region = top_face;
      fx.component = 0;
      fx.quantity = NonlinearMonitor::Quantity::Reaction;
      NonlinearMonitor fy = fx;
      fy.name = "fy";
      fy.component = 1;
      options.monitors = {fx, fy};
      NonlinearStaticAnalysis analysis(model, assembler, options);
      const NonlinearResult r = analysis.solve(0);
      all_completed = all_completed && r.completed;
      // Exact: uniaxial stress E strain in every block, strain =
      // (delta - g) / (total height) with g the closed gap.
      const Scalar closed = gap || mortar ? g0 : 0.0;
      const Scalar height = mortar ? h1 + h2 : h1;
      const Scalar strain = (delta - closed) / height;
      const Scalar p = e * strain;
      Scalar pressure_error = 0.0;
      Scalar max_gap = 0.0;
      int active = 0;
      for (const ContactNodeResult& c : r.contact_nodes) {
        if (c.status != ContactStatus::Open) ++active;
        max_gap = std::max(max_gap, std::abs(c.gap));
        if (!slip) pressure_error = std::max(pressure_error, std::abs(c.pressure - p) / p);
      }
      if (active != static_cast<int>(r.contact_nodes.size())) pressure_error = 1.0;
      Scalar displacement_error = 0.0;
      Scalar ratio = 0.0;
      if (!slip) {
        for (Index n = 0; n < model.mesh().num_nodes(); ++n) {
          const Vector3 x = model.mesh().node(n);
          const bool upper = mortar && n >= lower.num_nodes();
          Vector3 exact = Vector3::Zero();
          exact.y() = upper ? -closed - strain * (x.y() - g0) : -closed * (mortar ? 0.0 : 1.0) -
                                                                     strain * x.y();
          exact.x() = nu * strain * x.x();
          if (dim == 3) exact.z() = nu * strain * x.z();
          for (int k = 0; k < dim; ++k) {
            displacement_error = std::max(
                displacement_error, std::abs(r.displacement(n * dim + k) - exact(k)) / delta);
          }
        }
      } else {
        // Full slip: every node slips with mu p against the slip, and in
        // 2-D the tangential force is mu times the normal one.
        const Scalar rx = r.steps.back().monitors[0];
        const Scalar ry = r.steps.back().monitors[1];
        ratio = std::abs(rx / ry);
        int slipping = 0;
        for (const ContactNodeResult& c : r.contact_nodes) {
          if (c.status == ContactStatus::Slip) ++slipping;
          pressure_error = std::max(pressure_error,
                                    std::abs(c.traction.norm() - 0.3 * c.pressure) /
                                        std::max(c.pressure, 1.0e-300));
        }
        if (slipping != static_cast<int>(r.contact_nodes.size())) pressure_error = 1.0;
        if (dim == 2) displacement_error = std::abs(ratio - 0.3) / 0.3;
      }
      worst = std::max({worst, pressure_error, displacement_error});
      csv.raw_row({element_name(type), kind, std::to_string(active), fmt(pressure_error, 4),
                   fmt(displacement_error, 4), fmt(max_gap, 4), slip ? fmt(ratio, 12) : ""});
      json::Value c = json::Value::make_object();
      c.set("element", json::Value::make_string(element_name(type)));
      c.set("case", json::Value::make_string(kind));
      c.set("completed", json::Value::make_bool(r.completed));
      c.set("nodes_in_contact", json::Value::make_number(active));
      c.set("pressure_error", json::Value::make_number(pressure_error));
      c.set("displacement_error", json::Value::make_number(displacement_error));
      c.set("max_gap_m", json::Value::make_number(max_gap));
      if (slip) c.set("tangential_over_normal", json::Value::make_number(ratio));
      cases.push_back(c);
    }
  }
  csv.close();
  block.set("kind", json::Value::make_string(
                        "verification (exact homogeneous states on distorted meshes)"));
  block.set("cases", cases);
  block.set("note",
            json::Value::make_string(
                "Blocks 0.4 x 0.2 (x 0.3) m, E = 70 GPa, nu = 0.3, distorted interior "
                "nodes; the top pushed down by 1e-4 m. Rigid plane: frictionless, free "
                "to expand sideways, so uniaxial stress E (delta - g0) / H: the pressure "
                "at every node and the linear displacement field are exact (errors over "
                "the pressure and over delta). Mortar: a second block of 5 x 2 (x 4) "
                "cells on it, 2.5e-5 m above, its meshes not matching the lower one's; "
                "the same stress in both blocks. Full slip: the top also pushed sideways "
                "by 5e-4 m with mu = 0.3 - every node must slip with |t| = mu p, and in "
                "2-D the tangential force is then mu times the normal force (its relative "
                "error in the displacement_error column)."));
  summary.set("contact_patch", block);

  std::ostringstream note;
  note << "largest error " << fmt(worst, 3) << " over the pressures, displacements and the "
       << "full-slip force ratio of 14 cases (Q4, Tri3, Hex8, Tet4)";
  StudyOutcome outcome;
  outcome.name = "contact patch tests: rigid plane, mortar with non-matching meshes, full slip";
  outcome.kind = "verification";
  outcome.metric = "largest relative error of pressure, displacement and slip force";
  outcome.value = worst;
  outcome.tolerance = 1.0e-9;
  outcome.passed = all_completed && worst <= 1.0e-9;
  outcome.note = note.str();
  return outcome;
}

namespace {

/// One Hertz configuration solved on one mesh.
struct HertzResult {
  Scalar h = 0.0;               ///< element size of the slave surface in the contact zone [m]
  Scalar delta = 0.0;           ///< prescribed approach [m]
  Scalar load = 0.0;            ///< total load P: per unit length (2-D) [N/m] or [N] (3-D)
  Scalar a = 0.0;               ///< Hertz half-width (radius) of the load P [m]
  Scalar p0 = 0.0;              ///< Hertz peak pressure of P [Pa]
  Scalar centre_error = 0.0;    ///< |p(0) - p0| / p0
  /// RMS of p - p_Hertz over the nodes within 0.8 a, over the RMS of p_Hertz
  /// there (each node weighted by D_j): the smooth interior of the zone.
  Scalar interior_error = 0.0;
  Scalar rms_error = 0.0;       ///< the same over the whole slave surface
  Scalar last_in = 0.0;         ///< largest r of a node in contact [m]
  Scalar first_out = 0.0;       ///< smallest r of an open node beyond it [m]
  /// The edge of the contact zone is resolved to an element: the Hertz a
  /// lies within h of the interval between the last node in contact and the
  /// first open one.
  bool edge_resolved = false;
  int in_contact = 0;
  int iterations = 0;
  bool completed = false;
  std::vector<Scalar> r, pressure, hertz;
};

enum class HertzCase { ElasticOnRigid, RigidOnElastic, ElasticPair };

std::string case_name(int dim, HertzCase c) {
  const std::string body = dim == 2 ? "cylinder" : "sphere";
  switch (c) {
    case HertzCase::ElasticOnRigid: return "elastic " + body + " on a rigid flat";
    case HertzCase::RigidOnElastic: return "rigid " + body + " into an elastic block";
    case HertzCase::ElasticPair:
      return "elastic " + body + " on an elastic block (mortar, non-matching)";
  }
  return "";
}

/// The element size of the mortar pair's master (the curved cap) over its
/// slave's (the block's flat top): the two meshes do not match across the
/// interface. A curved master is a polygon (a polyhedron) whose facets the
/// slave nodes between its vertices see; meshed finer than the slave, the
/// error this adds falls at second order, coarser at first order only (the
/// coarse-master series of `hertz-line`).
constexpr Scalar kMasterRatio = 0.75;

/// Solve one configuration: R the radius, W = H the size of each body, h0 the
/// element size of the slave surface in the contact zone (uniform up to
/// `fine`), the bodies pushed together by delta.
HertzResult solve_hertz(int dim, HertzCase kind, Scalar radius, Scalar size, Scalar h0,
                        Scalar fine, Scalar delta, const IsotropicMaterial& m,
                        Scalar master_ratio = kMasterRatio) {
  const Scalar growth = dim == 2 ? 1.25 : 1.35;
  Mesh mesh;
  Index cap_nodes = 0;
  switch (kind) {
    case HertzCase::ElasticOnRigid:
      mesh = cap_mesh(dim, radius, size, size, h0, fine, growth);
      cap_nodes = mesh.num_nodes();
      break;
    case HertzCase::RigidOnElastic: mesh = block_mesh(dim, size, size, h0, fine, growth); break;
    case HertzCase::ElasticPair: {
      const Mesh cap = cap_mesh(dim, radius, size, size, master_ratio * h0, fine, growth);
      cap_nodes = cap.num_nodes();
      mesh = merge(cap, block_mesh(dim, size, size, h0, fine, growth));
      break;
    }
  }
  const auto on_cap_bottom = [&](const Vector3& x) {
    const Scalar yb = radius - std::sqrt(radius * radius - x.x() * x.x() - x.z() * x.z());
    return std::abs(x.y() - yb) <= 1.0e-12 * radius;
  };
  std::vector<Index> cap_bottom;
  std::vector<Index> block_top;
  for (Index n = 0; n < mesh.num_nodes(); ++n) {
    const Vector3 x = mesh.node(n);
    const bool cap = n < cap_nodes;
    if (cap && on_cap_bottom(x)) cap_bottom.push_back(n);
    if (!cap && std::abs(x.y()) <= 1.0e-12 * radius) block_top.push_back(n);
  }
  FemModel model(std::move(mesh), m, 1.0,
                 dim == 2 ? StressState::PlaneStrain : StressState::ThreeDimensional,
                 IntegrationOptions());
  // Symmetry planes x = 0 (and z = 0), the cap pushed down at its top, the
  // block held at its bottom.
  model.constraints().push_back(fix(box(-kInf, 0.0), 0));
  if (dim == 3) model.constraints().push_back(fix(box(-kInf, kInf, -kInf, kInf, -kInf, 0.0), 2));
  if (kind != HertzCase::RigidOnElastic) {
    model.constraints().push_back(fix(box(-kInf, kInf, size, kInf), 1, -delta));
  }
  if (kind != HertzCase::ElasticOnRigid) {
    model.constraints().push_back(fix(box(-kInf, kInf, -kInf, -size), 1));
  }
  LoadCaseSpec lc;
  lc.name = "press";
  lc.prescribed_displacement_only = true;
  model.load_case_specs().push_back(lc);
  model.finalize();
  Assembler assembler(model);
  ContactPairSpec pair;
  pair.name = "contact";
  Scalar e_star = 0.0;
  const Scalar e = m.youngs_modulus();
  const Scalar nu = m.poisson_ratio();
  switch (kind) {
    case HertzCase::ElasticOnRigid:
      pair.slave = nodes_of(cap_bottom);
      pair.obstacle.kind = RigidObstacle::Kind::Plane;
      pair.obstacle.direction = Vector3::UnitY();
      e_star = e / (1.0 - nu * nu);
      break;
    case HertzCase::RigidOnElastic:
      pair.slave = box(-kInf, kInf, 0.0, kInf);
      pair.obstacle.kind =
          dim == 2 ? RigidObstacle::Kind::Cylinder : RigidObstacle::Kind::Sphere;
      pair.obstacle.point = Vector3(0.0, radius, 0.0);
      pair.obstacle.radius = radius;
      pair.obstacle.motion = Vector3(0.0, -delta, 0.0);
      e_star = e / (1.0 - nu * nu);
      break;
    case HertzCase::ElasticPair:
      pair.rigid = false;
      pair.slave = nodes_of(block_top);  // flat: its normals are Hertz's direction
      pair.master = nodes_of(cap_bottom);
      e_star = e / (2.0 * (1.0 - nu * nu));
      break;
  }
  NonlinearStaticAnalysis analysis(model, assembler, contact_run({pair}));
  const NonlinearResult r = analysis.solve(0);
  HertzResult out;
  const std::vector<Scalar> xs = graded(h0, fine, size, growth);
  out.h = xs[1] - xs[0];
  out.delta = delta;
  out.completed = r.completed;
  out.iterations = r.total_iterations;
  if (!r.completed || r.contact_pairs.empty()) return out;
  // The full model's load from the half (quarter) model's contact force.
  const Scalar quarter = dim == 2 ? 2.0 : 4.0;
  out.load = quarter * std::abs(r.contact_pairs[0].force.y());
  if (dim == 2) {
    out.a = std::sqrt(4.0 * out.load * radius / (kPi * e_star));
    out.p0 = 2.0 * out.load / (kPi * out.a);
  } else {
    out.a = std::cbrt(3.0 * out.load * radius / (4.0 * e_star));
    out.p0 = 3.0 * out.load / (2.0 * kPi * out.a * out.a);
  }
  Scalar num = 0.0;
  Scalar den = 0.0;
  Scalar inner_num = 0.0;
  Scalar inner_den = 0.0;
  out.first_out = kInf;
  std::vector<std::pair<Scalar, const ContactNodeResult*>> sorted;
  for (const ContactNodeResult& c : r.contact_nodes) {
    const Vector3 x = model.mesh().node(c.node);
    sorted.emplace_back(std::hypot(x.x(), x.z()), &c);
  }
  std::sort(sorted.begin(), sorted.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
  for (const auto& [rr, c] : sorted) {
    const Scalar ph = rr < out.a ? out.p0 * std::sqrt(1.0 - rr * rr / (out.a * out.a)) : 0.0;
    const Scalar d2 = (c->pressure - ph) * (c->pressure - ph);
    num += c->weight * d2;
    den += c->weight * ph * ph;
    if (rr <= 0.8 * out.a) {
      inner_num += c->weight * d2;
      inner_den += c->weight * ph * ph;
    }
    if (rr <= 1.0e-12 * radius) out.centre_error = std::abs(c->pressure - out.p0) / out.p0;
    if (c->status != ContactStatus::Open) {
      ++out.in_contact;
      out.last_in = std::max(out.last_in, rr);
    }
    if (rr <= 3.0 * out.a) {
      out.r.push_back(rr);
      out.pressure.push_back(c->pressure);
      out.hertz.push_back(ph);
    }
  }
  for (const auto& [rr, c] : sorted) {
    if (c->status == ContactStatus::Open && rr > out.last_in) {
      out.first_out = std::min(out.first_out, rr);
    }
  }
  out.rms_error = std::sqrt(num / den);
  out.interior_error = std::sqrt(inner_num / inner_den);
  out.edge_resolved = out.last_in - out.h <= out.a && out.a <= out.first_out + out.h;
  return out;
}

/// The approach that gives a contact half-width near `a_target`: one solve
/// on a coarse mesh, the approach scaled by (a_target / a)^2 (P grows
/// about in proportion to the approach and a with its square root).
Scalar calibrate(int dim, HertzCase kind, Scalar radius, Scalar size, Scalar h0, Scalar fine,
                 Scalar a_target, const IsotropicMaterial& m,
                 Scalar master_ratio = kMasterRatio) {
  const Scalar guess = a_target * a_target / radius;
  const HertzResult first =
      solve_hertz(dim, kind, radius, size, h0, fine, guess, m, master_ratio);
  if (!first.completed || !(first.a > 0.0)) return guess;
  return guess * std::pow(a_target / first.a, 2.0);
}

std::vector<std::string> hertz_columns(int dim) {
  return {"case", "series", "R_over_a_target[-]", "L_over_a_target[-]", "h[m]", "a_over_h[-]",
          "approach[m]", dim == 2 ? "load[N/m]" : "load[N]", "a[m]", "p0[Pa]",
          "nodes_in_contact", "last_in[m]", "first_out[m]", "edge_within_an_element",
          "centre_error[-]", "interior_error[-]", "rms_error[-]", "newton_iterations"};
}

void write_hertz(CsvWriter& csv, CsvWriter* profile, const std::string& label,
                 const std::string& series, Scalar r_over_a, Scalar l_over_a,
                 const HertzResult& r) {
  csv.raw_row({label, series, fmt(r_over_a, 6), fmt(l_over_a, 6), fmt(r.h, 6),
               fmt(r.a / r.h, 6), fmt(r.delta, 8), fmt(r.load, 8), fmt(r.a, 8), fmt(r.p0, 8),
               std::to_string(r.in_contact), fmt(r.last_in, 8),
               std::isfinite(r.first_out) ? fmt(r.first_out, 8) : "",
               r.edge_resolved ? "yes" : "no", fmt(r.centre_error, 6),
               fmt(r.interior_error, 6), fmt(r.rms_error, 6), std::to_string(r.iterations)});
  if (profile == nullptr) return;
  for (std::size_t i = 0; i < r.r.size(); ++i) {
    profile->raw_row({label, fmt(r.h, 6), fmt(r.a / r.h, 6), fmt(r.r[i] / r.a, 10),
                      fmt(r.pressure[i] / r.p0, 10), fmt(r.hertz[i] / r.p0, 10)});
  }
}

json::Value hertz_json(const HertzResult& r) {
  json::Value o = json::Value::make_object();
  o.set("h_m", json::Value::make_number(r.h));
  o.set("approach_m", json::Value::make_number(r.delta));
  o.set("load", json::Value::make_number(r.load));
  o.set("a_m", json::Value::make_number(r.a));
  o.set("p0_Pa", json::Value::make_number(r.p0));
  o.set("a_over_h", json::Value::make_number(r.a / r.h));
  o.set("nodes_in_contact", json::Value::make_number(r.in_contact));
  o.set("last_node_in_contact_m", json::Value::make_number(r.last_in));
  if (std::isfinite(r.first_out)) {
    o.set("first_open_node_m", json::Value::make_number(r.first_out));
  }
  o.set("edge_within_an_element", json::Value::make_bool(r.edge_resolved));
  o.set("centre_pressure_error", json::Value::make_number(r.centre_error));
  o.set("interior_pressure_error", json::Value::make_number(r.interior_error));
  o.set("rms_pressure_error", json::Value::make_number(r.rms_error));
  o.set("newton_iterations", json::Value::make_number(r.iterations));
  o.set("completed", json::Value::make_bool(r.completed));
  return o;
}

}  // namespace

StudyOutcome study_hertz_line(const std::string& out_dir, json::Value& summary) {
  const IsotropicMaterial m(200.0e9, 0.3, 7850.0, "steel");
  const Scalar a_target = 1.0e-3;
  const Scalar fine = 1.5 * a_target;
  const std::vector<Scalar> divisions{5.0, 10.0, 20.0, 40.0};
  CsvWriter csv(path_join(out_dir, "hertz_line.csv"), hertz_columns(2));
  CsvWriter profile(path_join(out_dir, "hertz_line_profile.csv"),
                    {"case", "h[m]", "a_over_h[-]", "x_over_a[-]", "p_over_p0[-]",
                     "hertz_over_p0[-]"});
  json::Value block = json::Value::make_object();
  json::Value cases = json::Value::make_array();
  bool ok = true;
  Scalar worst = 0.0;       // centre and interior errors, finest meshes
  Scalar worst_flat = 0.0;  // the same, flat obstacle and mortar pair
  std::ostringstream note;
  // The mesh series: R = 50 a, bodies 25 a across, a / h from 5 to 40. A
  // curved rigid obstacle presses along its own normal, tilted by x / R
  // from the direction Hertz assumes: a difference of order a / R, hence its
  // own tolerance.
  const Scalar radius = 50.0 * a_target;
  const Scalar size = 25.0 * a_target;
  for (const HertzCase kind :
       {HertzCase::ElasticOnRigid, HertzCase::RigidOnElastic, HertzCase::ElasticPair}) {
    const std::string label = case_name(2, kind);
    const Scalar tolerance = kind == HertzCase::RigidOnElastic ? 3.0e-3 : 1.0e-3;
    const Scalar delta =
        calibrate(2, kind, radius, size, a_target / divisions.front(), fine, a_target, m);
    json::Value c = json::Value::make_object();
    c.set("case", json::Value::make_string(label));
    c.set("tolerance", json::Value::make_number(tolerance));
    json::Value meshes = json::Value::make_array();
    std::vector<HertzResult> runs;
    for (const Scalar d : divisions) {
      runs.push_back(solve_hertz(2, kind, radius, size, a_target / d, fine, delta, m));
      write_hertz(csv, &profile, label, "mesh", 50.0, 25.0, runs.back());
      meshes.push_back(hertz_json(runs.back()));
      ok = ok && runs.back().completed && runs.back().edge_resolved;
    }
    const HertzResult& finest = runs.back();
    const Scalar error = std::max(finest.centre_error, finest.interior_error);
    ok = ok && error <= tolerance && finest.interior_error < runs.front().interior_error;
    worst = std::max(worst, error);
    if (kind != HertzCase::RigidOnElastic) worst_flat = std::max(worst_flat, error);
    c.set("meshes", meshes);
    cases.push_back(c);
    note << label << ": centre " << fmt(finest.centre_error, 3) << ", interior "
         << fmt(finest.interior_error, 3) << ", RMS " << fmt(finest.rms_error, 3) << "; ";
  }
  block.set("cases", cases);

  // The floor the refinement reaches is the difference between the finite
  // model and Hertz's half-space theory. Two series at a / h = 40 show where
  // it comes from: the elastic cylinder on the rigid flat as the bodies grow
  // (R = 200 a, so a / R plays no part), and the rigid cylinder as a / R
  // shrinks (bodies 100 a across).
  const auto floor_series = [&](HertzCase kind, const std::string& series,
                                const std::vector<std::pair<Scalar, Scalar>>& setups,
                                json::Value& out_runs) {
    std::vector<Scalar> errors;
    for (const auto& [r_over_a, l_over_a] : setups) {
      const Scalar rr = r_over_a * a_target;
      const Scalar ll = l_over_a * a_target;
      const Scalar delta = calibrate(2, kind, rr, ll, a_target / 10.0, fine, a_target, m);
      const HertzResult r = solve_hertz(2, kind, rr, ll, a_target / 40.0, fine, delta, m);
      write_hertz(csv, nullptr, case_name(2, kind), series, r_over_a, l_over_a, r);
      json::Value j = hertz_json(r);
      j.set("R_over_a_target", json::Value::make_number(r_over_a));
      j.set("L_over_a_target", json::Value::make_number(l_over_a));
      out_runs.push_back(j);
      ok = ok && r.completed;
      errors.push_back(r.centre_error);
    }
    bool falls = true;
    for (std::size_t i = 1; i < errors.size(); ++i) falls = falls && errors[i] < errors[i - 1];
    return std::make_pair(falls, errors);
  };
  json::Value size_runs = json::Value::make_array();
  const auto [size_falls, size_errors] =
      floor_series(HertzCase::ElasticOnRigid, "body size",
                   {{200.0, 25.0}, {200.0, 50.0}, {200.0, 100.0}}, size_runs);
  json::Value curvature_runs = json::Value::make_array();
  const auto [curvature_falls, curvature_errors] =
      floor_series(HertzCase::RigidOnElastic, "curvature",
                   {{50.0, 100.0}, {100.0, 100.0}, {200.0, 100.0}}, curvature_runs);
  ok = ok && size_falls && curvature_falls;
  const Scalar curvature_order = observed_order(1.0 / 50.0, curvature_errors.front(),
                                                1.0 / 200.0, curvature_errors.back());
  json::Value floors = json::Value::make_object();
  floors.set("body_size", size_runs);
  floors.set("curvature", curvature_runs);
  floors.set("curvature_error_order_in_a_over_R", json::Value::make_number(curvature_order));
  block.set("model_floor", floors);

  // The mortar pair with its curved master meshed twice as coarse as its
  // slave: the slave nodes between the master's vertices see its chords, an
  // error of first order in the element size that must still fall with it.
  json::Value coarse_runs = json::Value::make_array();
  std::vector<Scalar> coarse_errors;
  {
    const HertzCase kind = HertzCase::ElasticPair;
    const Scalar ratio = 2.0;
    const Scalar delta =
        calibrate(2, kind, radius, size, a_target / 10.0, fine, a_target, m, ratio);
    std::vector<Scalar> hs;
    for (const Scalar d : {10.0, 20.0, 40.0}) {
      const HertzResult r =
          solve_hertz(2, kind, radius, size, a_target / d, fine, delta, m, ratio);
      write_hertz(csv, nullptr, case_name(2, kind), "master twice as coarse", 50.0, 25.0, r);
      coarse_runs.push_back(hertz_json(r));
      ok = ok && r.completed && r.edge_resolved;
      if (!coarse_errors.empty()) ok = ok && r.centre_error < coarse_errors.back();
      coarse_errors.push_back(r.centre_error);
      hs.push_back(r.h);
    }
    const Scalar order = observed_order(hs[1], coarse_errors[1], hs[2], coarse_errors[2]);
    json::Value coarse = json::Value::make_object();
    coarse.set("master_over_slave_element_size", json::Value::make_number(ratio));
    coarse.set("meshes", coarse_runs);
    coarse.set("centre_error_order_finest_pair", json::Value::make_number(order));
    block.set("coarse_master", coarse);
  }
  csv.close();
  profile.close();
  block.set("kind", json::Value::make_string(
                        "verification + validation (Hertz's half-space theory, approached as "
                        "the mesh is refined, down to the gap between the finite model and "
                        "the theory, which two further series identify)"));
  block.set("note",
            json::Value::make_string(
                "Plane strain, steel (E = 200 GPa, nu = 0.3); half models (symmetry at x = "
                "0) graded outward at a ratio of 1.25 from a uniform zone of 1.5 a_target, "
                "a_target = 1 mm, the approach chosen by one coarse solve so that a is near "
                "a_target. From the total load P per unit length: a = sqrt(4 P R / (pi "
                "E*)), p0 = 2 P / (pi a), E* = E / (1 - nu^2) against a rigid body and E / "
                "(2 (1 - nu^2)) for two of the same material. Mesh series: R = 50 a_target, "
                "bodies 25 a_target wide and deep, a / h from 5 to 40; the mortar pair's "
                "master (the cylinder) meshed 4/3 as finely as its slave (the block's "
                "flat top). Errors: at the centre node, and the RMS of p - p_Hertz over the "
                "nodes within 0.8 a (the smooth interior) and over the whole surface, each "
                "node weighted by its D_j; the edge is resolved when a lies within one "
                "element of the interval between the last node in contact and the first "
                "open one. Model floor, a / h = 40: the elastic cylinder on the rigid flat "
                "with R = 200 a_target and bodies 25, 50 and 100 a_target across (the finite "
                "body), and the rigid cylinder into a block 100 a_target across with R = "
                "50, 100 and 200 a_target (its contact force along the cylinder's normal, "
                "tilted by x / R from the vertical Hertz assumes). Coarse master: the mortar "
                "pair with its master meshed twice as coarse as its slave, a / h = 10, 20 "
                "and 40 - the slave nodes between the polygonal master's vertices see its "
                "chords, an error of first order in the element size."));
  summary.set("hertz_line", block);

  note << "model floor: centre error " << fmt(size_errors.front(), 3) << " -> "
       << fmt(size_errors.back(), 3) << " as the bodies grow from 25 a to 100 a, "
       << fmt(curvature_errors.front(), 3) << " -> " << fmt(curvature_errors.back(), 3)
       << " as a / R falls from 1/50 to 1/200 (order " << fmt(curvature_order, 3)
       << "); master twice as coarse as the slave: centre error " << fmt(coarse_errors.front(), 3)
       << " -> " << fmt(coarse_errors.back(), 3) << " from a / h = 10 to 40";
  StudyOutcome outcome;
  outcome.name = "Hertz line contact: cylinder on flat, rigid and elastic (Q4, plane strain)";
  outcome.kind = "verification + validation";
  outcome.metric =
      "largest centre / interior pressure error against Hertz, finest meshes (a / h ~ 42)";
  outcome.value = worst;
  outcome.tolerance = 3.0e-3;
  outcome.passed = ok && worst <= 3.0e-3 && worst_flat <= 1.0e-3;
  outcome.note = note.str();
  return outcome;
}

StudyOutcome study_hertz_point(const std::string& out_dir, json::Value& summary) {
  const IsotropicMaterial m(200.0e9, 0.3, 7850.0, "steel");
  const Scalar a_target = 1.0e-3;
  const Scalar radius = 50.0 * a_target;
  const Scalar size = 15.0 * a_target;
  const Scalar fine = 1.5 * a_target;
  CsvWriter csv(path_join(out_dir, "hertz_point.csv"), hertz_columns(3));
  CsvWriter profile(path_join(out_dir, "hertz_point_profile.csv"),
                    {"case", "h[m]", "a_over_h[-]", "r_over_a[-]", "p_over_p0[-]",
                     "hertz_over_p0[-]"});
  json::Value block = json::Value::make_object();
  json::Value cases = json::Value::make_array();
  bool ok = true;
  std::ostringstream note;
  // The sphere on the rigid flat to a / h = 9; the mortar pair, with twice
  // the unknowns and more, to a / h = 6.
  std::vector<std::vector<HertzResult>> all;
  for (const HertzCase kind : {HertzCase::ElasticOnRigid, HertzCase::ElasticPair}) {
    const bool single = kind == HertzCase::ElasticOnRigid;
    const std::string label = case_name(3, kind);
    const std::vector<Scalar> divisions =
        single ? std::vector<Scalar>{4.0, 6.0, 9.0} : std::vector<Scalar>{4.0, 6.0};
    const Scalar delta =
        calibrate(3, kind, radius, size, a_target / divisions.front(), fine, a_target, m);
    std::vector<HertzResult> runs;
    json::Value meshes = json::Value::make_array();
    for (const Scalar d : divisions) {
      runs.push_back(solve_hertz(3, kind, radius, size, a_target / d, fine, delta, m));
      write_hertz(csv, &profile, label, "mesh", 50.0, 15.0, runs.back());
      meshes.push_back(hertz_json(runs.back()));
      ok = ok && runs.back().completed && runs.back().edge_resolved;
      if (single) ok = ok && runs.back().centre_error <= 3.0e-3;
    }
    for (std::size_t i = 1; i < runs.size(); ++i) {
      ok = ok && runs[i].rms_error < runs[i - 1].rms_error;
    }
    const std::size_t k = runs.size() - 2;
    const Scalar order =
        observed_order(runs[k].h, runs[k].rms_error, runs.back().h, runs.back().rms_error);
    const Scalar interior_order = observed_order(runs[k].h, runs[k].interior_error,
                                                 runs.back().h, runs.back().interior_error);
    if (!single) {
      // The facets of the non-matching curved master add an error that
      // falls at second order or faster.
      ok = ok && interior_order >= 2.0 && runs.back().centre_error <= 1.0e-2;
    }
    json::Value c = json::Value::make_object();
    c.set("case", json::Value::make_string(label));
    c.set("meshes", meshes);
    c.set("rms_order_finest_pair", json::Value::make_number(order));
    c.set("interior_order_finest_pair", json::Value::make_number(interior_order));
    cases.push_back(c);
    note << label << ": RMS " << fmt(runs.back().rms_error, 3) << " at a / h = "
         << fmt(runs.back().a / runs.back().h, 3) << " (order " << fmt(order, 3) << "), interior "
         << fmt(runs.back().interior_error, 3) << " (order " << fmt(interior_order, 3)
         << "), centre " << fmt(runs.back().centre_error, 3) << "; ";
    all.push_back(std::move(runs));
  }
  csv.close();
  profile.close();
  const Scalar finest = all[0].back().rms_error;
  block.set("cases", cases);
  block.set("kind", json::Value::make_string(
                        "verification + validation (Hertz's half-space theory)"));
  block.set("note",
            json::Value::make_string(
                "Steel (E = 200 GPa, nu = 0.3), R = 50 a_target, a_target = 1 mm; quarter "
                "models (symmetry at x = 0 and z = 0), each body 15 a_target wide and "
                "deep, Hex8 graded outward at a ratio of 1.35 from a uniform zone of 1.5 "
                "a_target; the approach chosen by one coarse solve so that a is near "
                "a_target. From the total load P: a = (3 P R / (4 E*))^(1/3), p0 = 3 P / (2 "
                "pi a^2), E* = E / (1 - nu^2) on the rigid flat, E / (2 (1 - nu^2)) for the "
                "mortar pair, whose master (the sphere) is meshed 4/3 as finely as its "
                "slave (the block's flat top), the two meshes not matching. In 3-D the RMS "
                "error over the surface is set by the square-root edge of the pressure at "
                "every angle and falls steadily. Required: the RMS error falls with each "
                "refinement and the edge lies within one element of Hertz's on every mesh; "
                "on the rigid flat the centre error stays within 3e-3 and the RMS error at "
                "a / h = 9 within 0.05; on the mortar pair the interior error falls at "
                "second order or faster (the facets of the curved master) and the centre "
                "error at a / h = 6 is within 1e-2."));
  summary.set("hertz_point", block);

  StudyOutcome outcome;
  outcome.name = "Hertz point contact: sphere on a rigid flat and a mortar pair (Hex8)";
  outcome.kind = "verification + validation";
  outcome.metric = "RMS pressure error against Hertz, sphere on the rigid flat, a / h = 9";
  outcome.value = finest;
  outcome.tolerance = 0.05;
  outcome.passed = ok && finest <= 0.05;
  outcome.note = note.str();
  return outcome;
}

}  // namespace verify
}  // namespace sparlab
