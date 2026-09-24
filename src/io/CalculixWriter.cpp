#include "sparlab/io/CalculixWriter.hpp"

#include "sparlab/core/Exceptions.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>

namespace sparlab {
namespace {

std::string sanitise(const std::string& name) {
  std::string out;
  for (char c : name) {
    const bool ok = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
                    (c >= '0' && c <= '9') || c == '-' || c == '_';
    out.push_back(ok ? c : '_');
  }
  return out.empty() ? "unnamed" : out;
}

}  // namespace

std::string calculix_element_type(const FemModel& model) {
  switch (model.mesh().element_type()) {
    case ElementType::Quad4:
      return model.stress_state() == StressState::PlaneStrain ? "CPE4" : "CPS4";
    case ElementType::Hex8:
      return "C3D8";
    case ElementType::Tri3:
      return model.stress_state() == StressState::PlaneStrain ? "CPE3" : "CPS3";
    case ElementType::Tet4:
      return "C3D4";
  }
  throw IoError("no CalculiX element type for this mesh");
}

std::vector<std::string> write_calculix_decks(const FemModel& model, const std::string& stem,
                                              const std::string& case_name) {
  if (!model.finalized()) throw IoError("the model must be finalised before export");
  const Mesh& mesh = model.mesh();
  const int dim = mesh.dim();
  const std::string element = calculix_element_type(model);
  const std::vector<Vector>& loads = model.load_vectors();
  const std::vector<LoadCaseSpec>& specs = model.load_case_specs();

  std::vector<std::string> paths;
  for (std::size_t l = 0; l < loads.size(); ++l) {
    const std::string path = stem + "_" + sanitise(specs[l].name) + ".inp";
    std::ofstream out(path, std::ios::out | std::ios::trunc);
    if (!out) throw IoError("cannot open '" + path + "' for writing");
    out << std::setprecision(17);

    out << "*HEADING\n";
    out << "SparLab cross-validation export: " << case_name << " / " << specs[l].name
        << " (" << element << ", " << mesh.num_elements() << " elements)\n";

    out << "*NODE, NSET=NALL\n";
    for (Index n = 0; n < mesh.num_nodes(); ++n) {
      const Vector3 x = mesh.node(n);
      out << n + 1 << ", " << x.x() << ", " << x.y() << ", " << x.z() << "\n";
    }

    out << "*ELEMENT, TYPE=" << element << ", ELSET=EALL\n";
    for (Index e = 0; e < mesh.num_elements(); ++e) {
      const Index* nodes = mesh.element_nodes(e);
      out << e + 1;
      for (int a = 0; a < mesh.nodes_per_elem(); ++a) out << ", " << nodes[a] + 1;
      out << "\n";
    }

    out << "*MATERIAL, NAME=MAT\n*ELASTIC\n"
        << model.material().youngs_modulus() << ", " << model.material().poisson_ratio()
        << "\n";
    out << "*SOLID SECTION, ELSET=EALL, MATERIAL=MAT\n";
    if (dim == 2) out << model.thickness() << "\n";

    out << "*STEP\n*STATIC\n";
    out << "*BOUNDARY\n";
    for (Index d : model.dofs().constrained_dofs()) {
      const Index node = d / dim;
      const int component = static_cast<int>(d % dim) + 1;
      out << node + 1 << ", " << component << ", " << component << ", "
          << model.dofs().prescribed_value(d) << "\n";
    }
    out << "*CLOAD\n";
    for (Index n = 0; n < mesh.num_nodes(); ++n) {
      for (int k = 0; k < dim; ++k) {
        const Scalar f = loads[l](n * dim + k);
        if (f != 0.0) out << n + 1 << ", " << k + 1 << ", " << f << "\n";
      }
    }
    out << "*NODE FILE\nU\n*EL FILE\nS\n*END STEP\n";
    out.flush();
    if (!out) throw IoError("failed while writing '" + path + "'");
    paths.push_back(path);
  }
  return paths;
}

}  // namespace sparlab
