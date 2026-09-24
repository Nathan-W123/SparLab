#include "sparlab/elements/Element.hpp"

#include "sparlab/core/Exceptions.hpp"
#include "sparlab/elements/Hex8.hpp"
#include "sparlab/elements/Quad4.hpp"
#include "sparlab/elements/Tet4.hpp"
#include "sparlab/elements/Tri3.hpp"

#include <sstream>

namespace sparlab {

const std::vector<int>& Element::face_nodes(int local_face) const {
  const std::vector<std::vector<int>>& faces = element_local_faces(type());
  if (local_face < 0 || local_face >= static_cast<int>(faces.size())) {
    std::ostringstream os;
    os << to_string(type()) << " local face index " << local_face << " is outside [0, "
       << faces.size() - 1 << "]";
    throw MeshError(os.str());
  }
  return faces[static_cast<std::size_t>(local_face)];
}

std::unique_ptr<Element> make_element(ElementType type) {
  switch (type) {
    case ElementType::Quad4: return std::make_unique<Quad4Element>();
    case ElementType::Hex8: return std::make_unique<Hex8Element>();
    case ElementType::Tri3: return std::make_unique<Tri3Element>();
    case ElementType::Tet4: return std::make_unique<Tet4Element>();
  }
  throw ConfigError("no element implementation registered for the requested type");
}

}  // namespace sparlab
