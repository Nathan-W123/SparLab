/// \file CalculixWriter.hpp
/// \brief Export of a finalised model as CalculiX (Abaqus-style) input decks.
///
/// One `.inp` per load case is written, each a complete linear static job:
/// nodes (1-based), elements (`CPS4` for plane stress, `CPE4` for plane
/// strain, `C3D8` for solids, in the same local node order SparLab uses),
/// an isotropic `*ELASTIC` material, a `*SOLID SECTION` (with the thickness
/// on its data line for plane elements), the prescribed DOFs as `*BOUNDARY`
/// cards, the assembled nodal forces as `*CLOAD` cards, and `*NODE FILE, U`
/// so the solver writes nodal displacements to its `.frd` result file.
///
/// The exported problem is the *same discrete problem* SparLab solves
/// (identical mesh, element type, integration order, material and consistent
/// nodal loads), so agreement is expected to solver precision for `C3D8` and
/// to the level at which CalculiX's plane elements - which it expands into
/// solid elements through the thickness - reproduce a plane-stress Q4. Both
/// are measured, not assumed; see python/scripts/cross_validate.py.
#pragma once

#include "sparlab/core/Types.hpp"
#include "sparlab/fem/FemModel.hpp"

#include <string>
#include <vector>

namespace sparlab {

/// Write `<stem>_<load case>.inp` for every load case of `model`.
/// \return the paths written, in load-case order.
/// \throws IoError when a file cannot be written.
std::vector<std::string> write_calculix_decks(const FemModel& model, const std::string& stem,
                                              const std::string& case_name);

/// CalculiX element keyword for the model's element type and stress state.
std::string calculix_element_type(const FemModel& model);

}  // namespace sparlab
