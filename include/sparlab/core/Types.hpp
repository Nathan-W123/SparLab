/// \file Types.hpp
/// \brief Fundamental scalar / matrix aliases and dimensional constants.
///
/// All quantities in SparLab are expressed in strict SI base units:
///   length [m], force [N], stress and modulus [Pa], density [kg/m^3],
///   mass [kg], frequency [Hz], compliance [J] = [N m].
/// See docs/conventions.md for the full convention table.
#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cstddef>
#include <vector>

namespace sparlab {

/// Floating point type used throughout the library.
using Scalar = double;

/// Integer type used for sparse-matrix indices (must match Eigen's default).
using StorageIndex = int;

/// Global node / element / DOF index type.
using Index = int;

using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;
using Vector2 = Eigen::Vector2d;
using Vector3 = Eigen::Vector3d;
using Matrix2 = Eigen::Matrix2d;
using Matrix3 = Eigen::Matrix3d;

using SparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::ColMajor, StorageIndex>;
using Triplet = Eigen::Triplet<Scalar, StorageIndex>;
using TripletList = std::vector<Triplet>;

/// Spatial dimension of the current formulation (2-D plane problems).
inline constexpr int kDim = 2;

/// Translational degrees of freedom carried by every node (u_x, u_y).
inline constexpr int kDofsPerNode = 2;

/// Number of independent in-plane stress/strain components in Voigt notation
/// (sigma_xx, sigma_yy, sigma_xy).
inline constexpr int kVoigt = 3;

/// Cartesian component labels used when parsing configuration files.
enum class Component : int { X = 0, Y = 1 };

/// Two-dimensional idealisation used to build the constitutive matrix.
enum class StressState {
  PlaneStress,  ///< sigma_zz = 0 (thin sheet). Default and fully verified path.
  PlaneStrain   ///< epsilon_zz = 0 (thick section / plane of symmetry).
};

}  // namespace sparlab
