/// \file test_material.cpp
/// \brief Constitutive matrices and material input validation.
#include "sparlab/core/Exceptions.hpp"
#include "sparlab/fem/StressRecovery.hpp"
#include "sparlab/material/IsotropicMaterial.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <Eigen/Eigenvalues>

using namespace sparlab;
using Catch::Approx;

TEST_CASE("plane-stress constitutive matrix matches the closed form", "[material]") {
  const Scalar e = 70.0e9;
  const Scalar nu = 0.33;
  const IsotropicMaterial m(e, nu, 2700.0);
  const Matrix3 d = m.plane_stress_matrix();

  const Scalar f = e / (1.0 - nu * nu);
  REQUIRE(d(0, 0) == Approx(f));
  REQUIRE(d(1, 1) == Approx(f));
  REQUIRE(d(0, 1) == Approx(f * nu));
  REQUIRE(d(1, 0) == Approx(f * nu));
  REQUIRE(d(2, 2) == Approx(f * (1.0 - nu) / 2.0));
  REQUIRE(d(0, 2) == Approx(0.0));
  REQUIRE(d(1, 2) == Approx(0.0));

  // The shear entry must equal the shear modulus: G = E / (2(1+nu)).
  REQUIRE(d(2, 2) == Approx(m.shear_modulus()));

  // Symmetric and positive definite.
  REQUIRE((d - d.transpose()).cwiseAbs().maxCoeff() == Approx(0.0).margin(1.0e-6));
  Eigen::SelfAdjointEigenSolver<Matrix3> es(d);
  REQUIRE(es.eigenvalues().minCoeff() > 0.0);
}

TEST_CASE("plane-strain constitutive matrix matches the closed form", "[material]") {
  const Scalar e = 200.0e9;
  const Scalar nu = 0.29;
  const IsotropicMaterial m(e, nu, 7850.0);
  const Matrix3 d = m.plane_strain_matrix();

  const Scalar f = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
  REQUIRE(d(0, 0) == Approx(f * (1.0 - nu)));
  REQUIRE(d(0, 1) == Approx(f * nu));
  REQUIRE(d(2, 2) == Approx(f * (1.0 - 2.0 * nu) / 2.0));
  REQUIRE(d(2, 2) == Approx(m.shear_modulus()));

  Eigen::SelfAdjointEigenSolver<Matrix3> es(d);
  REQUIRE(es.eigenvalues().minCoeff() > 0.0);

  // Plane strain is stiffer than plane stress in the normal directions.
  REQUIRE(d(0, 0) > m.plane_stress_matrix()(0, 0));
}

TEST_CASE("constitutive() dispatches on the stress state", "[material]") {
  const IsotropicMaterial m(1.0e9, 0.25, 1000.0);
  REQUIRE((m.constitutive(StressState::PlaneStress) - m.plane_stress_matrix())
              .cwiseAbs()
              .maxCoeff() == Approx(0.0));
  REQUIRE((m.constitutive(StressState::PlaneStrain) - m.plane_strain_matrix())
              .cwiseAbs()
              .maxCoeff() == Approx(0.0));
}

TEST_CASE("invalid material parameters are rejected with actionable errors",
          "[material][diagnostics]") {
  REQUIRE_THROWS_AS(IsotropicMaterial(0.0, 0.3, 1.0), ConfigError);
  REQUIRE_THROWS_AS(IsotropicMaterial(-1.0, 0.3, 1.0), ConfigError);
  REQUIRE_THROWS_AS(IsotropicMaterial(1.0, 0.5, 1.0), ConfigError);
  REQUIRE_THROWS_AS(IsotropicMaterial(1.0, 0.6, 1.0), ConfigError);
  REQUIRE_THROWS_AS(IsotropicMaterial(1.0, -1.5, 1.0), ConfigError);
  REQUIRE_THROWS_AS(IsotropicMaterial(1.0, 0.3, -1.0), ConfigError);
}

TEST_CASE("von Mises stress uses the correct out-of-plane component",
          "[material][stress]") {
  const Scalar nu = 0.3;

  SECTION("uniaxial tension gives the axial stress") {
    const Vector3 s(100.0e6, 0.0, 0.0);
    REQUIRE(von_mises(s, StressState::PlaneStress, nu) == Approx(100.0e6));
  }

  SECTION("pure shear gives sqrt(3) tau") {
    const Vector3 s(0.0, 0.0, 50.0e6);
    REQUIRE(von_mises(s, StressState::PlaneStress, nu) ==
            Approx(std::sqrt(3.0) * 50.0e6));
  }

  SECTION("equal biaxial tension gives the same value in plane stress") {
    const Vector3 s(80.0e6, 80.0e6, 0.0);
    REQUIRE(von_mises(s, StressState::PlaneStress, nu) == Approx(80.0e6));
  }

  SECTION("plane strain includes sigma_zz = nu (sxx + syy)") {
    const Vector3 s(80.0e6, 80.0e6, 0.0);
    // sigma_zz = 0.3 * 160e6 = 48e6, so the deviator shrinks.
    const Scalar sz = nu * (s(0) + s(1));
    const Scalar expected = std::abs(s(0) - sz);
    REQUIRE(von_mises(s, StressState::PlaneStrain, nu) == Approx(expected));
    REQUIRE(von_mises(s, StressState::PlaneStrain, nu) <
            von_mises(s, StressState::PlaneStress, nu));
  }
}

TEST_CASE("principal stresses bracket the normal components", "[stress]") {
  const Vector3 s(120.0e6, -40.0e6, 30.0e6);
  Scalar s1 = 0.0;
  Scalar s2 = 0.0;
  principal_stresses(s, s1, s2);
  REQUIRE(s1 >= s2);
  // Invariants: trace and determinant of the 2-D stress tensor.
  REQUIRE(s1 + s2 == Approx(s(0) + s(1)));
  REQUIRE(s1 * s2 == Approx(s(0) * s(1) - s(2) * s(2)));
}
