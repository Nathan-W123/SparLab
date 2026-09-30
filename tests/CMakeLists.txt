# Catch2 v3 is preferred; the system package is used when available and
# otherwise fetched, so a fresh clone can build its tests offline *or* online.
find_package(Catch2 3 QUIET)

if(NOT Catch2_FOUND)
  message(STATUS "SparLab: Catch2 3 not found locally, fetching it")
  include(FetchContent)
  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.5.2
  )
  FetchContent_MakeAvailable(Catch2)
endif()

set(SPARLAB_TEST_SOURCES
  test_quadrature_shape.cpp
  test_material.cpp
  test_mesh.cpp
  test_element.cpp
  test_assembly_solver.cpp
  test_static_verification.cpp
  test_modal.cpp
  test_topopt.cpp
  test_io.cpp
  test_solid.cpp
  test_geometry.cpp
  test_mma.cpp
  test_unstructured.cpp
  test_multigrid.cpp
  test_projection.cpp
  test_tet10.cpp
  test_buckling.cpp
  test_manufacturing.cpp
  test_dofs.cpp
  test_loads_thermal.cpp
  test_design_loads.cpp
  test_shell_topology.cpp
  test_nonlinear_part_check.cpp
  test_nonlinear.cpp
  test_plasticity.cpp
  test_dynamics.cpp
  test_contact.cpp
  test_shell.cpp
  test_beam.cpp
)

add_executable(sparlab_tests ${SPARLAB_TEST_SOURCES})
target_link_libraries(sparlab_tests PRIVATE sparlab Catch2::Catch2WithMain)
target_include_directories(sparlab_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(SPARLAB_ENABLE_WARNINGS AND NOT MSVC)
  target_compile_options(sparlab_tests PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(Catch)
catch_discover_tests(sparlab_tests
  TEST_PREFIX "sparlab."
  PROPERTIES TIMEOUT 900
)

# The whole suite once more in a single process, which catches state one test
# case leaves behind for the next; the cases above each run in their own.
# Debug builds leave it out: they already run every case under Eigen's
# assertions, and the suite in one process there outlasts its timeout (over
# half an hour on a CI runner) for no coverage the Release run lacks.
string(TOUPPER "${CMAKE_BUILD_TYPE}" SPARLAB_TEST_BUILD_TYPE)
if(NOT SPARLAB_TEST_BUILD_TYPE STREQUAL "DEBUG")
  add_test(NAME sparlab.all COMMAND sparlab_tests)
  set_tests_properties(sparlab.all PROPERTIES TIMEOUT 1800)
endif()
