# SparLab convenience targets.
#
# Everything the README promises is one `make` target. Each target is a thin
# wrapper over a script in scripts/ or a direct binary invocation, so nothing is
# hidden and every command can also be typed by hand.
#
#   make configure   configure the CMake build (Release by default)
#   make build       build the library, apps and tests
#   make test        run the whole Catch2 suite through ctest
#   make verify      run the verification / validation studies
#   make cross-validation  compare nodal displacements (and buckling load
#                    factors) with CalculiX and scikit-fem
#   make tet10-study Tet4 against Tet10 on the engine mount meshed from CAD
#                    (needs the gmsh Python package)
#   make benchmark   run one benchmark          (CASE=cantilever_beam)
#   make benchmarks  run every benchmark deck (2-D and 3-D)
#   make study       run the aerospace parametric design study
#   make scaling     run the runtime scaling benchmarks (2-D and 3-D)
#   make figures     regenerate every figure and the animation
#   make meshes      regenerate the Gmsh meshes of the real-geometry decks
#                    (needs the gmsh Python package; the meshes are committed)
#   make results     refresh the result tables committed under docs/results
#   make all         build, test, verify, cross-validation, tet10-study,
#                    benchmarks, study, scaling, figures, results
#   make clean       remove the build directory
#   make distclean   also remove results/

BUILD_DIR   ?= build
BUILD_TYPE  ?= Release
JOBS        ?= $(shell nproc 2>/dev/null || echo 4)
GENERATOR   ?= $(shell command -v ninja >/dev/null 2>&1 && echo Ninja || echo "Unix Makefiles")
CASE        ?= cantilever_beam
RESULTS_DIR ?= results
PYTHON      ?= python3

CMAKE_FLAGS ?= -DCMAKE_BUILD_TYPE=$(BUILD_TYPE)

.PHONY: all configure build test verify cross-validation tet10-study benchmark \
        benchmarks study scaling figures meshes results clean distclean help

help:
	@sed -n '3,26p' Makefile

configure:
	cmake -S . -B $(BUILD_DIR) -G "$(GENERATOR)" $(CMAKE_FLAGS)

build: configure
	cmake --build $(BUILD_DIR) -j $(JOBS)

test: build
	cd $(BUILD_DIR) && ctest --output-on-failure -j $(JOBS)

verify: build
	./scripts/run_verification.sh $(RESULTS_DIR)/verification

cross-validation: build
	./scripts/run_cross_validation.sh $(RESULTS_DIR)

tet10-study: build
	./scripts/run_tet10_study.sh $(RESULTS_DIR)

benchmark: build
	./scripts/run_benchmark.sh $(CASE) $(RESULTS_DIR)

benchmarks: build
	./scripts/run_all_benchmarks.sh $(RESULTS_DIR)

study: build
	./scripts/run_aerospace_study.sh $(RESULTS_DIR)/study

scaling: build
	./scripts/run_scaling.sh $(RESULTS_DIR)/benchmark

figures:
	./scripts/make_figures.sh $(RESULTS_DIR) docs/figures

meshes:
	$(PYTHON) python/scripts/make_meshes.py --output configs/meshes --tet10

results:
	$(PYTHON) python/scripts/write_result_tables.py --results $(RESULTS_DIR) \
	    --output docs/results

all: build test verify cross-validation tet10-study benchmarks study scaling figures \
     results

clean:
	rm -rf $(BUILD_DIR)

distclean: clean
	rm -rf $(RESULTS_DIR)
