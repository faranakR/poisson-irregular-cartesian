cmake_minimum_required(VERSION 3.16)
project(Poisson_Solver
        VERSION 1.0.0
        DESCRIPTION "Heat Equation Solver on Irregular Domains using Cartesian Grids"
        LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    message(STATUS "Build type not specified, defaulting to Release")
endif()

# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    # Common flags
    add_compile_options(-Wall -Wextra -Wpedantic)

    # Debug flags
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -DDEBUG")

    # Release flags
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG")

    # Release with debug info
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O2 -g")

elseif(MSVC)
    # MSVC flags
    add_compile_options(/W4)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Od /DEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
endif()

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/src)

# Source files - organize by component (for future additions)
set(CORE_HEADERS
        src/Array2D.h
        # Future headers will be added here:
        # src/Grid2D.h
        # src/Gradient2D.h
        # src/Laplace2D.h
        # src/Phi2D.h
        # src/Integration2D.h
        # src/Boundary2D.h
        # src/LinearSolver2D.h
        # src/Examples2D.h
)

set(CORE_SOURCES
        # Future source files will be added here if headers have implementations:
        # src/Grid2D.cpp
        # src/Phi2D.cpp
        # etc.
)

# Main executable
add_executable(heat_solver
        main.cpp
        ${CORE_HEADERS}
        ${CORE_SOURCES})

target_compile_features(heat_solver PRIVATE cxx_std_20)

# Test executable
add_executable(test_array2d
        test_array2d.cpp
        src/Array2D.h)

target_compile_features(test_array2d PRIVATE cxx_std_20)

# Optional: Add a custom target to run tests
add_custom_target(run_tests
        COMMAND test_array2d
        DEPENDS test_array2d
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running Array2D tests...")

# Optional: Enable testing with CTest
enable_testing()
add_test(NAME Array2D_Tests COMMAND test_array2d)

# Installation rules (optional)
install(TARGETS heat_solver
        RUNTIME DESTINATION bin)

install(FILES ${CORE_HEADERS}
        DESTINATION include/heat_equation_solver)

# Print configuration summary
message(STATUS "")
message(STATUS "=== Configuration Summary ===")
message(STATUS "Project:        ${PROJECT_NAME} v${PROJECT_VERSION}")
message(STATUS "Build type:     ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Standard:   ${CMAKE_CXX_STANDARD}")
message(STATUS "Compiler:       ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")
message(STATUS "Targets:")
message(STATUS "  heat_solver  - Main heat equation solver")
message(STATUS "  test_array2d - Array2D class test suite")
message(STATUS "  run_tests    - Run all tests")
message(STATUS "")
message(STATUS "To build: cmake --build .")
message(STATUS "To test:  cmake --build . --target run_tests")
message(STATUS "To install: cmake --install .")
message(STATUS "=============================")
message(STATUS "")