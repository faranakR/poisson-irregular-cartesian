# Poisson & Heat on Irregular Domains (Cartesian Grids)

A C++20 codebase for solving Poisson and Heat equations on irregular 2D domains using uniform Cartesian grids and level-set geometry. The design emphasizes clear separation of concerns: data structures, finite-difference operators, domain/level-set management, boundary handling, linear algebra, and time integration.

---

## Overview

- **Grid2D**: uniform node-centered grids with inclusive endpoints.
- **Array2D**: lightweight 2D container with `A(i,j)` indexing and basic linear-algebra style operators.
- **FiniteDifference2D**: interior finite-difference operators (gradients, Laplacian, upwind building blocks) and level-set geometry helpers (e.g., curvature).
- **Irregular boundaries**: handled via level sets, ghost values, and cut-cell stencils (planned).
- **Roadmap**: Poisson → Heat (explicit/implicit) → Two-phase Stefan with a moving interface.

---

## Repository Layout
.
├── src/
│ ├── Array2D.h # Dense 2D array (row-major), operators, indexing
│ ├── Grid2D.h # Uniform Cartesian grid, meshgrid, utilities
│ └── FiniteDifference2D.h # Interior-only FD ops, level-set helpers
├── tests/
│ ├── test_array2d.cpp
│ ├── test_grid2d.cpp
│ └── test_finitedifference2d.cpp
├── CMakeLists.txt
├── main.cpp # Architecture notes / scratch
└── README.md


---

## Architecture (layers)

1. **Core Data**  
   `Array2D`, `Grid2D`
2. **Mathematical Operators**  
   `FiniteDifference2D` (interior-only: no boundary conditions, no ghost cells)
3. **Domain Management** *(planned)*  
   `LevelSet2D` (advection, reinitialization), `DomainDecomposition2D`
4. **Boundary Conditions** *(planned)*  
   Outer BCs (Dirichlet/Neumann/Robin/Periodic), interface BCs (φ=0), ghost/cut-cells
5. **Linear Algebra** *(planned)*  
   Sparse assembly, PCG/GMRES, preconditioning
6. **Time Integration** *(planned)*  
   Explicit, implicit, Crank–Nicolson
7. **Problem Solvers** *(planned)*  
   Poisson, Heat, Stefan (two-phase, moving interface)

> Intentional design choice: `FiniteDifference2D` does **not** enforce boundary conditions. That responsibility lives above, so interior operators remain reusable across Poisson/Heat/Stefan.

---

## Build

Requirements: CMake ≥ 3.16, a C++20 compiler.

```bash
cmake -S . -B build
cmake --build build -j


Run tests via CTest:

ctest --test-dir build --output-on-failure


Or run individual executables (names may vary by generator):

./build/test_array2d
./build/test_grid2d
./build/test_finite_difference2d


CLion users: select targets in the Run/Debug configuration dropdown.

What’s Implemented and Tested

Array2D
Construction, indexing, shape queries, elementwise arithmetic, safety checks.

Grid2D
Axis generation, meshgrid() (MATLAB/NumPy convention), utilities (coords, nearest, contains).

FiniteDifference2D
Forward/backward/central differences (interior points), gradient, 5-point Laplacian, upwind building blocks, sign function, mean curvature (discrete approximation).

Test coverage includes exact polynomials, harmonic functions, and geometry checks on manufactured data.

Usage Example
#include "src/Array2D.h"
#include "src/Grid2D.h"
#include "src/FiniteDifference2D.h"
#include <iostream>

int main() {
    // Grid: 51x51 on [0,1]×[0,1]
    Grid2D<double> grid(/*nx=*/51, /*ny=*/51, /*xmin=*/0.0, /*xmax=*/1.0,
                        /*ymin=*/0.0, /*ymax=*/1.0);

    // Field u(x,y) = x^2 + y^2
    Array2D<double> u(grid.nx(), grid.ny(), 0.0);
    for (std::size_t i = 0; i < grid.nx(); ++i)
        for (std::size_t j = 0; j < grid.ny(); ++j)
            u(i,j) = grid.x(i)*grid.x(i) + grid.y(j)*grid.y(j);

    // Interior Laplacian (boundary entries are zero by design here)
    auto lap = FiniteDifference2D<double>::laplacian5_interior(u, grid);

    // Inspect an interior point (avoid boundaries)
    std::size_t I = 25, J = 25;
    std::cout << "Laplacian at (" << I << "," << J << ") ≈ " << lap(I,J) << "\n";
    // Expected: ∇²u = 2 + 2 = 4
}

Roadmap

Level set layer (LevelSet2D): advection, reinitialization, normals/curvature utilities

Boundary manager for:

Outer BCs (Dirichlet, Neumann, Robin, Periodic)

Interface BCs on φ=0 via ghost and cut-cell methods

Symmetric cut-cell discretizations suitable for PCG

Heat equation (explicit / implicit / Crank–Nicolson)

Stefan problem (two-phase) with interface motion from jump conditions

Notes

Operators in FiniteDifference2D act on interior points. Boundary rows/cols are intentionally left for higher-level boundary handling.

Level-set, ghost fluid, and cut-cell functionality will be introduced as separate modules to keep numerical operators reusable and verifiable in isolation.


If you want this split into sections with badges, references, or a minimal “Getting Started” tutorial later, we can add that once Layer 3 is in place.

You said:
give me all in one .md file?
ChatGPT said:
# Poisson & Heat on Irregular Domains (Cartesian Grids)

A C++20 codebase for solving Poisson and Heat equations on **irregular 2D domains** using **uniform Cartesian grids** and **level-set** geometry. The design emphasizes a clear separation of concerns: data structures, finite-difference operators, domain/level-set management, boundary handling, linear algebra, and time integration.

---

## Overview

- **Grid2D** — uniform node-centered grids with inclusive endpoints.
- **Array2D** — lightweight 2D container with `A(i,j)` indexing and basic linear-algebra-style operators.
- **FiniteDifference2D** — interior finite-difference operators (gradients, Laplacian, upwind building blocks) and level-set geometry helpers (e.g., curvature).
- **Irregular boundaries** — handled via level sets, ghost values, and cut-cell stencils (planned).
- **Roadmap** — Poisson → Heat (explicit/implicit) → Two-phase Stefan with a moving interface.

---

## Repository Layout



.
├── src/
│ ├── Array2D.h # Dense 2D array (row-major), operators, indexing
│ ├── Grid2D.h # Uniform Cartesian grid, meshgrid, utilities
│ └── FiniteDifference2D.h # Interior-only FD ops, level-set helpers
├── tests/
│ ├── test_array2d.cpp
│ ├── test_grid2d.cpp
│ └── test_finitedifference2d.cpp
├── CMakeLists.txt
├── main.cpp # Architecture notes / scratch
└── README.md


---

## Architecture (Layers)

1. **Core Data**  
   `Array2D`, `Grid2D`
2. **Mathematical Operators**  
   `FiniteDifference2D` (interior-only: no boundary conditions, no ghost cells)
3. **Domain Management** *(planned)*  
   `LevelSet2D` (advection, reinitialization), `DomainDecomposition2D`
4. **Boundary Conditions** *(planned)*  
   Outer BCs (Dirichlet/Neumann/Robin/Periodic), interface BCs (φ=0), ghost/cut-cells
5. **Linear Algebra** *(planned)*  
   Sparse assembly, PCG/GMRES, preconditioning
6. **Time Integration** *(planned)*  
   Explicit, implicit, Crank–Nicolson
7. **Problem Solvers** *(planned)*  
   Poisson, Heat, Stefan (two-phase, moving interface)

> **Design choice:** `FiniteDifference2D` does **not** enforce boundary conditions. That responsibility lives above, so interior operators remain reusable across Poisson/Heat/Stefan.

---

## Build

**Requirements:** CMake ≥ 3.16, a C++20 compiler.

```bash
cmake -S . -B build
cmake --build build -j


Run tests via CTest:

ctest --test-dir build --output-on-failure


Or run individual executables (names may vary by generator):

./build/test_array2d
./build/test_grid2d
./build/test_finite_difference2d


CLion users: select targets in the Run/Debug configuration dropdown.

What’s Implemented and Tested

Array2D
Construction, indexing, shape queries, elementwise arithmetic, safety checks.

Grid2D
Axis generation, meshgrid() (MATLAB/NumPy convention), utilities (coords, nearest, contains).

FiniteDifference2D
Forward/backward/central differences (interior points), gradient, 5-point Laplacian, upwind building blocks, sign function, mean curvature (discrete approximation).

Test coverage includes exact polynomials, harmonic functions, and geometry checks on manufactured data.

Usage Example
#include "src/Array2D.h"
#include "src/Grid2D.h"
#include "src/FiniteDifference2D.h"
#include <iostream>

int main() {
    // Grid: 51x51 on [0,1]×[0,1]
    Grid2D<double> grid(/*nx=*/51, /*ny=*/51,
                        /*xmin=*/0.0, /*xmax=*/1.0,
                        /*ymin=*/0.0, /*ymax=*/1.0);

    // Field u(x,y) = x^2 + y^2
    Array2D<double> u(grid.nx(), grid.ny(), 0.0);
    for (std::size_t i = 0; i < grid.nx(); ++i)
        for (std::size_t j = 0; j < grid.ny(); ++j)
            u(i,j) = grid.x(i)*grid.x(i) + grid.y(j)*grid.y(j);

    // Interior Laplacian (boundary entries are zero by design here)
    auto lap = FiniteDifference2D<double>::laplacian5_interior(u, grid);

    // Inspect an interior point (avoid boundaries)
    std::size_t I = 25, J = 25;
    std::cout << "Laplacian at (" << I << "," << J << ") ≈ " << lap(I,J) << "\n";
    // Expected: ∇²u = 2 + 2 = 4
}

Roadmap

Level set layer (LevelSet2D): advection, reinitialization, normals/curvature utilities.

Boundary manager for:

Outer BCs (Dirichlet, Neumann, Robin, Periodic)

Interface BCs on φ=0 via ghost and cut-cell methods

Symmetric cut-cell discretizations suitable for PCG.

Heat equation (explicit / implicit / Crank–Nicolson).

Stefan problem (two-phase) with interface motion from jump conditions.

Notes

Operators in FiniteDifference2D act on interior points. Boundary rows/cols are intentionally left for higher-level boundary handling.

Level-set, ghost fluid, and cut-cell functionality will be introduced as separate modules to keep numerical operators reusable and verifiable in isolation.
