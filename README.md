# Poisson & Heat on Irregular Domains (Cartesian Grids)

A C++20 codebase for solving Poisson and Heat equations on irregular 2D domains using uniform Cartesian grids and level-set geometry. The design emphasizes clear separation of concerns: data structures, finite-difference operators, domain/level-set management, boundary handling, linear algebra, and time integration.

---

## Overview

- **Grid2D**: uniform node-centered grids with inclusive endpoints.
- **Array2D**: lightweight 2D container with `A(i,j)` indexing and basic linear-algebra style operators.
- **FiniteDifference2D**: interior finite-difference operators (gradients, Laplacian, upwind building blocks) and level-set geometry helpers (e.g., curvature).
- **Irregular boundaries**: handled via level sets, ghost values, and cut-cell stencils (planned).
- **Roadmap**: Poisson → Heat (explicit/implicit) → Two-phase Stefan with a moving interface.



