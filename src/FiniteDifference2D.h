#ifndef FINITEDIFFERENCE2D_H
#define FINITEDIFFERENCE2D_H

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include "Array2D.h"
#include "Grid2D.h"

/**
 * @file FiniteDifference2D.h
 * @brief Pure mathematical finite difference operators for 2D grids
 * 
 * This class provides ONLY the mathematical discretization operators.
 * It knows nothing about boundary conditions, ghost cells, or cut cells.
 * Higher layers will modify these operators for specific problems.
 * 
 * Design Philosophy:
 * - Returns INVALID values (NaN or 0) at boundaries
 * - Caller is responsible for boundary treatment
 * - All operators work on interior points only
 * - No ghost cell handling here
 */

template<typename T>
class FiniteDifference2D {
    static_assert(std::is_floating_point<T>::value, 
                  "FiniteDifference2D requires floating-point type");

public:
    // ============== BASIC DIFFERENCE OPERATORS ==============
    
    /**
     * @brief Forward difference in X: (u[i+1,j] - u[i,j])/dx
     * @note Returns 0 at right boundary (i = nx-1)
     */
    static Array2D<T> forward_diff_x(const Array2D<T>& u, const Grid2D<T>& grid) {
        const auto nx = grid.nx();
        const auto ny = grid.ny();
        const T inv_dx = T{1} / grid.dx();
        Array2D<T> du(nx, ny, T{0});
        
        // Interior points only
        for (size_t i = 0; i < nx-1; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                du(i,j) = (u(i+1,j) - u(i,j)) * inv_dx;
            }
        }
        // Right boundary left as 0 - caller must handle
        return du;
    }
    
    /**
     * @brief Backward difference in X: (u[i,j] - u[i-1,j])/dx
     * @note Returns 0 at left boundary (i = 0)
     */
    static Array2D<T> backward_diff_x(const Array2D<T>& u, const Grid2D<T>& grid) {
        const auto nx = grid.nx();
        const auto ny = grid.ny();
        const T inv_dx = T{1} / grid.dx();
        Array2D<T> du(nx, ny, T{0});
        
        for (size_t i = 1; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                du(i,j) = (u(i,j) - u(i-1,j)) * inv_dx;
            }
        }
        return du;
    }
    
    /**
     * @brief Central difference in X: (u[i+1,j] - u[i-1,j])/(2*dx)
     * @note Returns 0 at boundaries (i = 0 and i = nx-1)
     */
    static Array2D<T> central_diff_x(const Array2D<T>& u, const Grid2D<T>& grid) {
        const auto nx = grid.nx();
        const auto ny = grid.ny();
        const T inv_2dx = T{0.5} / grid.dx();
        Array2D<T> du(nx, ny, T{0});
        
        // Interior points only
        for (size_t i = 1; i < nx-1; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                du(i,j) = (u(i+1,j) - u(i-1,j)) * inv_2dx;
            }
        }
        return du;
    }
    
    // Y-direction operators (similar structure)
    static Array2D<T> forward_diff_y(const Array2D<T>& u, const Grid2D<T>& grid) {
        const auto nx = grid.nx();
        const auto ny = grid.ny();
        const T inv_dy = T{1} / grid.dy();
        Array2D<T> du(nx, ny, T{0});
        
        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny-1; ++j) {
                du(i,j) = (u(i,j+1) - u(i,j)) * inv_dy;
            }
        }
        return du;
    }
    
    static Array2D<T> backward_diff_y(const Array2D<T>& u, const Grid2D<T>& grid) {
        const auto nx = grid.nx();
        const auto ny = grid.ny();
        const T inv_dy = T{1} / grid.dy();
        Array2D<T> du(nx, ny, T{0});
        
        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 1; j < ny; ++j) {
                du(i,j) = (u(i,j) - u(i,j-1)) * inv_dy;
            }
        }
        return du;
    }
    
    static Array2D<T> central_diff_y(const Array2D<T>& u, const Grid2D<T>& grid) {
        const auto nx = grid.nx();
        const auto ny = grid.ny();
        const T inv_2dy = T{0.5} / grid.dy();
        Array2D<T> du(nx, ny, T{0});
        
        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 1; j < ny-1; ++j) {
                du(i,j) = (u(i,j+1) - u(i,j-1)) * inv_2dy;
            }
        }
        return du;
    }
    
    // ============== GRADIENT AND DIVERGENCE ==============
    
    /**
     * @brief Compute gradient ∇u = (∂u/∂x, ∂u/∂y) using central differences
     * @return Pair of arrays {ux, uy}
     */
    static std::pair<Array2D<T>, Array2D<T>> gradient(const Array2D<T>& u, 
                                                       const Grid2D<T>& grid) {
        return {central_diff_x(u, grid), central_diff_y(u, grid)};
    }
    
    /**
     * @brief Compute divergence ∇·F = ∂Fx/∂x + ∂Fy/∂y
     */
    static Array2D<T> divergence(const Array2D<T>& Fx, const Array2D<T>& Fy,
                                 const Grid2D<T>& grid) {
        auto div_x = central_diff_x(Fx, grid);
        auto div_y = central_diff_y(Fy, grid);
        return div_x + div_y;
    }
    
    // ============== LAPLACIAN OPERATORS ==============
    
    /**
     * @brief 5-point Laplacian for INTERIOR points only
     * @note Returns 0 at ALL boundary points
     * @note Caller must apply boundary conditions
     */
    static Array2D<T> laplacian5_interior(const Array2D<T>& u, const Grid2D<T>& grid) {
        const auto nx = grid.nx();
        const auto ny = grid.ny();
        const T inv_dx2 = T{1} / (grid.dx() * grid.dx());
        const T inv_dy2 = T{1} / (grid.dy() * grid.dy());
        
        Array2D<T> lap(nx, ny, T{0});
        
        // ONLY interior points - no boundary handling
        for (size_t i = 1; i < nx-1; ++i) {
            for (size_t j = 1; j < ny-1; ++j) {
                lap(i,j) = (u(i+1,j) - T{2}*u(i,j) + u(i-1,j)) * inv_dx2 +
                          (u(i,j+1) - T{2}*u(i,j) + u(i,j-1)) * inv_dy2;
            }
        }
        
        return lap;
    }
    
    /**
     * @brief 9-point Laplacian for improved rotational invariance
     * @note Interior points only
     */
    static Array2D<T> laplacian9_interior(const Array2D<T>& u, const Grid2D<T>& grid) {
        const auto nx = grid.nx();
        const auto ny = grid.ny();
        const T dx = grid.dx();
        const T dy = grid.dy();
        const T h2 = dx * dy;
        
        Array2D<T> lap(nx, ny, T{0});
        
        // Weights for 9-point stencil (Mehrstellen)
        const T a = T{2} / (T{3} * dx * dx);
        const T b = T{2} / (T{3} * dy * dy);
        const T c = T{1} / (T{6} * h2);
        const T d = -T{2} * (a + b);
        
        for (size_t i = 1; i < nx-1; ++i) {
            for (size_t j = 1; j < ny-1; ++j) {
                lap(i,j) = d * u(i,j) +
                          a * (u(i+1,j) + u(i-1,j)) +
                          b * (u(i,j+1) + u(i,j-1)) +
                          c * (u(i+1,j+1) + u(i-1,j+1) + 
                               u(i+1,j-1) + u(i-1,j-1));
            }
        }
        
        return lap;
    }
    
    // ============== LEVEL SET SPECIFIC OPERATORS ==============
    
    /**
     * @brief Compute |∇u| using upwind differences (for Hamilton-Jacobi)
     * @note Essential for level set reinitialization
     */
    static Array2D<T> gradient_magnitude_upwind(const Array2D<T>& u,
                                                const Grid2D<T>& grid) {
        const auto nx = grid.nx();
        const auto ny = grid.ny();
        const T dx = grid.dx();
        const T dy = grid.dy();
        
        Array2D<T> mag(nx, ny, T{0});
        
        for (size_t i = 1; i < nx-1; ++i) {
            for (size_t j = 1; j < ny-1; ++j) {
                // Godunov's scheme for |∇u|
                T ux_plus = (u(i+1,j) - u(i,j)) / dx;
                T ux_minus = (u(i,j) - u(i-1,j)) / dx;
                T uy_plus = (u(i,j+1) - u(i,j)) / dy;
                T uy_minus = (u(i,j) - u(i,j-1)) / dy;
                
                T a_plus = std::max(ux_minus, T{0});
                T b_minus = std::min(ux_plus, T{0});
                T c_plus = std::max(uy_minus, T{0});
                T d_minus = std::min(uy_plus, T{0});
                
                T grad_x_sq = std::max(a_plus*a_plus, b_minus*b_minus);
                T grad_y_sq = std::max(c_plus*c_plus, d_minus*d_minus);
                
                mag(i,j) = std::sqrt(grad_x_sq + grad_y_sq);
            }
        }
        
        return mag;
    }
    
    /**
     * @brief Compute mean curvature κ = ∇·(∇φ/|∇φ|)
     * @note For level set regularization and surface tension
     */
    static Array2D<T> mean_curvature(const Array2D<T>& phi,
                                     const Grid2D<T>& grid,
                                     T epsilon = T{1e-6}) {
        const auto nx = grid.nx();
        const auto ny = grid.ny();
        const T dx = grid.dx();
        const T dy = grid.dy();
        
        Array2D<T> kappa(nx, ny, T{0});
        
        for (size_t i = 1; i < nx-1; ++i) {
            for (size_t j = 1; j < ny-1; ++j) {
                // Central differences for derivatives
                T phi_x = (phi(i+1,j) - phi(i-1,j)) / (T{2} * dx);
                T phi_y = (phi(i,j+1) - phi(i,j-1)) / (T{2} * dy);
                T phi_xx = (phi(i+1,j) - T{2}*phi(i,j) + phi(i-1,j)) / (dx * dx);
                T phi_yy = (phi(i,j+1) - T{2}*phi(i,j) + phi(i,j-1)) / (dy * dy);
                T phi_xy = (phi(i+1,j+1) - phi(i-1,j+1) - 
                           phi(i+1,j-1) + phi(i-1,j-1)) / (T{4} * dx * dy);
                
                T denom = std::pow(phi_x*phi_x + phi_y*phi_y + epsilon, T{1.5});
                
                kappa(i,j) = (phi_xx * (phi_y*phi_y + epsilon) + 
                             phi_yy * (phi_x*phi_x + epsilon) - 
                             T{2} * phi_x * phi_y * phi_xy) / denom;
            }
        }
        
        return kappa;
    }
    
    /**
     * @brief Sign function for level set reinitialization
     */
    static Array2D<T> sign_function(const Array2D<T>& phi, 
                                    const Grid2D<T>& grid) {
        const auto nx = grid.nx();
        const auto ny = grid.ny();
        const T epsilon = std::min(grid.dx(), grid.dy());
        
        Array2D<T> S(nx, ny);
        
        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                S(i,j) = phi(i,j) / std::sqrt(phi(i,j)*phi(i,j) + epsilon*epsilon);
            }
        }
        
        return S;
    }
    
    // ============== UPWIND SCHEMES ==============
    
    /**
     * @brief First-order upwind for advection ∂u/∂t + c·∇u = 0
     * @param u Field to differentiate
     * @param cx X-component of velocity
     * @param cy Y-component of velocity
     * @return c·∇u using upwind differencing
     */
    static Array2D<T> upwind_advection(const Array2D<T>& u,
                                       const Array2D<T>& cx,
                                       const Array2D<T>& cy,
                                       const Grid2D<T>& grid) {
        const auto nx = grid.nx();
        const auto ny = grid.ny();
        const T dx = grid.dx();
        const T dy = grid.dy();
        
        Array2D<T> advect(nx, ny, T{0});
        
        for (size_t i = 1; i < nx-1; ++i) {
            for (size_t j = 1; j < ny-1; ++j) {
                T ux, uy;
                
                // Upwind in x-direction
                if (cx(i,j) > 0) {
                    ux = (u(i,j) - u(i-1,j)) / dx;  // Backward
                } else if (cx(i,j) < 0) {
                    ux = (u(i+1,j) - u(i,j)) / dx;  // Forward
                } else {
                    ux = (u(i+1,j) - u(i-1,j)) / (T{2} * dx);  // Central
                }
                
                // Upwind in y-direction
                if (cy(i,j) > 0) {
                    uy = (u(i,j) - u(i,j-1)) / dy;  // Backward
                } else if (cy(i,j) < 0) {
                    uy = (u(i,j+1) - u(i,j)) / dy;  // Forward
                } else {
                    uy = (u(i,j+1) - u(i,j-1)) / (T{2} * dy);  // Central
                }
                
                advect(i,j) = cx(i,j) * ux + cy(i,j) * uy;
            }
        }
        
        return advect;
    }
    
    // ============== HIGH-RESOLUTION SCHEMES ==============
    
    /**
     * @brief ENO2 (2nd-order Essentially Non-Oscillatory) in X direction
     * @note Requires at least 2 ghost cells on each side
     * @note Returns 0 near boundaries
     */
    static Array2D<T> eno2_x(const Array2D<T>& u,
                             const Array2D<T>& velocity,
                             const Grid2D<T>& grid) {
        const auto nx = grid.nx();
        const auto ny = grid.ny();
        const T dx = grid.dx();
        
        Array2D<T> du(nx, ny, T{0});
        
        // Need at least 4 points for ENO2
        for (size_t i = 2; i < nx-2; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                if (velocity(i,j) > 0) {
                    // Left-biased stencil
                    T D0 = (u(i,j) - u(i-1,j)) / dx;
                    T D1 = (u(i-1,j) - u(i-2,j)) / dx;
                    T D2 = (u(i+1,j) - u(i,j)) / dx;
                    
                    // Choose smoother stencil
                    T DD0 = std::abs(D2 - D0);
                    T DD1 = std::abs(D0 - D1);
                    
                    if (DD1 < DD0) {
                        // Use u[i-2], u[i-1], u[i]
                        du(i,j) = D0 + T{0.5} * (D0 - D1);
                    } else {
                        // Use u[i-1], u[i], u[i+1]
                        du(i,j) = D0 + T{0.5} * (D2 - D0);
                    }
                } else if (velocity(i,j) < 0) {
                    // Right-biased stencil (similar logic)
                    T D0 = (u(i+1,j) - u(i,j)) / dx;
                    T D1 = (u(i+2,j) - u(i+1,j)) / dx;
                    T D2 = (u(i,j) - u(i-1,j)) / dx;
                    
                    T DD0 = std::abs(D0 - D2);
                    T DD1 = std::abs(D1 - D0);
                    
                    if (DD1 < DD0) {
                        du(i,j) = D0 + T{0.5} * (D0 - D1);
                    } else {
                        du(i,j) = D0 + T{0.5} * (D2 - D0);
                    }
                } else {
                    // Central difference for zero velocity
                    du(i,j) = (u(i+1,j) - u(i-1,j)) / (T{2} * dx);
                }
            }
        }
        
        return du;
    }
    
    // ============== STENCIL STRUCTURE (for matrix assembly) ==============
    
    /**
     * @brief Structure representing a finite difference stencil
     * @note Used by higher layers to build matrices
     */
    struct Stencil {
        T center;      // Coefficient for u[i,j]
        T east;        // Coefficient for u[i+1,j]
        T west;        // Coefficient for u[i-1,j]
        T north;       // Coefficient for u[i,j+1]
        T south;       // Coefficient for u[i,j-1]
        T rhs_contrib; // Contribution to right-hand side
        
        Stencil() : center(0), east(0), west(0), north(0), south(0), rhs_contrib(0) {}
    };
    
    /**
     * @brief Get Laplacian stencil at point (i,j)
     * @note Returns standard 5-point stencil for interior points
     * @note Returns zero stencil for boundary points
     */
    static Stencil get_laplacian_stencil(size_t i, size_t j,
                                         const Grid2D<T>& grid) {
        Stencil s;
        
        const auto nx = grid.nx();
        const auto ny = grid.ny();
        
        // Check if interior point
        if (i > 0 && i < nx-1 && j > 0 && j < ny-1) {
            const T inv_dx2 = T{1} / (grid.dx() * grid.dx());
            const T inv_dy2 = T{1} / (grid.dy() * grid.dy());
            
            s.center = -T{2} * (inv_dx2 + inv_dy2);
            s.east = inv_dx2;
            s.west = inv_dx2;
            s.north = inv_dy2;
            s.south = inv_dy2;
        }
        // Boundary points return zero stencil - caller must handle
        
        return s;
    }
    
    /**
     * @brief Check if point is on computational boundary
     */
    static bool is_boundary_point(size_t i, size_t j,
                                  const Grid2D<T>& grid) {
        return (i == 0 || i == grid.nx()-1 || 
                j == 0 || j == grid.ny()-1);
    }
};

#endif // FINITEDIFFERENCE2D_H