#ifndef GRID2D_H
#define GRID2D_H

#include <vector>
#include <algorithm>
#include <utility>
#include <type_traits>
#include <cstddef>
#include <stdexcept>
#include "Array2D.h"

/**
 * @brief Uniform Cartesian 2D grid with inclusive endpoints.
 *
 * Nodes: x[i] = xmin + i*dx,  i = 0..nx-1,  dx = (xmax-xmin)/(nx-1)
 *        y[j] = ymin + j*dy,  j = 0..ny-1,  dy = (ymax-ymin)/(ny-1)
 *
 * Stores axis coordinates as 1D arrays (no duplication).
 * You can materialize full meshgrid arrays via meshgridX()/meshgridY().
 */

template<typename T>
class Grid2D {
    static_assert(std::is_floating_point<T>::v, "Grid2D<T>: T must be floating-point type.");

public:
    using size_type = std::size_t;

    /** Construct a node-centered uniform grid with inclusive endpoints **/
    explicit Grid2D(T xmin = static_cast<T>(-1),
           T xmax = static_cast<T>( 1),
           T ymin = static_cast<T>(-1),
           T ymax = static_cast<T>( 1),
           size_type nx = 100,
           size_type ny = 100) :
    xmin_(xmin), xmax_(xmax), ymin_(ymin), ymax_(ymax),
    nx_(nx), ny_(ny), x_(nx), y_(ny)
    {
        if (!(xmax_ > xmin_)) throw std::invalid_argument("Grid2D: xmax must be greater than xmin.");
        if (!(ymax_ > ymin_)) throw std::invalid_argument("Grid2D: ymax must be greater than ymin.");
        if (nx_ < 2 || ny_ < 2) throw std::invalid_argument("Grid2D: nx_ or ny_ must be greater than 2.");

        const T dxv = dx();
        const T dyv = dy();
        for (size_type i = 0; i < nx_; ++i) {x_[i] = xmin_ + dxv * static_cast<T>(i);}
        for (size_type j = 0; j < ny_; ++j) {y_[j] = ymin_ + dyv * static_cast<T>(j);}
    }

    [[nodiscard]] size_type nx()   const noexcept {return nx_;}
    [[nodiscard]] size_type ny()   const noexcept {return ny_;}
    [[nodiscard]] T         xmin() const noexcept {return xmin_;}
    [[nodiscard]] T         xmax() const noexcept {return xmax_;}
    [[nodiscard]] T         ymin() const noexcept {return ymin_;}
    [[nodiscard]] T         ymax() const noexcept {return ymax_;}
    [[nodiscard]] T         dx()   const noexcept {return (xmax_ - xmin_) / static_cast<T>(nx_ - 1);}
    [[nodiscard]] T         dy()   const noexcept {return (ymax_ - ymin_) / static_cast<T>(ny_ - 1);}

    /// x-coordinate at node i (0..nx-1). Throws on out-of-range.
    const T& x(size_type i) const {return x_.at(i);}
    /// y-coordinate at node j (0..ny-1). Throws on out-of-range.
    const T& y(size_type j) const {return y_.at(j);}

    /// Read-only access to full axes
    const std::vector<T>& xvec() const {return x_;}
    const std::vector<T>& yvec() const {return y_;}

    /// Generate 2D meshgrid arrays following MATLAB/NumPy convention
    /// X(i,j) = x_[i] for all j (X varies along rows, constant along columns)
    /// Y(i,j) = y_[j] for all i (Y varies along columns, constant along rows)
    ///
    /// Example for 3x4 grid:
    /// X = [x0 x0 x0 x0]    Y = [y0 y1 y2 y3]
    ///     [x1 x1 x1 x1]        [y0 y1 y2 y3]
    ///     [x2 x2 x2 x2]        [y0 y1 y2 y3]
    Array2D<T> meshgridX() const {
        Array2D<T> X(nx_, ny_, T{});
        for (size_type i = 0; i < nx_; ++i) {
            for (size_type j = 0; j < ny_; ++j) {
                X(i, j) = x_[i];  // x varies with row index i
            }
        }
        return X;
    }

    Array2D<T> meshgridY() const {
        Array2D<T> Y(nx_, ny_, T{});
        for (size_type i = 0; i < nx_; ++i) {
            for (size_type j = 0; j < ny_; ++j) {
                Y(i, j) = y_[j];  // y varies with column index j
            }
        }
        return Y;
    }

    /// Alternative: Generate both meshgrids at once (more efficient)
    std::pair<Array2D<T>, Array2D<T>> meshgrid() const {
        Array2D<T> X(nx_, ny_, T{});
        Array2D<T> Y(nx_, ny_, T{});
        for (size_type i = 0; i < nx_; ++i) {
            for (size_type j = 0; j < ny_; ++j) {
                X(i, j) = x_[i];
                Y(i, j) = y_[j];
            }
        }
        return {X, Y};
    }

    /// Utility: Get (x,y) coordinates at grid point (i,j)
    std::pair<T, T> coords(size_type i, size_type j) const {
        if (i >= nx_ || j >= ny_) {
            throw std::out_of_range("Grid2D::coords: indices out of range");
        }
        return {x_[i], y_[j]};
    }

    /// Check if a point (x,y) is inside the grid domain
    bool contains(T x, T y) const noexcept{
        return (x >= xmin_ && x <= xmax_ && y >= ymin_ && y <= ymax_);
    }

    /// Find nearest grid indices for a point (x,y)
    /// Returns {i, j} where (x_[i], y_[j]) is closest to (x,y)
    std::pair<size_type, size_type> nearest(T x, T y) const {
        // Clamp to domain
        x = std::max(xmin_, std::min(xmax_, x));
        y = std::max(ymin_, std::min(ymax_, y));

        // Find nearest indices
        auto i = static_cast<size_type>(
            std::round((x - xmin_) / dx())
        );
        auto j = static_cast<size_type>(
            std::round((y - ymin_) / dy())
        );

        // Ensure within bounds
        i = std::min(i, nx_ - 1);
        j = std::min(j, ny_ - 1);

        return {i, j};
    }

private:
    // domain
    T xmin_, xmax_, ymin_, ymax_;
    // resolution
    size_type nx_, ny_;
    // axes
    std::vector<T> x_, y_;
};

#endif //GRID2D_H

