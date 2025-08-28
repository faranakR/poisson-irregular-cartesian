#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <vector>
#include <array>
#include <stdexcept>
#include <cstddef>     // size_t
#include <utility>     // std::swap
#include <type_traits> // std::is_floating_point, std::is_same_v
#include <complex>     // std::complex, std::conj
#include <limits>      // epsilon
#include <ostream>     // std::ostream
#include <iomanip>     // std::setw, std::setprecision
#include <cmath>       // std::abs

/**
 * @brief A simple, numerics-friendly 2D array with flat row-major storage.
 * @tparam T element type (e.g., double, float, int, std::complex<double>)
 *
 *  - Contiguous row-major buffer: (i,j) ↦ i*cols + j
 *  - 0-based indexing via operator()(i,j)
 *  - Compound ops (+=, -=, *=, /=) are members; plain ops are free functions
 *  - Includes transpose, inverse (square floating-point), and solve(A,B) ≈ MATLAB A\B
 */
template<typename T>
class Array2D {
public:
    // -------------------- CTORS & SHAPE --------------------
    Array2D() noexcept : nrows_(0), ncols_(0) {}
    Array2D(std::size_t r, std::size_t c, const T& v = T{})
        : nrows_(r), ncols_(c), buf_(r*c, v) {}

    [[nodiscard]] std::array<std::size_t,2> shape()  const noexcept { return {nrows_, ncols_}; }
    [[nodiscard]] std::size_t               rows()   const noexcept { return nrows_; }
    [[nodiscard]] std::size_t               cols()   const noexcept { return ncols_; }
    [[nodiscard]] std::size_t               length() const noexcept { return buf_.size(); }
    [[nodiscard]] const T*                  data()   const noexcept { return buf_.data(); }
    [[nodiscard]] T*                        data()         noexcept { return buf_.data(); }

    // -------------------- ELEMENT ACCESS -------------------
    const T& operator()(std::size_t i, std::size_t j) const {
        if (i >= nrows_ || j >= ncols_) throw std::out_of_range("Array2D: index out of range");
        return buf_[i*ncols_ + j];
    }
    T& operator()(std::size_t i, std::size_t j) {
        if (i >= nrows_ || j >= ncols_) throw std::out_of_range("Array2D: index out of range");
        return buf_[i*ncols_ + j];
    }

    // -------------------- MODIFIERS ------------------------
    void resize(std::size_t r, std::size_t c, const T& v = T{}) {
        nrows_ = r; ncols_ = c; buf_.assign(r*c, v);
    }
    void fill(const T& v) { std::fill(buf_.begin(), buf_.end(), v); }

    // -------------------- COMPOUND OPS ---------------------
    Array2D& operator+=(const Array2D& other) {
        if (shape() != other.shape()) throw std::invalid_argument("Array2D: shape mismatch in +=");
        for (std::size_t i = 0; i < buf_.size(); ++i) buf_[i] += other.buf_[i];
        return *this;
    }
    Array2D& operator-=(const Array2D& other) {
        if (shape() != other.shape()) throw std::invalid_argument("Array2D: shape mismatch in -=");
        for (std::size_t i = 0; i < buf_.size(); ++i) buf_[i] -= other.buf_[i];
        return *this;
    }
    Array2D& operator+=(const T s) { for (auto& x : buf_) x += s; return *this; }
    Array2D& operator-=(const T s) { for (auto& x : buf_) x -= s; return *this; }
    Array2D& operator*=(const T s) { for (auto& x : buf_) x *= s; return *this; }
    Array2D& operator/=(const T s) { for (auto& x : buf_) x /= s; return *this; }

    // -------------------- LINEAR ALGEBRA -------------------
    /**
     * @brief In-place transpose for square matrices (O(1) extra memory).
     * @throws std::logic_error if matrix is not square.
     */
    void transpose_inplace() {
        if (nrows_ != ncols_) throw std::logic_error("transpose_inplace: requires square matrix");
        for (std::size_t i = 0; i < nrows_; ++i) {
            for (std::size_t j = i+1; j < ncols_; ++j) {
                std::swap((*this)(i,j), (*this)(j,i));
            }
        }
        // nrows_ == ncols_, so no swap needed
    }

    /**
     * @brief Matrix inverse (Gauss–Jordan with partial pivoting).
     * @return A^{-1}
     * @throws std::invalid_argument if not square or T is not floating-point.
     * @throws std::runtime_error if singular/ill-conditioned pivot encountered.
     *
     * @note Available only when T is floating-point (float/double/long double).
     */
    Array2D inverse() const {
        static_assert(std::is_floating_point<T>::value,
                      "inverse() requires floating-point T");
        if (nrows_ != ncols_) throw std::invalid_argument("inverse: matrix must be square");

        const std::size_t n = nrows_;
        Array2D M = *this;          // working copy
        Array2D Inv(n, n, T{});     // identity
        for (std::size_t i = 0; i < n; ++i) Inv(i,i) = T{1};

        const T eps = std::numeric_limits<T>::epsilon();

        for (std::size_t col = 0; col < n; ++col) {
            // partial pivot on column 'col'
            std::size_t pivot = col;
            T maxAbs = std::abs(M(col, col));
            for (std::size_t r = col+1; r < n; ++r) {
                T v = std::abs(M(r, col));
                if (v > maxAbs) { maxAbs = v; pivot = r; }
            }
            if (!(maxAbs > eps)) {
                throw std::runtime_error("inverse: singular or ill-conditioned matrix");
            }
            if (pivot != col) {
                for (std::size_t j = 0; j < n; ++j) {
                    std::swap(M(col,j),  M(pivot,j));
                    std::swap(Inv(col,j), Inv(pivot,j));
                }
            }
            // normalize pivot row
            const T d = M(col, col);
            for (std::size_t j = 0; j < n; ++j) {
                M(col, j)   /= d;
                Inv(col, j) /= d;
            }
            // eliminate other rows
            for (std::size_t r = 0; r < n; ++r) {
                if (r == col) continue;
                const T f = M(r, col);
                if (f == T{}) continue;
                for (std::size_t j = 0; j < n; ++j) {
                    M(r, j)   -= f * M(col, j);
                    Inv(r, j) -= f * Inv(col, j);
                }
            }
        }
        return Inv;
    }

private:
    std::size_t  nrows_;
    std::size_t  ncols_;
    std::vector<T> buf_;

    // grant free operators access to buf_ efficiently via public API when needed
    template<typename U> friend Array2D<U> operator+(Array2D<U>, const Array2D<U>&);
    template<typename U> friend Array2D<U> operator-(Array2D<U>, const Array2D<U>&);
    template<typename U> friend Array2D<U> operator*(const Array2D<U>&, const Array2D<U>&);
};

// ========================= Free utilities =========================

/**
 * @brief Transpose copy: Aᵀ (shape cols×rows).
 */
template<typename T>
Array2D<T> transpose(const Array2D<T>& A) {
    Array2D<T> AT(A.cols(), A.rows(), T{});
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j)
            AT(j, i) = A(i, j);
    return AT;
}

/**
 * @brief Conjugate transpose (Hermitian): Aᴴ = conj(A)ᵀ.
 * Falls back to plain transpose if T is real.
 */
template<typename T>
Array2D<T> conj_transpose(const Array2D<T>& A) {
    Array2D<T> AH(A.cols(), A.rows(), T{});
    if constexpr (std::is_same_v<T, std::complex<float>> ||
                  std::is_same_v<T, std::complex<double>> ||
                  std::is_same_v<T, std::complex<long double>>) {
        for (std::size_t i = 0; i < A.rows(); ++i)
            for (std::size_t j = 0; j < A.cols(); ++j)
                AH(j, i) = std::conj(A(i, j));
    } else {
        // real T: same as transpose
        for (std::size_t i = 0; i < A.rows(); ++i)
            for (std::size_t j = 0; j < A.cols(); ++j)
                AH(j, i) = A(i, j);
    }
    return AH;
}

/**
 * @brief MATLAB-style backslash: solve A*X = B for X.
 * @return X with shape (A.rows() × B.cols()).
 *
 * @details
 *  - Performs Gaussian elimination w/ partial pivoting on an augmented system.
 *  - For B being a vector (n×1), pass it as Array2D<T>(n,1).
 *  - Throws if shapes don’t match or if singular.
 *  - For serious workloads, prefer a proper LU/QR from a BLAS/LAPACK/Eigen backend.
 */
template<typename T>
Array2D<T> solve(const Array2D<T>& A, const Array2D<T>& B) {
    static_assert(std::is_floating_point<T>::value,
                  "solve() requires floating-point T");
    if (A.rows() != A.cols())
        throw std::invalid_argument("solve: A must be square");
    if (B.rows() != A.rows())
        throw std::invalid_argument("solve: B.rows must equal A.rows");

    const std::size_t n = A.rows();
    const std::size_t m = B.cols();

    // Build augmented [A | B]
    Array2D<T> M(n, n + m, T{});
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) M(i, j) = A(i, j);
        for (std::size_t j = 0; j < m; ++j) M(i, n + j) = B(i, j);
    }

    const T eps = std::numeric_limits<T>::epsilon();

    // Forward elimination with partial pivoting
    for (std::size_t col = 0; col < n; ++col) {
        // Pivot
        std::size_t pivot = col;
        T maxAbs = std::abs(M(col, col));
        for (std::size_t r = col + 1; r < n; ++r) {
            T v = std::abs(M(r, col));
            if (v > maxAbs) { maxAbs = v; pivot = r; }
        }
        if (!(maxAbs > eps)) throw std::runtime_error("solve: singular or ill-conditioned");

        if (pivot != col) {
            for (std::size_t j = col; j < n + m; ++j) std::swap(M(col, j), M(pivot, j));
        }
        // Normalize pivot row
        const T d = M(col, col);
        for (std::size_t j = col; j < n + m; ++j) M(col, j) /= d;

        // Eliminate below
        for (std::size_t r = col + 1; r < n; ++r) {
            const T f = M(r, col);
            if (f == T{}) continue;
            for (std::size_t j = col; j < n + m; ++j) M(r, j) -= f * M(col, j);
        }
    }

    // Back substitution
    for (std::size_t col = n; col-- > 0; ) {
        for (std::size_t r = 0; r < col; ++r) {
            const T f = M(r, col);
            if (f == T{}) continue;
            for (std::size_t j = col; j < n + m; ++j) M(r, j) -= f * M(col, j);
        }
    }

    // Extract X from [I | X]
    Array2D<T> X(n, m, T{});
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < m; ++j)
            X(i, j) = M(i, n + j);

    return X;
}

// ========================= Free operators =========================

// element-wise matrix ± matrix
template<typename T>
Array2D<T> operator+(Array2D<T> lhs, const Array2D<T>& rhs) { lhs += rhs; return lhs; }
template<typename T>
Array2D<T> operator-(Array2D<T> lhs, const Array2D<T>& rhs) { lhs -= rhs; return lhs; }

// matrix ± scalar
template<typename T>
Array2D<T> operator+(Array2D<T> lhs, const T s) { lhs += s; return lhs; }
template<typename T>
Array2D<T> operator-(Array2D<T> lhs, const T s) { lhs -= s; return lhs; }
template<typename T>
Array2D<T> operator+(const T s, Array2D<T> rhs) { rhs += s; return rhs; }
template<typename T>
Array2D<T> operator-(const T s, const Array2D<T>& rhs) {
    Array2D<T> out(rhs.rows(), rhs.cols(), s);
    out -= rhs; return out;
}

// scalar * matrix (scaling)
template<typename T>
Array2D<T> operator*(Array2D<T> lhs, const T s) { lhs *= s; return lhs; }
template<typename T>
Array2D<T> operator*(const T s, Array2D<T> rhs) { rhs *= s; return rhs; }
template<typename T>
Array2D<T> operator/(Array2D<T> lhs, const T s) { lhs /= s; return lhs; }

/**
 * @brief Matrix–matrix multiplication: C = A * B (dot-product triple loop).
 * @throws std::invalid_argument if A.cols() != B.rows()
 */
template<typename T>
Array2D<T> operator*(const Array2D<T>& A, const Array2D<T>& B) {
    if (A.cols() != B.rows()) throw std::invalid_argument("Array2D: shape mismatch in matmul");
    Array2D<T> C(A.rows(), B.cols(), T{});
    for (std::size_t i = 0; i < A.rows(); ++i) {
        for (std::size_t j = 0; j < B.cols(); ++j) {
            T sum = T{};
            for (std::size_t k = 0; k < A.cols(); ++k) sum += A(i,k) * B(k,j);
            C(i,j) = sum;
        }
    }
    return C;
}

// ========================= Pretty printing ========================

/**
 * @brief Stream insertion (MATLAB-ish): prints rows on separate lines.
 *
 * Example:
 *   [ 1 2 3
 *     4 5 6 ]
 */
template<typename T>
std::ostream& operator<<(std::ostream& os, const Array2D<T>& A) {
    os << "[";
    for (std::size_t i = 0; i < A.rows(); ++i) {
        if (i > 0) os << " ";
        for (std::size_t j = 0; j < A.cols(); ++j) {
            os << A(i,j);
            if (j + 1 < A.cols()) os << " ";
        }
        if (i + 1 < A.rows()) os << "\n";
    }
    os << "]";
    return os;
}

#endif // ARRAY2D_H
