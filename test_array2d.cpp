/**
 * @file test_array2d.cpp
 * @brief Comprehensive test suite for Array2D class
 *
 * Tests all functionality including:
 * - Construction and basic operations
 * - Element access and modification
 * - Arithmetic operations
 * - Linear algebra operations
 * - Edge cases and error handling
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <cassert>
#include "src/Array2D.h"

// Helper function to check floating point equality
template<typename T>
bool approx_equal(T a, T b, T tol = 1e-10) {
    return std::abs(a - b) < tol;
}

// Overload for mixed types (automatically converts to double)
template<typename T, typename U>
bool approx_equal(T a, U b, double tol = 1e-10) {
    return std::abs(static_cast<double>(a) - static_cast<double>(b)) < tol;
}

// Test result tracking
int tests_passed = 0;
int tests_failed = 0;

#define TEST_CASE(name) \
    std::cout << "\n[TEST] " << name << "... "; \
    try {

#define TEST_END \
        std::cout << "PASSED\n"; \
        tests_passed++; \
    } catch (const std::exception& e) { \
        std::cout << "FAILED: " << e.what() << "\n"; \
        tests_failed++; \
    }

#define ASSERT(condition) \
    if (!(condition)) { \
        throw std::runtime_error("Assertion failed: " #condition); \
    }

#define ASSERT_APPROX(a, b) \
    if (!approx_equal(a, b)) { \
        throw std::runtime_error("Values not approximately equal"); \
    }

int main() {
    std::cout << "========================================\n";
    std::cout << "    Array2D Class Test Suite\n";
    std::cout << "========================================\n";

    // ============= Construction and Basic Properties =============

    TEST_CASE("Default constructor")
        Array2D<double> A;
        ASSERT(A.rows() == 0);
        ASSERT(A.cols() == 0);
        ASSERT(A.length() == 0);
    TEST_END

    TEST_CASE("Parameterized constructor")
        Array2D<double> B(3, 4, 2.5);
        ASSERT(B.rows() == 3);
        ASSERT(B.cols() == 4);
        ASSERT(B.length() == 12);
        auto shape = B.shape();
        ASSERT(shape[0] == 3);
        ASSERT(shape[1] == 4);
        // Check all elements initialized to 2.5
        for (size_t i = 0; i < B.rows(); ++i) {
            for (size_t j = 0; j < B.cols(); ++j) {
                ASSERT_APPROX(B(i,j), 2.5);
            }
        }
    TEST_END

    // ============= Element Access =============

    TEST_CASE("Element access and modification")
        Array2D<double> C(2, 3, 0.0);
        C(0,0) = 1.0; C(0,1) = 2.0; C(0,2) = 3.0;
        C(1,0) = 4.0; C(1,1) = 5.0; C(1,2) = 6.0;

        ASSERT_APPROX(C(0,0), 1.0);
        ASSERT_APPROX(C(0,1), 2.0);
        ASSERT_APPROX(C(0,2), 3.0);
        ASSERT_APPROX(C(1,0), 4.0);
        ASSERT_APPROX(C(1,1), 5.0);
        ASSERT_APPROX(C(1,2), 6.0);
    TEST_END

    TEST_CASE("Out of bounds access throws")
        Array2D<double> D(2, 2);
        bool caught = false;
        try {
            D(2, 0) = 1.0;  // Out of bounds
        } catch (const std::out_of_range&) {
            caught = true;
        }
        ASSERT(caught);
    TEST_END

    // ============= Modifiers =============

    TEST_CASE("Resize operation")
        Array2D<double> E(2, 2, 1.0);
        E.resize(3, 4, 5.0);
        ASSERT(E.rows() == 3);
        ASSERT(E.cols() == 4);
        ASSERT_APPROX(E(0,0), 5.0);
        ASSERT_APPROX(E(2,3), 5.0);
    TEST_END

    TEST_CASE("Fill operation")
        Array2D<double> F(3, 3, 1.0);
        F.fill(7.5);
        for (size_t i = 0; i < F.rows(); ++i) {
            for (size_t j = 0; j < F.cols(); ++j) {
                ASSERT_APPROX(F(i,j), 7.5);
            }
        }
    TEST_END

    // ============= Arithmetic Operations =============

    TEST_CASE("Matrix addition (+=)")
        Array2D<double> G1(2, 2, 1.0);
        Array2D<double> G2(2, 2, 2.0);
        G1 += G2;
        for (size_t i = 0; i < G1.rows(); ++i) {
            for (size_t j = 0; j < G1.cols(); ++j) {
                ASSERT_APPROX(G1(i,j), 3.0);
            }
        }
    TEST_END

    TEST_CASE("Matrix subtraction (-=)")
        Array2D<double> H1(2, 2, 5.0);
        Array2D<double> H2(2, 2, 2.0);
        H1 -= H2;
        for (size_t i = 0; i < H1.rows(); ++i) {
            for (size_t j = 0; j < H1.cols(); ++j) {
                ASSERT_APPROX(H1(i,j), 3.0);
            }
        }
    TEST_END

    TEST_CASE("Scalar operations (+=, -=, *=, /=)")
        Array2D<double> I(2, 2, 4.0);
        I += 2.0;  // All elements now 6.0
        ASSERT_APPROX(I(0,0), 6.0);

        I -= 1.0;  // All elements now 5.0
        ASSERT_APPROX(I(0,0), 5.0);

        I *= 2.0;  // All elements now 10.0
        ASSERT_APPROX(I(0,0), 10.0);

        I /= 5.0;  // All elements now 2.0
        ASSERT_APPROX(I(0,0), 2.0);
    TEST_END

    TEST_CASE("Free operators (matrix + matrix)")
        Array2D<double> J1(2, 2, 1.0);
        Array2D<double> J2(2, 2, 2.0);
        Array2D<double> J3 = J1 + J2;
        ASSERT_APPROX(J3(0,0), 3.0);
        ASSERT_APPROX(J3(1,1), 3.0);
    TEST_END

    TEST_CASE("Free operators (matrix - matrix)")
        Array2D<double> K1(2, 2, 5.0);
        Array2D<double> K2(2, 2, 2.0);
        Array2D<double> K3 = K1 - K2;
        ASSERT_APPROX(K3(0,0), 3.0);
    TEST_END

    TEST_CASE("Free operators (scalar * matrix)")
        Array2D<double> L(2, 2, 3.0);
        Array2D<double> L2 = 2.0 * L;
        Array2D<double> L3 = L * 2.0;
        ASSERT_APPROX(L2(0,0), 6.0);
        ASSERT_APPROX(L3(0,0), 6.0);
    TEST_END

    // ============= Matrix Multiplication =============

    TEST_CASE("Matrix multiplication")
        Array2D<double> M1(2, 3);
        M1(0,0) = 1; M1(0,1) = 2; M1(0,2) = 3;
        M1(1,0) = 4; M1(1,1) = 5; M1(1,2) = 6;

        Array2D<double> M2(3, 2);
        M2(0,0) = 7;  M2(0,1) = 8;
        M2(1,0) = 9;  M2(1,1) = 10;
        M2(2,0) = 11; M2(2,1) = 12;

        Array2D<double> M3 = M1 * M2;
        ASSERT(M3.rows() == 2);
        ASSERT(M3.cols() == 2);
        // M3 = [1*7+2*9+3*11, 1*8+2*10+3*12] = [58, 64]
        //      [4*7+5*9+6*11, 4*8+5*10+6*12]   [139, 154]
        ASSERT_APPROX(M3(0,0), 58.0);
        ASSERT_APPROX(M3(0,1), 64.0);
        ASSERT_APPROX(M3(1,0), 139.0);
        ASSERT_APPROX(M3(1,1), 154.0);
    TEST_END

    TEST_CASE("Identity matrix multiplication")
        Array2D<double> N(3, 3, 0.0);
        N(0,0) = 1; N(1,1) = 1; N(2,2) = 1; // Identity

        Array2D<double> N2(3, 3);
        N2(0,0) = 1; N2(0,1) = 2; N2(0,2) = 3;
        N2(1,0) = 4; N2(1,1) = 5; N2(1,2) = 6;
        N2(2,0) = 7; N2(2,1) = 8; N2(2,2) = 9;

        Array2D<double> N3 = N * N2;
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                ASSERT_APPROX(N3(i,j), N2(i,j));
            }
        }
    TEST_END

    // ============= Linear Algebra Operations =============

    TEST_CASE("Transpose (non-square)")
        Array2D<double> O(2, 3);
        O(0,0) = 1; O(0,1) = 2; O(0,2) = 3;
        O(1,0) = 4; O(1,1) = 5; O(1,2) = 6;

        Array2D<double> OT = transpose(O);
        ASSERT(OT.rows() == 3);
        ASSERT(OT.cols() == 2);
        ASSERT_APPROX(OT(0,0), 1.0);
        ASSERT_APPROX(OT(0,1), 4.0);
        ASSERT_APPROX(OT(1,0), 2.0);
        ASSERT_APPROX(OT(1,1), 5.0);
        ASSERT_APPROX(OT(2,0), 3.0);
        ASSERT_APPROX(OT(2,1), 6.0);
    TEST_END

    TEST_CASE("In-place transpose (square matrix)")
        Array2D<double> P(3, 3);
        P(0,0) = 1; P(0,1) = 2; P(0,2) = 3;
        P(1,0) = 4; P(1,1) = 5; P(1,2) = 6;
        P(2,0) = 7; P(2,1) = 8; P(2,2) = 9;

        P.transpose_inplace();
        ASSERT_APPROX(P(0,1), 4.0);
        ASSERT_APPROX(P(0,2), 7.0);
        ASSERT_APPROX(P(1,0), 2.0);
        ASSERT_APPROX(P(1,2), 8.0);
        ASSERT_APPROX(P(2,0), 3.0);
        ASSERT_APPROX(P(2,1), 6.0);
    TEST_END

    TEST_CASE("Matrix inverse (2x2)")
        Array2D<double> Q(2, 2);
        Q(0,0) = 4; Q(0,1) = 7;
        Q(1,0) = 2; Q(1,1) = 6;
        // det(Q) = 4*6 - 7*2 = 24 - 14 = 10
        // Q^(-1) = (1/10) * [6, -7; -2, 4]

        Array2D<double> Q_inv = Q.inverse();
        ASSERT_APPROX(Q_inv(0,0), 0.6);
        ASSERT_APPROX(Q_inv(0,1), -0.7);
        ASSERT_APPROX(Q_inv(1,0), -0.2);
        ASSERT_APPROX(Q_inv(1,1), 0.4);

        // Verify Q * Q_inv = I
        Array2D<double> I = Q * Q_inv;
        ASSERT_APPROX(I(0,0), 1.0);
        ASSERT_APPROX(I(0,1), 0.0);
        ASSERT_APPROX(I(1,0), 0.0);
        ASSERT_APPROX(I(1,1), 1.0);
    TEST_END

    TEST_CASE("Matrix inverse (3x3)")
        Array2D<double> R(3, 3);
        R(0,0) = 2; R(0,1) = 1; R(0,2) = 1;
        R(1,0) = 1; R(1,1) = 3; R(1,2) = 2;
        R(2,0) = 1; R(2,1) = 0; R(2,2) = 0;

        Array2D<double> R_inv = R.inverse();

        // Verify R * R_inv = I
        Array2D<double> I = R * R_inv;
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                double expected = (i == j) ? 1.0 : 0.0;
                ASSERT_APPROX(I(i,j), expected);
            }
        }
    TEST_END

    TEST_CASE("Solve linear system Ax=b")
        // System: 2x + y = 5
        //         x + 3y = 7
        // Solution: x = 8/5, y = 9/5
        Array2D<double> A(2, 2);
        A(0,0) = 2; A(0,1) = 1;
        A(1,0) = 1; A(1,1) = 3;

        Array2D<double> b(2, 1);
        b(0,0) = 5;
        b(1,0) = 7;

        Array2D<double> x = solve(A, b);
        ASSERT_APPROX(x(0,0), 8.0/5.0);
        ASSERT_APPROX(x(1,0), 9.0/5.0);

        // Verify Ax = b
        Array2D<double> check = A * x;
        ASSERT_APPROX(check(0,0), 5.0);
        ASSERT_APPROX(check(1,0), 7.0);
    TEST_END

    TEST_CASE("Solve multiple RHS")
        Array2D<double> A(2, 2);
        A(0,0) = 3; A(0,1) = 1;
        A(1,0) = 1; A(1,1) = 2;

        Array2D<double> B(2, 2);
        B(0,0) = 9; B(0,1) = 8;
        B(1,0) = 8; B(1,1) = 7;

        Array2D<double> X = solve(A, B);

        // Verify AX = B
        Array2D<double> check = A * X;
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                ASSERT_APPROX(check(i,j), B(i,j));
            }
        }
    TEST_END

    // ============= Complex Numbers Support =============

    TEST_CASE("Complex matrix operations")
        using Complex = std::complex<double>;
        Array2D<Complex> C1(2, 2);
        C1(0,0) = Complex(1, 2);
        C1(0,1) = Complex(3, 4);
        C1(1,0) = Complex(5, 6);
        C1(1,1) = Complex(7, 8);

        Array2D<Complex> C2 = conj_transpose(C1);
        // Conjugate transpose should give:
        // (1,-2)  (5,-6)
        // (3,-4)  (7,-8)
        ASSERT(C2(0,0) == Complex(1, -2));
        ASSERT(C2(0,1) == Complex(5, -6));
        ASSERT(C2(1,0) == Complex(3, -4));
        ASSERT(C2(1,1) == Complex(7, -8));
    TEST_END

    // ============= Edge Cases =============

    TEST_CASE("Empty matrix operations")
        Array2D<double> E1;
        Array2D<double> E2;
        // These should not crash
        E1.fill(5.0);
        ASSERT(E1.rows() == 0);
        ASSERT(E1.cols() == 0);
    TEST_END

    TEST_CASE("Single element matrix")
        Array2D<double> S(1, 1, 42.0);
        ASSERT_APPROX(S(0,0), 42.0);
        S *= 2.0;
        ASSERT_APPROX(S(0,0), 84.0);

        Array2D<double> S_inv = S.inverse();
        ASSERT_APPROX(S_inv(0,0), 1.0/84.0);
    TEST_END

    TEST_CASE("Shape mismatch detection")
        Array2D<double> M1(2, 3);
        Array2D<double> M2(2, 2);
        bool caught = false;
        try {
            M1 += M2;  // Should throw
        } catch (const std::invalid_argument&) {
            caught = true;
        }
        ASSERT(caught);
    TEST_END

    // ============= Performance / Stress Test =============

    TEST_CASE("Large matrix operations")
        const size_t n = 100;
        Array2D<double> Big1(n, n, 1.0);
        Array2D<double> Big2(n, n, 2.0);

        // Addition
        Array2D<double> Big3 = Big1 + Big2;
        ASSERT_APPROX(Big3(0,0), 3.0);
        ASSERT_APPROX(Big3(n-1,n-1), 3.0);

        // Matrix multiplication (will be slow for large n)
        Array2D<double> I(n, n, 0.0);
        for (size_t i = 0; i < n; ++i) I(i,i) = 1.0;

        Array2D<double> Big4 = I * Big1;
        ASSERT_APPROX(Big4(0,0), 1.0);
        ASSERT_APPROX(Big4(n-1,n-1), 1.0);
    TEST_END

    // ============= Output Test =============

    TEST_CASE("Stream output operator")
        Array2D<double> Out(2, 3);
        Out(0,0) = 1; Out(0,1) = 2; Out(0,2) = 3;
        Out(1,0) = 4; Out(1,1) = 5; Out(1,2) = 6;

        std::cout << "\nExample matrix output:\n" << Out << "\n";
        std::cout << "(Visual inspection only) ";
    TEST_END

    // ============= Summary =============

    std::cout << "\n========================================\n";
    std::cout << "Test Results Summary:\n";
    std::cout << "  Passed: " << tests_passed << "\n";
    std::cout << "  Failed: " << tests_failed << "\n";

    if (tests_failed == 0) {
        std::cout << "\n✓ All tests PASSED! Array2D class is working correctly.\n";
    } else {
        std::cout << "\n✗ Some tests FAILED. Please review the implementation.\n";
    }
    std::cout << "========================================\n";

    return (tests_failed == 0) ? 0 : 1;
}