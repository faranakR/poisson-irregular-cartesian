/**
 * @file test_finitedifference2d.cpp
 * @brief Test suite for FiniteDifference2D operators
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include "../src/FiniteDifference2D.h"
#include "../src/Grid2D.h"
#include "../src/Array2D.h"

template<typename T>
bool approx_equal(T a, T b, T tol = 1e-10) {
    return std::abs(a - b) <= tol;
}

int tests_passed = 0;
int tests_failed = 0;

#define TEST_CASE(name) \
    std::cout << "\n[TEST] " << name << "... "; \
    try {

#define END_TEST \
        std::cout << "PASSED\n"; \
        tests_passed++; \
    } catch (const std::exception& e) { \
        std::cout << "FAILED: " << e.what() << "\n"; \
        tests_failed++; \
    }

#define ASSERT(condition) \
    if (!(condition)) throw std::runtime_error("Assertion failed: " #condition);

#define ASSERT_APPROX(a, b) \
    if (!approx_equal((a), (b))) throw std::runtime_error("Values not approximately equal");

int main() {
    std::cout << "========================================\n";
    std::cout << "  FiniteDifference2D Test Suite\n";
    std::cout << "========================================\n";

    // 5×5 grid on [0,1]×[0,1]
    Grid2D<double> grid(0.0, 1.0, 0.0, 1.0, 5, 5);

    // ============= Test Basic Differences =============
    TEST_CASE("Forward difference X")
        Array2D<double> u(5, 5, 0.0);
        // u = x
        for (size_t i = 0; i < 5; ++i)
            for (size_t j = 0; j < 5; ++j)
                u(i,j) = grid.x(i);

        auto du_dx = FiniteDifference2D<double>::forward_diff_x(u, grid);

        // du/dx = 1 in interior
        for (size_t i = 0; i < 4; ++i)
            for (size_t j = 0; j < 5; ++j)
                ASSERT_APPROX(du_dx(i,j), 1.0);

        // right boundary (no BC handling in math layer)
        for (size_t j = 0; j < 5; ++j)
            ASSERT_APPROX(du_dx(4,j), 0.0);
    END_TEST

    TEST_CASE("Central difference Y")
        Array2D<double> v(5, 5, 0.0);
        // v = y^2  => dv/dy = 2y
        for (size_t i = 0; i < 5; ++i)
            for (size_t j = 0; j < 5; ++j) {
                double y = grid.y(j);
                v(i,j) = y * y;
            }

        auto dv_dy = FiniteDifference2D<double>::central_diff_y(v, grid);

        // interior in y
        for (size_t i = 0; i < 5; ++i)
            for (size_t j = 1; j < 4; ++j) {
                double y = grid.y(j);
                ASSERT_APPROX(dv_dy(i,j), 2.0 * y);
            }

        // boundaries set to 0 by math layer
        for (size_t i = 0; i < 5; ++i) {
            ASSERT_APPROX(dv_dy(i,0), 0.0);
            ASSERT_APPROX(dv_dy(i,4), 0.0);
        }
    END_TEST

    // ============= Test Gradient =============
    TEST_CASE("Gradient of linear function")
        Array2D<double> w(5, 5, 0.0);
        // w = 2x + 3y
        for (size_t i = 0; i < 5; ++i)
            for (size_t j = 0; j < 5; ++j)
                w(i,j) = 2.0 * grid.x(i) + 3.0 * grid.y(j);

        auto [wx, wy] = FiniteDifference2D<double>::gradient(w, grid);

        for (size_t i = 1; i < 4; ++i)
            for (size_t j = 1; j < 4; ++j) {
                ASSERT_APPROX(wx(i,j), 2.0);
                ASSERT_APPROX(wy(i,j), 3.0);
            }
    END_TEST

    // ============= Test Laplacian =============
    TEST_CASE("Laplacian of harmonic function")
        Array2D<double> h(5, 5, 0.0);
        // h = x^2 - y^2 => ∇²h = 0
        for (size_t i = 0; i < 5; ++i)
            for (size_t j = 0; j < 5; ++j) {
                double x = grid.x(i);
                double y = grid.y(j);
                h(i,j) = x*x - y*y;
            }

        auto lap_h = FiniteDifference2D<double>::laplacian5_interior(h, grid);

        for (size_t i = 1; i < 4; ++i)
            for (size_t j = 1; j < 4; ++j)
                ASSERT_APPROX(lap_h(i,j), 0.0);
    END_TEST

    TEST_CASE("Laplacian of quadratic function")
        Array2D<double> q(5, 5, 0.0);
        // q = x^2 + y^2 => ∇²q = 4
        for (size_t i = 0; i < 5; ++i)
            for (size_t j = 0; j < 5; ++j) {
                double x = grid.x(i);
                double y = grid.y(j);
                q(i,j) = x*x + y*y;
            }

        auto lap_q = FiniteDifference2D<double>::laplacian5_interior(q, grid);

        for (size_t i = 1; i < 4; ++i)
            for (size_t j = 1; j < 4; ++j)
                ASSERT_APPROX(lap_q(i,j), 4.0);
    END_TEST

    // ============= Test Upwind Schemes =============
    TEST_CASE("Upwind advection with constant velocity")
        Array2D<double> uadv(5, 5, 0.0);
        Array2D<double> cx(5, 5, 1.0);  // +x
        Array2D<double> cy(5, 5, 0.0);  // 0y

        // u = x  => c·∇u = 1
        for (size_t i = 0; i < 5; ++i)
            for (size_t j = 0; j < 5; ++j)
                uadv(i,j) = grid.x(i);

        auto adv = FiniteDifference2D<double>::upwind_advection(uadv, cx, cy, grid);

        for (size_t i = 1; i < 4; ++i)
            for (size_t j = 1; j < 4; ++j)
                ASSERT_APPROX(adv(i,j), 1.0);
    END_TEST

    // ============= Test Level Set Operators =============
    TEST_CASE("Sign function")
        Array2D<double> phi(5, 5, 0.0);
        for (size_t i = 0; i < 5; ++i)
            for (size_t j = 0; j < 5; ++j)
                phi(i,j) = grid.x(i) - 0.5;  // negative for x<0.5

        auto S = FiniteDifference2D<double>::sign_function(phi, grid);

        for (size_t i = 0; i < 5; ++i)
            for (size_t j = 0; j < 5; ++j) {
                if (phi(i,j) > 0) {ASSERT(S(i,j) > 0);}
                else if (phi(i,j) < 0) {ASSERT(S(i,j) < 0);}
            }
    END_TEST

    TEST_CASE("Mean curvature of circle")
        // 21×21 grid on [-1,1]×[-1,1]
        Grid2D<double> fine_grid(-1.0, 1.0, -1.0, 1.0, 21, 21);
        Array2D<double> phic(21, 21, 0.0);

        double R = 0.5;
        for (size_t i = 0; i < 21; ++i)
            for (size_t j = 0; j < 21; ++j) {
                double x = fine_grid.x(i);
                double y = fine_grid.y(j);
                phic(i,j) = std::sqrt(x*x + y*y) - R;
            }

        auto kappa = FiniteDifference2D<double>::mean_curvature(phic, fine_grid);

        // κ ≈ 1/R near interface
        for (size_t i = 5; i < 16; ++i)
            for (size_t j = 5; j < 16; ++j)
                if (std::abs(phic(i,j)) < 0.1)
                    ASSERT(std::abs(kappa(i,j) - 1.0/R) < 0.5);
    END_TEST

    // ============= Test Stencil Structure =============
    TEST_CASE("Laplacian stencil for interior point")
        auto st = FiniteDifference2D<double>::get_laplacian_stencil(2, 2, grid);

        double dx = grid.dx(), dy = grid.dy();
        double inv_dx2 = 1.0 / (dx * dx);
        double inv_dy2 = 1.0 / (dy * dy);

        ASSERT_APPROX(st.center, -2.0 * (inv_dx2 + inv_dy2));
        ASSERT_APPROX(st.east,   inv_dx2);
        ASSERT_APPROX(st.west,   inv_dx2);
        ASSERT_APPROX(st.north,  inv_dy2);
        ASSERT_APPROX(st.south,  inv_dy2);
    END_TEST

    TEST_CASE("Laplacian stencil for boundary point")
        auto st_b = FiniteDifference2D<double>::get_laplacian_stencil(0, 2, grid);

        ASSERT_APPROX(st_b.center, 0.0);
        ASSERT_APPROX(st_b.east,   0.0);
        ASSERT_APPROX(st_b.west,   0.0);
        ASSERT_APPROX(st_b.north,  0.0);
        ASSERT_APPROX(st_b.south,  0.0);
    END_TEST

    TEST_CASE("Boundary point detection")
        ASSERT(FiniteDifference2D<double>::is_boundary_point(0, 2, grid));
        ASSERT(FiniteDifference2D<double>::is_boundary_point(4, 2, grid));
        ASSERT(FiniteDifference2D<double>::is_boundary_point(2, 0, grid));
        ASSERT(FiniteDifference2D<double>::is_boundary_point(2, 4, grid));
        ASSERT(!FiniteDifference2D<double>::is_boundary_point(2, 2, grid));
    END_TEST

    // ============= Summary =============
    std::cout << "\n========================================\n";
    std::cout << "Test Results Summary:\n";
    std::cout << "  Passed: " << tests_passed << "\n";
    std::cout << "  Failed: " << tests_failed << "\n";
    std::cout << ((tests_failed == 0) ? "\n✓ All tests PASSED!\n" : "\n✗ Some tests FAILED.\n");
    std::cout << "========================================\n";

    return (tests_failed == 0) ? 0 : 1;
}
