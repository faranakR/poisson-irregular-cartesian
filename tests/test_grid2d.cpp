/**
 * @file test_grid2d.cpp
 * @brief Comprehensive test suite for Grid2D class
 * @author Faranak Rajabi
 * @date 8/28/25
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

#include "../src/Grid2D.h"

// Helper functions to check equality of floating-point types
template<typename T>
bool approx_equal(T a, T b, T tol = 1e-10) {
    return std::abs(a - b) <= tol;
}

// Overload for mixed types, convert to double
template<typename T, typename U>
bool approx_equal(T a, U b, double tol = 1e-10) {
    return std::abs(static_cast<double>(a) - static_cast<double>(b)) <= tol;
}

// Tracking tests success
int tests_passed {0};
int tests_failed {0};

#define TEST_CASE(name)\
    std::cout << "\n[TEST]: " << name << "... "; \
    try {

#define END_TEST \
        std::cout << "PASSED\n"; \
        tests_passed++; } \
    catch (const std::exception& e){ \
        std::cout << "FAILED: " << e.what() << "\n"; \
        tests_failed++;\
    }

#define ASSERT(condition) \
    if (!(condition)) throw std::runtime_error("Assertion failed: "#condition);

#define ASSERT_APPROX(a, b) \
    if (!(approx_equal(a, b))) throw std::runtime_error("Values not approximately equal");

int main() {
    std::cout << "========================================\n";
    std::cout << "       Grid2D Class Test Suite\n";
    std::cout << "========================================\n";

    // ============= Construction Tests =============

    TEST_CASE("Default construction")
        Grid2D<double> g1;
        ASSERT(g1.nx() == 100);
        ASSERT(g1.ny() == 100);
        ASSERT_APPROX(g1.xmin(), -1.0);
        ASSERT_APPROX(g1.xmax(), 1.0);
        ASSERT_APPROX(g1.ymin(), -1.0);
        ASSERT_APPROX(g1.ymax(), 1.0);
    END_TEST

    TEST_CASE("Custom construction")
        Grid2D<double> g2(0.0, 2.0, -1.0, 3.0, 5, 7);
        ASSERT(g2.nx() == 5);
        ASSERT(g2.ny() == 7);
        ASSERT_APPROX(g2.xmin(), 0.0);
        ASSERT_APPROX(g2.xmax(), 2.0);
        ASSERT_APPROX(g2.ymin(), -1.0);
        ASSERT_APPROX(g2.ymax(), 3.0);
    END_TEST

    TEST_CASE("Grid spacing calculation")
        Grid2D<double> g3(0.0, 4.0, 0.0, 6.0, 5, 4);
        // dx = (4.0 - 0.0) / (5 - 1) = 1.0
        // dy = (6.0 - 0.0) / (4 - 1) = 2.0
        ASSERT_APPROX(g3.dx(), 1.0);
        ASSERT_APPROX(g3.dy(), 2.0);
    END_TEST

    TEST_CASE("Invalid construction - xmax <= xmin")
        bool caught = false;
        try {
            Grid2D<double> g(2.0, 1.0, 0.0, 1.0, 10, 10);
        } catch (const std::invalid_argument&) {
            caught = true;
        }
        ASSERT(caught);
    END_TEST

    TEST_CASE("Invalid construction - ymax <= ymin")
        bool caught = false;
        try {
            Grid2D<double> g(0.0, 1.0, 2.0, 1.0, 10, 10);
        } catch (const std::invalid_argument&) {
            caught = true;
        }
        ASSERT(caught);
    END_TEST

    TEST_CASE("Invalid construction - nx < 2")
        bool caught = false;
        try {
            Grid2D<double> g(0.0, 1.0, 0.0, 1.0, 1, 10);
        } catch (const std::invalid_argument&) {
            caught = true;
        }
        ASSERT(caught);
    END_TEST

    TEST_CASE("Invalid construction - ny < 2")
        bool caught = false;
        try {
            Grid2D<double> g(0.0, 1.0, 0.0, 1.0, 10, 1);
        } catch (const std::invalid_argument&) {
            caught = true;
        }
        ASSERT(caught);
    END_TEST

    // ============= Coordinate Access Tests =============

    TEST_CASE("X-axis coordinates")
        Grid2D<double> g4(0.0, 3.0, 0.0, 1.0, 4, 2);
        // x[0] = 0, x[1] = 1, x[2] = 2, x[3] = 3
        ASSERT_APPROX(g4.x(0), 0.0);
        ASSERT_APPROX(g4.x(1), 1.0);
        ASSERT_APPROX(g4.x(2), 2.0);
        ASSERT_APPROX(g4.x(3), 3.0);
    END_TEST

    TEST_CASE("Y-axis coordinates")
        Grid2D<double> g5(0.0, 1.0, 0.0, 2.0, 2, 3);
        // y[0] = 0, y[1] = 1, y[2] = 2
        ASSERT_APPROX(g5.y(0), 0.0);
        ASSERT_APPROX(g5.y(1), 1.0);
        ASSERT_APPROX(g5.y(2), 2.0);
    END_TEST

    TEST_CASE("Out of bounds x access")
        Grid2D<double> g6(0.0, 1.0, 0.0, 1.0, 3, 3);
        bool caught = false;
        try {
            auto val = g6.x(3); // Out of bounds
        } catch (const std::out_of_range&) {
            caught = true;
        }
        ASSERT(caught);
    END_TEST

    TEST_CASE("Out of bounds y access")
        Grid2D<double> g7(0.0, 1.0, 0.0, 1.0, 3, 3);
        bool caught = false;
        try {
            auto val = g7.y(3); // Out of bounds
        } catch (const std::out_of_range&) {
            caught = true;
        }
        ASSERT(caught);
    END_TEST

    TEST_CASE("Vector access")
        Grid2D<double> g8(0.0, 2.0, 0.0, 3.0, 3, 4);
        auto xvec = g8.xvec();
        auto yvec = g8.yvec();
        ASSERT(xvec.size() == 3);
        ASSERT(yvec.size() == 4);
        ASSERT_APPROX(xvec[0], 0.0);
        ASSERT_APPROX(xvec[1], 1.0);
        ASSERT_APPROX(xvec[2], 2.0);
        ASSERT_APPROX(yvec[0], 0.0);
        ASSERT_APPROX(yvec[1], 1.0);
        ASSERT_APPROX(yvec[2], 2.0);
        ASSERT_APPROX(yvec[3], 3.0);
    END_TEST

    // ============= Meshgrid Tests =============

    TEST_CASE("meshgridX generation")
        Grid2D<double> g9(0.0, 2.0, 0.0, 1.0, 3, 2);
        auto X = g9.meshgridX();
        ASSERT(X.rows() == 3);
        ASSERT(X.cols() == 2);
        // X should be [0 0; 1 1; 2 2]
        ASSERT_APPROX(X(0,0), 0.0);
        ASSERT_APPROX(X(0,1), 0.0);
        ASSERT_APPROX(X(1,0), 1.0);
        ASSERT_APPROX(X(1,1), 1.0);
        ASSERT_APPROX(X(2,0), 2.0);
        ASSERT_APPROX(X(2,1), 2.0);
    END_TEST

    TEST_CASE("meshgridY generation")
        Grid2D<double> g10(0.0, 2.0, 0.0, 1.0, 3, 2);
        auto Y = g10.meshgridY();
        ASSERT(Y.rows() == 3);
        ASSERT(Y.cols() == 2);
        // Y should be [0 1; 0 1; 0 1]
        ASSERT_APPROX(Y(0,0), 0.0);
        ASSERT_APPROX(Y(0,1), 1.0);
        ASSERT_APPROX(Y(1,0), 0.0);
        ASSERT_APPROX(Y(1,1), 1.0);
        ASSERT_APPROX(Y(2,0), 0.0);
        ASSERT_APPROX(Y(2,1), 1.0);
    END_TEST

    TEST_CASE("Combined meshgrid generation")
        Grid2D<double> g11(0.0, 1.0, 0.0, 2.0, 2, 3);
        auto [X, Y] = g11.meshgrid();

        ASSERT(X.rows() == 2);
        ASSERT(X.cols() == 3);
        ASSERT(Y.rows() == 2);
        ASSERT(Y.cols() == 3);

        // Check X values
        ASSERT_APPROX(X(0,0), 0.0);
        ASSERT_APPROX(X(0,1), 0.0);
        ASSERT_APPROX(X(0,2), 0.0);
        ASSERT_APPROX(X(1,0), 1.0);
        ASSERT_APPROX(X(1,1), 1.0);
        ASSERT_APPROX(X(1,2), 1.0);

        // Check Y values
        ASSERT_APPROX(Y(0,0), 0.0);
        ASSERT_APPROX(Y(0,1), 1.0);
        ASSERT_APPROX(Y(0,2), 2.0);
        ASSERT_APPROX(Y(1,0), 0.0);
        ASSERT_APPROX(Y(1,1), 1.0);
        ASSERT_APPROX(Y(1,2), 2.0);
    END_TEST

    // ============= Utility Function Tests =============

    TEST_CASE("coords function")
        Grid2D<double> g12(0.0, 2.0, 0.0, 3.0, 3, 4);
        auto [x1, y1] = g12.coords(0, 0);
        ASSERT_APPROX(x1, 0.0);
        ASSERT_APPROX(y1, 0.0);

        auto [x2, y2] = g12.coords(1, 2);
        ASSERT_APPROX(x2, 1.0);
        ASSERT_APPROX(y2, 2.0);

        auto [x3, y3] = g12.coords(2, 3);
        ASSERT_APPROX(x3, 2.0);
        ASSERT_APPROX(y3, 3.0);
    END_TEST

    TEST_CASE("coords out of bounds")
        Grid2D<double> g13(0.0, 1.0, 0.0, 1.0, 2, 2);
        bool caught = false;
        try {
            auto [x, y] = g13.coords(2, 1); // i out of bounds
        } catch (const std::out_of_range&) {
            caught = true;
        }
        ASSERT(caught);
    END_TEST

    TEST_CASE("contains function - inside domain")
        Grid2D<double> g14(0.0, 2.0, -1.0, 1.0, 5, 5);
        ASSERT(g14.contains(0.5, 0.0) == true);
        ASSERT(g14.contains(1.0, -0.5) == true);
        ASSERT(g14.contains(0.0, 0.0) == true);
        ASSERT(g14.contains(2.0, 1.0) == true);  // On boundary
    END_TEST

    TEST_CASE("contains function - outside domain")
        Grid2D<double> g15(0.0, 2.0, -1.0, 1.0, 5, 5);
        ASSERT(g15.contains(-0.1, 0.0) == false);
        ASSERT(g15.contains(2.1, 0.0) == false);
        ASSERT(g15.contains(1.0, -1.1) == false);
        ASSERT(g15.contains(1.0, 1.1) == false);
        ASSERT(g15.contains(3.0, 2.0) == false);
    END_TEST

    TEST_CASE("nearest function - point inside domain")
        Grid2D<double> g16(0.0, 4.0, 0.0, 3.0, 5, 4);
        // Grid x: 0, 1, 2, 3, 4
        // Grid y: 0, 1, 2, 3

        auto [i1, j1] = g16.nearest(0.4, 0.3);
        ASSERT(i1 == 0);  // Nearest to x=0
        ASSERT(j1 == 0);  // Nearest to y=0

        auto [i2, j2] = g16.nearest(1.6, 1.7);
        ASSERT(i2 == 2);  // Nearest to x=2
        ASSERT(j2 == 2);  // Nearest to y=2

        auto [i3, j3] = g16.nearest(3.8, 2.9);
        ASSERT(i3 == 4);  // Nearest to x=4
        ASSERT(j3 == 3);  // Nearest to y=3
    END_TEST

    TEST_CASE("nearest function - point outside domain")
        Grid2D<double> g17(0.0, 2.0, 0.0, 2.0, 3, 3);
        // Grid x: 0, 1, 2
        // Grid y: 0, 1, 2

        // Point outside gets clamped to boundary
        auto [i1, j1] = g17.nearest(-1.0, -1.0);
        ASSERT(i1 == 0);  // Clamped to x=0
        ASSERT(j1 == 0);  // Clamped to y=0

        auto [i2, j2] = g17.nearest(3.0, 3.0);
        ASSERT(i2 == 2);  // Clamped to x=2
        ASSERT(j2 == 2);  // Clamped to y=2

        auto [i3, j3] = g17.nearest(1.0, 3.0);
        ASSERT(i3 == 1);  // x=1 is in domain
        ASSERT(j3 == 2);  // Clamped to y=2
    END_TEST

    // ============= Edge Cases =============

    TEST_CASE("Minimum size grid (2x2)")
        Grid2D<double> g18(0.0, 1.0, 0.0, 1.0, 2, 2);
        ASSERT(g18.nx() == 2);
        ASSERT(g18.ny() == 2);
        ASSERT_APPROX(g18.dx(), 1.0);
        ASSERT_APPROX(g18.dy(), 1.0);

        auto X = g18.meshgridX();
        ASSERT(X.rows() == 2);
        ASSERT(X.cols() == 2);
    END_TEST

    TEST_CASE("Asymmetric grid")
        Grid2D<double> g19(-5.0, 5.0, 0.0, 1.0, 11, 2);
        ASSERT_APPROX(g19.dx(), 1.0);
        ASSERT_APPROX(g19.dy(), 1.0);
        ASSERT_APPROX(g19.x(5), 0.0);  // Middle point
    END_TEST

    TEST_CASE("Float type grid")
        Grid2D<float> g20(0.0f, 1.0f, 0.0f, 1.0f, 3, 3);
        ASSERT(g20.nx() == 3);
        ASSERT(g20.ny() == 3);
        float dx = g20.dx();
        float dy = g20.dy();
        ASSERT_APPROX(dx, 0.5f);
        ASSERT_APPROX(dy, 0.5f);
    END_TEST

    TEST_CASE("Large grid performance")
        Grid2D<double> g21(0.0, 100.0, 0.0, 100.0, 1000, 1000);
        ASSERT(g21.nx() == 1000);
        ASSERT(g21.ny() == 1000);

        // Test that large meshgrid generation doesn't crash
        auto [X, Y] = g21.meshgrid();
        ASSERT(X.rows() == 1000);
        ASSERT(X.cols() == 1000);
        ASSERT(Y.rows() == 1000);
        ASSERT(Y.cols() == 1000);

        // Check corners
        ASSERT_APPROX(X(0,0), 0.0);
        ASSERT_APPROX(Y(0,0), 0.0);
        ASSERT_APPROX(X(999,999), 100.0);
        ASSERT_APPROX(Y(999,999), 100.0);
    END_TEST

    // ============= Integration Test =============

    TEST_CASE("Complete workflow test")
        // Create a grid for heat equation domain
        Grid2D<double> grid(-1.0, 1.0, -1.0, 1.0, 21, 21);

        // Check grid properties
        ASSERT_APPROX(grid.dx(), 0.1);
        ASSERT_APPROX(grid.dy(), 0.1);

        // Generate meshgrids
        auto [X, Y] = grid.meshgrid();

        // Check that meshgrids have correct shape
        ASSERT(X.rows() == 21);
        ASSERT(X.cols() == 21);

        // Test point location
        auto [i, j] = grid.nearest(0.0, 0.0);
        ASSERT(i == 10);  // Center of 21-point grid
        ASSERT(j == 10);

        // Verify center coordinates
        auto [cx, cy] = grid.coords(10, 10);
        ASSERT_APPROX(cx, 0.0);
        ASSERT_APPROX(cy, 0.0);

        // Check domain containment
        ASSERT(grid.contains(0.0, 0.0) == true);
        ASSERT(grid.contains(-0.5, 0.5) == true);
        ASSERT(grid.contains(1.5, 0.0) == false);
    END_TEST

    // ============= Summary =============

    std::cout << "\n========================================\n";
    std::cout << "Test Results Summary:\n";
    std::cout << "  Passed: " << tests_passed << "\n";
    std::cout << "  Failed: " << tests_failed << "\n";

    if (tests_failed == 0) {
        std::cout << "\n✓ All tests PASSED! Grid2D class is working correctly.\n";
    } else {
        std::cout << "\n✗ Some tests FAILED. Please review the implementation.\n";
    }
    std::cout << "========================================\n";

    return (tests_failed == 0) ? 0 : 1;
}