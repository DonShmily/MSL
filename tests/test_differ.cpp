/**
 * Differentiation module tests
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include "differ/differ.hpp"


using namespace msl;
using namespace msl::differ;

void test_diff()
{
    std::cout << "Test 1: First Difference\n";
    std::cout << "=========================\n";

    std::vector<double> y = {1, 4, 9, 16, 25}; // x²
    auto dy = diff(y);

    std::cout << "y = [1, 4, 9, 16, 25]  (x²)\n";
    std::cout << "diff(y) = [";
    for (auto val : dy)
        std::cout << val << " ";
    std::cout << "]\nExpected: [3, 5, 7, 9]\n\n";
}

void test_gradient_uniform()
{
    std::cout << "Test 2: Gradient (Uniform Spacing)\n";
    std::cout << "===================================\n";

    // Test with sin(x)
    const size_t n = 11;
    std::vector<double> x(n), y(n);
    double dx = 0.1;

    for (size_t i = 0; i < n; ++i)
    {
        x[i] = i * dx;
        y[i] = std::sin(x[i]);
    }

    auto grad = central_gradient(y, dx);

    std::cout << "Function: y = sin(x)\n";
    std::cout << "Derivative should be: cos(x)\n\n";

    std::cout << std::setw(8) << "x" << std::setw(12) << "sin(x)"
              << std::setw(12) << "dy/dx" << std::setw(12) << "cos(x)"
              << std::setw(12) << "Error\n";
    std::cout << std::string(56, '-') << "\n";

    for (size_t i = 0; i < n; ++i)
    {
        double exact = std::cos(x[i]);
        double error = std::abs(grad[i] - exact);
        std::cout << std::setw(8) << std::setprecision(3) << x[i]
                  << std::setw(12) << std::setprecision(6) << y[i]
                  << std::setw(12) << grad[i] << std::setw(12) << exact
                  << std::setw(12) << std::scientific << error << "\n";
    }
    std::cout << "\n";
}

void test_gradient2()
{
    std::cout << "Test 3: Second Derivative\n";
    std::cout << "==========================\n";

    // Test with y = x²
    std::vector<double> x, y;
    double dx = 0.5;

    for (double xi = 0; xi <= 5.0; xi += dx)
    {
        x.push_back(xi);
        y.push_back(xi * xi);
    }

    auto d2y = central_gradient2(y, dx);

    std::cout << "Function: y = x²\n";
    std::cout << "Second derivative should be: 2\n\n";

    std::cout << "Sample results:\n";
    for (size_t i = 0; i < std::min(size_t(6), d2y.size()); ++i)
    {
        std::cout << "  x=" << x[i] << ": d²y/dx² = " << d2y[i] << "\n";
    }
    std::cout << "\n";
}

void test_matrix_diff()
{
    std::cout << "Test 4: Matrix Difference\n";
    std::cout << "==========================\n";

    // Create 4x3 matrix
    matrixd mat(4, 3);
    for (size_t i = 0; i < 4; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            mat(i, j) = (i + 1) * 10 + (j + 1);
        }
    }

    std::cout << "Original matrix (4x3):\n";
    for (size_t i = 0; i < 4; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            std::cout << std::setw(5) << mat(i, j);
        }
        std::cout << "\n";
    }

    auto diff_rows = diff(mat, 0); // Vertical
    std::cout << "\nRow-wise diff (3x3):\n";
    for (size_t i = 0; i < diff_rows.rows(); ++i)
    {
        for (size_t j = 0; j < diff_rows.cols(); ++j)
        {
            std::cout << std::setw(5) << diff_rows(i, j);
        }
        std::cout << "\n";
    }

    auto diff_cols = diff(mat, 1); // Horizontal
    std::cout << "\nColumn-wise diff (4x2):\n";
    for (size_t i = 0; i < diff_cols.rows(); ++i)
    {
        for (size_t j = 0; j < diff_cols.cols(); ++j)
        {
            std::cout << std::setw(5) << diff_cols(i, j);
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void test_matrix_gradient()
{
    std::cout << "Test 5: Matrix Gradient (2D)\n";
    std::cout << "=============================\n";

    // Create a 2D function: z = x² + y²
    const size_t nx = 5, ny = 5;
    matrixd z(ny, nx);
    double dx = 0.5, dy = 0.5;

    for (size_t i = 0; i < ny; ++i)
    {
        for (size_t j = 0; j < nx; ++j)
        {
            double x = j * dx;
            double y = i * dy;
            z(i, j) = x * x + y * y;
        }
    }

    // Compute 2D gradient
    auto [grad_x, grad_y] = central_gradient2d(z, dx, dy);

    std::cout << "Function: z = x² + y²\n";
    std::cout << "∂z/∂x should be: 2x\n";
    std::cout << "∂z/∂y should be: 2y\n\n";

    std::cout << "At center point (x=1, y=1):\n";
    size_t i_center = 2, j_center = 2;
    double x_center = j_center * dx;
    double y_center = i_center * dy;

    std::cout << "  ∂z/∂x = " << grad_x(i_center, j_center)
              << " (expected: " << 2 * x_center << ")\n";
    std::cout << "  ∂z/∂y = " << grad_y(i_center, j_center)
              << " (expected: " << 2 * y_center << ")\n\n";
}

void test_laplacian()
{
    std::cout << "Test 6: Laplacian\n";
    std::cout << "==================\n";

    // Test with f = x² + y² (∇²f = 4)
    const size_t n = 7;
    matrixd f(n, n);
    double dx = 0.5, dy = 0.5;

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            double x = j * dx - 1.5;
            double y = i * dy - 1.5;
            f(i, j) = x * x + y * y;
        }
    }

    auto lap = laplacian(f, dx, dy);

    std::cout << "Function: f = x² + y²\n";
    std::cout << "Laplacian should be: 4\n\n";

    std::cout << "Laplacian at center points:\n";
    for (size_t i = 2; i < 5; ++i)
    {
        for (size_t j = 2; j < 5; ++j)
        {
            std::cout << "  (" << i << "," << j << "): " << std::setprecision(4)
                      << lap(i, j) << "\n";
        }
    }
    std::cout << "\n";
}

void test_vector_field()
{
    std::cout << "Test 7: Divergence & Curl\n";
    std::cout << "==========================\n";

    // Test with simple vector field
    const size_t n = 5;
    matrixd Fx(n, n), Fy(n, n);
    double dx = 0.5, dy = 0.5;

    // Field: F = (x, y) - should have div = 2, curl = 0
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            double x = j * dx;
            double y = i * dy;
            Fx(i, j) = x;
            Fy(i, j) = y;
        }
    }

    auto div_F = divergence(Fx, Fy, dx, dy);
    auto curl_F = curl(Fx, Fy, dx, dy);

    std::cout << "Vector field: F = (x, y)\n";
    std::cout << "Expected: div(F) = 2, curl(F) = 0\n\n";

    std::cout << "At center point:\n";
    size_t center = n / 2;
    std::cout << "  div(F) = " << div_F(center, center) << "\n";
    std::cout << "  curl(F) = " << curl_F(center, center) << "\n\n";
}

void test_accuracy_comparison()
{
    std::cout << "Test 8: Accuracy vs Step Size\n";
    std::cout << "==============================\n";

    // Test gradient accuracy for different step sizes
    auto f = [](double x) { return std::exp(x); };
    auto df = [](double x) { return std::exp(x); };

    double x0 = 1.0;
    std::vector<double> step_sizes = {0.1, 0.01, 0.001, 0.0001};

    std::cout << "Function: f(x) = e^x at x = 1\n";
    std::cout << "Exact derivative: " << df(x0) << "\n\n";

    std::cout << std::setw(12) << "Step size" << std::setw(15) << "Gradient"
              << std::setw(15) << "Error\n";
    std::cout << std::string(42, '-') << "\n";

    for (double h : step_sizes)
    {
        std::vector<double> x = {x0, x0 + h};
        std::vector<double> y = {f(x0), f(x0 + h)};

        auto grad = forward_gradient(y, x);
        double error = std::abs(grad[1] - df(x0));

        std::cout << std::setw(12) << h << std::setw(15) << grad[1]
                  << std::setw(15) << std::scientific << error << "\n";
    }

    std::cout << "\n";
}

int test_differ()
{
    test_diff();
    test_gradient_uniform();
    test_gradient2();
    test_matrix_diff();
    test_matrix_gradient();
    test_laplacian();
    test_vector_field();
    test_accuracy_comparison();

    return 0;
}
