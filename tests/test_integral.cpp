/**
 * Integration module tests
 */

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>

#include "integral.hpp"

using namespace msl;
using namespace msl::integral;

void test_cumtrapz_vector()
{
    std::cout << "Test 1: Cumulative Trapezoidal (Vector)\n";
    std::cout << "========================================\n";

    // Integrate y = x from 0 to 4 (result should be x²/2)
    std::vector<double> x = {0, 1, 2, 3, 4};
    std::vector<double> y = {0, 1, 2, 3, 4};

    auto result = cumtrapz(x, y);

    std::cout << "Integrating y=x:\n";
    for (size_t i = 0; i < x.size(); ++i)
    {
        double exact = 0.5 * x[i] * x[i];
        double error = std::abs(result[i] - exact);
        std::cout << "  x=" << x[i] << ": integral=" << std::setprecision(4)
                  << result[i] << ", exact=" << exact
                  << ", error=" << std::scientific << error << "\n";
    }
    std::cout << "\n";
}

void test_cumtrapz_matrix()
{
    std::cout << "Test 2: Cumulative Trapezoidal (Matrix)\n";
    std::cout << "========================================\n";

    // Matrix with 3 columns: x, 2x, x²
    matrixd mat(5, 3);
    std::vector<double> x = {0, 1, 2, 3, 4};

    for (size_t i = 0; i < 5; ++i)
    {
        mat(i, 0) = x[i];        // x
        mat(i, 1) = 2.0 * x[i];  // 2x
        mat(i, 2) = x[i] * x[i]; // x²
    }

    auto result = cumtrapz(x, mat);

    std::cout << "Integrating [x, 2x, x²]:\n";
    std::cout << "At x=4:\n";
    std::cout << "  Column 0 (x):   " << result(4, 0) << " (exact: 8)\n";
    std::cout << "  Column 1 (2x):  " << result(4, 1) << " (exact: 16)\n";
    std::cout << "  Column 2 (x²):  " << result(4, 2) << " (exact: 21.33)\n\n";
}

void test_trapz_total()
{
    std::cout << "Test 3: Total Integral (Trapezoidal)\n";
    std::cout << "=====================================\n";

    // Integrate sin(x) from 0 to π (exact: 2)
    const size_t n = 101;
    std::vector<double> x(n), y(n);

    for (size_t i = 0; i < n; ++i)
    {
        x[i] = i * std::numbers::pi / (n - 1);
        y[i] = std::sin(x[i]);
    }

    double result = trapz(x, y);
    double exact = 2.0;
    double error = std::abs(result - exact);

    std::cout << "∫sin(x)dx from 0 to π:\n";
    std::cout << "  Result: " << std::setprecision(8) << result << "\n";
    std::cout << "  Exact:  " << exact << "\n";
    std::cout << "  Error:  " << std::scientific << error << "\n\n";
}

void test_simps()
{
    std::cout << "Test 4: Simpson's Rule\n";
    std::cout << "======================\n";

    // Integrate x³ from 0 to 2 (exact: 4)
    std::vector<double> x = {0, 0.5, 1.0, 1.5, 2.0};
    std::vector<double> y(x.size());

    for (size_t i = 0; i < x.size(); ++i)
    {
        y[i] = x[i] * x[i] * x[i];
    }

    double trap_result = trapz(x, y);
    double simp_result = simpson(y, 0.5);
    double exact = 4.0;

    std::cout << "∫x³dx from 0 to 2:\n";
    std::cout << "  Trapezoidal: " << trap_result
              << " (error: " << std::abs(trap_result - exact) << ")\n";
    std::cout << "  Simpson's:   " << simp_result
              << " (error: " << std::abs(simp_result - exact) << ")\n";
    std::cout << "  Exact:       " << exact << "\n\n";
}

void test_matrix_integration()
{
    std::cout << "Test 5: Matrix Column Integration\n";
    std::cout << "==================================\n";

    // Create matrix with sine waves at different frequencies
    const size_t n = 100;
    const size_t n_cols = 3;
    matrixd mat(n, n_cols);

    std::vector<double> freqs = {1.0, 2.0, 3.0};
    double dx = 2.0 * std::numbers::pi / (n - 1);

    for (size_t j = 0; j < n_cols; ++j)
    {
        for (size_t i = 0; i < n; ++i)
        {
            double x = i * dx;
            mat(i, j) = std::sin(freqs[j] * x);
        }
    }

    auto results = trapz(mat, dx);

    std::cout << "Integrating sin(fx) from 0 to 2π:\n";
    for (size_t j = 0; j < n_cols; ++j)
    {
        std::cout << "  f=" << freqs[j] << ": " << results[j]
                  << " (exact: ~0)\n";
    }
    std::cout << "\n";
}

void test_adaptive_quad()
{
    std::cout << "Test 6: Adaptive Quadrature\n";
    std::cout << "===========================\n";

    // Test with difficult function
    auto f1 = [](double x) { return std::exp(-x * x); };       // Gaussian
    auto f2 = [](double x) { return 1.0 / (1.0 + x * x); };    // Lorentzian
    auto f3 = [](double x) { return std::sin(10.0 * x) / x; }; // Oscillatory

    std::cout << "∫exp(-x²)dx from 0 to 3:\n";
    double r1 = quad(f1, 0.0, 3.0);
    std::cout << "  Result: " << r1 << " (exact: ~0.886)\n\n";

    std::cout << "∫1/(1+x²)dx from 0 to 1:\n";
    double r2 = quad(f2, 0.0, 1.0);
    double exact2 = std::numbers::pi / 4.0;
    std::cout << "  Result: " << r2 << "\n";
    std::cout << "  Exact:  " << exact2 << "\n";
    std::cout << "  Error:  " << std::abs(r2 - exact2) << "\n\n";
}

void test_integral_comparison()
{
    std::cout << "Test 7: Method Comparison\n";
    std::cout << "=========================\n";

    // Compare accuracy for e^x from 0 to 1 (exact: e - 1)
    const double exact = std::numbers::e - 1.0;

    std::vector<size_t> n_points = {5, 11, 21, 51};

    std::cout << std::setw(10) << "N points" << std::setw(15) << "Trapezoidal"
              << std::setw(15) << "Simpson's" << std::setw(15) << "Romberg\n";
    std::cout << std::string(55, '-') << "\n";

    for (auto n : n_points)
    {
        std::vector<double> y(n);
        double dx = 1.0 / (n - 1);

        for (size_t i = 0; i < n; ++i)
        {
            y[i] = std::exp(i * dx);
        }

        double trap = trapz(y, dx);

        std::cout << std::setw(10) << n << std::setw(15) << std::scientific
                  << std::abs(trap - exact);

        if (n >= 3)
        {
            double simp = simpson(y, dx);
            std::cout << std::setw(15) << std::abs(simp - exact);
        }
        else
        {
            std::cout << std::setw(15) << "-";
        }

        // Check if n = 2^k + 1
        size_t check = 1;
        while (check < n)
            check *= 2;
        if (check + 1 == n)
        {
            double romb = romberg(y, dx);
            std::cout << std::setw(15) << std::abs(romb - exact);
        }
        else
        {
            std::cout << std::setw(15) << "-";
        }

        std::cout << "\n";
    }
    std::cout << "\n";
}

int test_integral()
{
    try
    {
        test_cumtrapz_vector();
        test_cumtrapz_matrix();
        test_trapz_total();
        test_simps();
        test_matrix_integration();
        test_adaptive_quad();
        test_integral_comparison();

        std::cout << "All integration tests completed!\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}