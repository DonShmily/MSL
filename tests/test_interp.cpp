/**
 * Simple interpolation tests
 */

#include <cmath>
#include <iostream>
#include "interp/interp_1d_akima.hpp"
#include "interp/interp_1d_linear.hpp"
#include "interp/interp_1d_polynomial.hpp"
#include "interp/interp_1d_spline.hpp"


using namespace msl::interp;

void test_linear()
{
    std::cout << "Linear Interpolation Test\n";
    std::cout << "--------------------------\n";

    // Data points
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> y = {0.0, 1.0, 4.0, 9.0, 16.0}; // x^2

    Linear interp(x, y);

    // Test at known points
    std::cout << "At x=0.0: " << interp(0.0) << " (expected: 0.0)\n";
    std::cout << "At x=1.5: " << interp(1.5) << " (expected: ~2.5)\n";
    std::cout << "At x=2.5: " << interp(2.5) << " (expected: ~6.5)\n\n";
}

void test_cubic()
{
    std::cout << "Cubic Spline Test\n";
    std::cout << "-----------------\n";

    // Sine wave data
    std::vector<double> x, y;
    for (int i = 0; i <= 10; ++i)
    {
        double xi = i * 0.5;
        x.push_back(xi);
        y.push_back(std::sin(xi));
    }

    CubicSpline spline(x, y);
    spline.set_extrapolation(msl::interp::ExtrapolationMode::Linear);

    // Interpolate at finer grid
    std::cout << "Interpolating sin(x):\n";
    for (double xi = 0; xi <= 6; xi += 1.0)
    {
        double y_interp = spline(xi);
        double y_exact = std::sin(xi);
        double error = std::abs(y_interp - y_exact);
        std::cout << "  x=" << xi << ": interp=" << y_interp
                  << ", exact=" << y_exact << ", error=" << error << "\n";
    }
    std::cout << "\n";
}

void test_interp_comparison()
{
    std::cout << "Method Comparison (Runge function)\n";
    std::cout << "-----------------------------------\n";

    // Runge function: 1 / (1 + 25*x^2)
    auto runge = [](double x) { return 1.0 / (1.0 + 25.0 * x * x); };

    // Sample points
    std::vector<double> x_data, y_data;
    for (int i = 0; i <= 10; ++i)
    {
        double x = -1.0 + 2.0 * i / 10.0;
        x_data.push_back(x);
        y_data.push_back(runge(x));
    }

    // Test point
    double x_test = 0.5;
    double y_exact = runge(x_test);

    Linear linear(x_data, y_data);
    CubicSpline cubic(x_data, y_data);
    AkimaSpline akima(x_data, y_data);
    PolynomialInterpolator poly(x_data, y_data);

    std::cout << "At x=" << x_test << ":\n";
    std::cout << "  Exact:  " << y_exact << "\n";
    std::cout << "  Linear: " << linear(x_test)
              << " (error: " << std::abs(linear(x_test) - y_exact) << ")\n";
    std::cout << "  Cubic:  " << cubic(x_test)
              << " (error: " << std::abs(cubic(x_test) - y_exact) << ")\n";
    std::cout << "  Akima:  " << akima(x_test)
              << " (error: " << std::abs(akima(x_test) - y_exact) << ")\n";
    std::cout << "  Polynomial: " << poly(x_test)
              << " (error: " << std::abs(poly(x_test) - y_exact) << ")\n\n";
}

int test_interp()
{
    try
    {
        test_linear();
        test_cubic();
        test_interp_comparison();

        std::cout << "All tests completed!\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
