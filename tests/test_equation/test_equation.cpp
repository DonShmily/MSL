#include <cmath>
#include <iostream>
#include <stdexcept>

#include "equation.hpp"

using namespace msl;

static int g_failures = 0;
static int g_total = 0;

#define EXPECT_TRUE(cond)                                                      \
    do {                                                                       \
        ++g_total;                                                             \
        if (!(cond)) {                                                         \
            ++g_failures;                                                      \
            std::cerr << "[FAIL] " << __FILE__ << ":" << __LINE__ << " - "     \
                      << #cond << "\n";                                        \
        }                                                                      \
    } while (0)

#define EXPECT_EQ(a, b) EXPECT_TRUE((a) == (b))

#define EXPECT_NEAR(a, b, eps)                                                 \
    do {                                                                       \
        ++g_total;                                                             \
        if (std::fabs((a) - (b)) > (eps)) {                                    \
            ++g_failures;                                                      \
            std::cerr << "[FAIL] " << __FILE__ << ":" << __LINE__ << " - "     \
                      << #a << " ~= " << #b << " (got: " << (a) << " vs "      \
                      << (b) << ")\n";                                         \
        }                                                                      \
    } while (0)

int test_equation() {
    auto square_minus_two = [](double x) { return x * x - 2.0; };
    auto square_minus_two_deriv = [](double x) { return 2.0 * x; };
    double sqrt2 = std::sqrt(2.0);

    auto bisection_result =
        equation::bisection(square_minus_two, 0.0, 2.0, 1e-12);
    EXPECT_TRUE(bisection_result);
    EXPECT_EQ(bisection_result.status, equation::root_status::converged);
    EXPECT_NEAR(bisection_result.root, sqrt2, 1e-10);
    EXPECT_NEAR(bisection_result.residual, 0.0, 1e-10);
    EXPECT_NEAR(bisection_result.value(), sqrt2, 1e-10);
    EXPECT_NEAR(static_cast<double>(bisection_result), sqrt2, 1e-10);

    double bisection_root =
        equation::bisection_root(square_minus_two, 0.0, 2.0, 1e-12);
    EXPECT_NEAR(bisection_root, sqrt2, 1e-10);

    auto newton_result =
        equation::newton(square_minus_two, square_minus_two_deriv, 1.0, 1e-12);
    EXPECT_TRUE(newton_result.converged());
    EXPECT_NEAR(newton_result.root, sqrt2, 1e-12);
    EXPECT_NEAR(equation::newton_root(
                    square_minus_two, square_minus_two_deriv, 1.0, 1e-12),
                sqrt2,
                1e-12);

    auto secant_result = equation::secant(square_minus_two, 1.0, 2.0, 1e-12);
    EXPECT_TRUE(secant_result.converged());
    EXPECT_NEAR(secant_result.root, sqrt2, 1e-10);
    EXPECT_NEAR(
        equation::secant_root(square_minus_two, 1.0, 2.0, 1e-12), sqrt2, 1e-10);

    auto falsi_result =
        equation::regula_falsi(square_minus_two, 0.0, 2.0, 1e-12);
    EXPECT_TRUE(falsi_result.converged());
    EXPECT_NEAR(falsi_result.root, sqrt2, 1e-10);
    EXPECT_NEAR(equation::regula_falsi_root(square_minus_two, 0.0, 2.0, 1e-12),
                sqrt2,
                1e-10);

    auto fixed_point = equation::bisection(
        [](double x) { return std::cos(x) - x; }, 0.0, 1.0, 1e-12);
    EXPECT_TRUE(fixed_point.converged());
    EXPECT_NEAR(fixed_point.root, 0.7390851332151607, 1e-10);

    auto invalid_interval =
        equation::bisection([](double x) { return x * x + 1.0; }, -1.0, 1.0);
    EXPECT_EQ(invalid_interval.status, equation::root_status::invalid_interval);
    bool value_threw = false;
    try {
        (void)invalid_interval.value();
    } catch (const std::runtime_error &) {
        value_threw = true;
    }
    EXPECT_TRUE(value_threw);

    auto zero_derivative = equation::newton(
        [](double x) { return x * x + 1.0; }, [](double) { return 0.0; }, 1.0);
    EXPECT_EQ(zero_derivative.status, equation::root_status::zero_derivative);

    auto max_iterations =
        equation::secant(square_minus_two, 1.0, 2.0, 1e-16, 1);
    EXPECT_EQ(max_iterations.status, equation::root_status::max_iterations);

    bool threw = false;
    try {
        (void)equation::bisection(square_minus_two, 2.0, 0.0);
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE(threw);

    std::cout << "Total checks: " << g_total << ", failures: " << g_failures
              << "\n";
    return g_failures == 0 ? 0 : 1;
}

int main() { return test_equation(); }
