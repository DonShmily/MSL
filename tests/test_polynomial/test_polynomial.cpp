#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "polynomial.hpp"
#include "polynomial/polynomial.hpp"

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

int test_polynomial() {
    polynomial::Polynomial direct(std::vector<double>{1.0, 2.0, 3.0});
    EXPECT_EQ(direct.degree(), 2);
    EXPECT_NEAR(direct(2.0), 17.0, 1e-12);
    EXPECT_NEAR(direct.derivative(2.0), 14.0, 1e-12);

    std::vector<double> query{0.0, 1.0, 2.0};
    auto values = direct(query);
    EXPECT_EQ(values.size(), query.size());
    EXPECT_NEAR(values[0], 1.0, 1e-12);
    EXPECT_NEAR(values[1], 6.0, 1e-12);
    EXPECT_NEAR(values[2], 17.0, 1e-12);

    std::vector<double> out(query.size());
    polynomial::polyval(std::vector<double>{1.0, 2.0, 3.0}, query, out);
    EXPECT_NEAR(out[2], 17.0, 1e-12);

    std::vector<double> x{-2.0, -1.0, 0.0, 1.0, 2.0};
    std::vector<double> y;
    for (double xi : x) {
        y.push_back(1.0 - 2.0 * xi + 0.5 * xi * xi);
    }

    auto coeffs = polynomial::polyfit(x, y, 2);
    EXPECT_EQ(coeffs.size(), 3);
    EXPECT_NEAR(coeffs[0], 1.0, 1e-10);
    EXPECT_NEAR(coeffs[1], -2.0, 1e-10);
    EXPECT_NEAR(coeffs[2], 0.5, 1e-10);

    std::vector<double> coeffs_out(3);
    polynomial::polyfit(x, y, coeffs_out, 2);
    EXPECT_NEAR(coeffs_out[0], 1.0, 1e-10);
    EXPECT_NEAR(coeffs_out[1], -2.0, 1e-10);
    EXPECT_NEAR(coeffs_out[2], 0.5, 1e-10);

    auto constant_fit = polynomial::Polynomial::from_fit(y, 0);
    EXPECT_EQ(constant_fit.degree(), 0);

    bool threw = false;
    try {
        polynomial::Polynomial bad;
        (void)bad.degree();
    } catch (const std::runtime_error &) {
        threw = true;
    }
    EXPECT_TRUE(threw);

    std::cout << "Total checks: " << g_total << ", failures: " << g_failures
              << "\n";
    return g_failures == 0 ? 0 : 1;
}

int main() { return test_polynomial(); }
