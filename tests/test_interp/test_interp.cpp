#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "interp.hpp"

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

static void expect_vector_near(const std::vector<double> &actual,
                               const std::vector<double> &expected,
                               double eps) {
    EXPECT_EQ(actual.size(), expected.size());
    for (size_t i = 0; i < actual.size() && i < expected.size(); ++i) {
        EXPECT_NEAR(actual[i], expected[i], eps);
    }
}

int test_interp() {
    std::vector<double> x{0.0, 1.0, 2.0, 3.0};
    std::vector<double> y{1.0, 3.0, 5.0, 7.0};
    std::vector<double> xq{-1.0, 0.5, 1.5, 3.0, 4.0};
    std::vector<double> linear_expected{-1.0, 2.0, 4.0, 7.0, 9.0};

    expect_vector_near(interp::interp1_linear(x, y, xq),
                       linear_expected,
                       1e-12);
    expect_vector_near(interp::interp1_cubic(x, y, xq),
                       linear_expected,
                       1e-10);
    expect_vector_near(interp::interp1_pchip(x, y, xq),
                       linear_expected,
                       1e-12);
    expect_vector_near(interp::interp1_akima(x, y, xq),
                       linear_expected,
                       1e-12);

    std::vector<double> near_xq{0.25, 1.4, 2.6};
    std::vector<double> near_expected{1.0, 3.0, 7.0};
    expect_vector_near(interp::interp1_near(x, y, near_xq),
                       near_expected,
                       1e-12);

    std::vector<double> out(xq.size());
    interp::interp1_linear(x, y, xq, out);
    expect_vector_near(out, linear_expected, 1e-12);

    bool threw = false;
    try {
        (void)interp::interp1_linear(
            std::vector<double>{0.0, 0.0},
            std::vector<double>{1.0, 2.0},
            std::vector<double>{0.5});
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE(threw);

    std::cout << "Total checks: " << g_total << ", failures: " << g_failures
              << "\n";
    return g_failures == 0 ? 0 : 1;
}

int main() { return test_interp(); }
