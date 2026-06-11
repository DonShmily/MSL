#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "integral.hpp"
#include "matrix.hpp"

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

int test_integral() {
    std::vector<double> y{0.0, 1.0, 4.0, 9.0};
    auto ct = integral::cumtrapz(y, 1.0);
    EXPECT_EQ(ct.size(), y.size());
    EXPECT_NEAR(ct[0], 0.0, 1e-12);
    EXPECT_NEAR(ct[1], 0.5, 1e-12);
    EXPECT_NEAR(ct[2], 3.0, 1e-12);
    EXPECT_NEAR(ct[3], 9.5, 1e-12);
    EXPECT_NEAR(integral::trapz(y, 1.0), 9.5, 1e-12);

    std::vector<double> x{0.0, 0.5, 2.0, 3.0};
    auto ct_nonuniform = integral::cumtrapz(x, y);
    EXPECT_NEAR(ct_nonuniform[1], 0.25, 1e-12);
    EXPECT_NEAR(ct_nonuniform[2], 4.0, 1e-12);
    EXPECT_NEAR(ct_nonuniform[3], 10.5, 1e-12);
    EXPECT_NEAR(integral::trapz(x, y), 10.5, 1e-12);

    std::vector<double> parabola{0.0, 1.0, 4.0};
    EXPECT_NEAR(integral::simpson(parabola, 1.0), 8.0 / 3.0, 1e-12);
    auto cs = integral::cumsimpson(parabola, 1.0);
    EXPECT_EQ(cs.size(), parabola.size());
    EXPECT_NEAR(cs[0], 0.0, 1e-12);
    EXPECT_NEAR(cs[2], 8.0 / 3.0, 1e-12);

    matrix::matrixd mat(4, 2);
    for (size_t i = 0; i < mat.rows(); ++i) {
        mat(i, 0) = static_cast<double>(i);
        mat(i, 1) = 2.0 * static_cast<double>(i);
    }
    auto col_trapz = integral::trapz(mat, 1.0);
    EXPECT_EQ(col_trapz.size(), 2);
    EXPECT_NEAR(col_trapz[0], 4.5, 1e-12);
    EXPECT_NEAR(col_trapz[1], 9.0, 1e-12);

    auto mat_cum = integral::cumtrapz(mat, 1.0);
    EXPECT_EQ(mat_cum.rows(), 4);
    EXPECT_EQ(mat_cum.cols(), 2);
    EXPECT_NEAR(mat_cum(3, 0), 4.5, 1e-12);
    EXPECT_NEAR(mat_cum(3, 1), 9.0, 1e-12);

    bool threw = false;
    try {
        (void)integral::trapz(std::vector<double>{1.0}, 1.0);
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE(threw);

    std::cout << "Total checks: " << g_total << ", failures: " << g_failures
              << "\n";
    return g_failures == 0 ? 0 : 1;
}

int main() { return test_integral(); }
