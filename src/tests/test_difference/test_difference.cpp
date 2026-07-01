#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "difference.hpp"
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

std::filesystem::path project_root() {
    auto path = std::filesystem::current_path();
    while (!path.empty()) {
        if (std::filesystem::exists(path / "msl")
            && std::filesystem::exists(path / "tests")) {
            return path;
        }
        auto parent = path.parent_path();
        if (parent == path) {
            break;
        }
        path = parent;
    }
    return std::filesystem::current_path();
}

int test_difference() {
    std::vector<double> y{1.0, 4.0, 9.0, 16.0};
    auto d = difference::diff(y);
    EXPECT_EQ(d.size(), 3);
    EXPECT_NEAR(d[0], 3.0, 1e-12);
    EXPECT_NEAR(d[1], 5.0, 1e-12);
    EXPECT_NEAR(d[2], 7.0, 1e-12);

    auto fg = difference::forward_gradient(y, 0.5);
    EXPECT_EQ(fg.size(), 3);
    EXPECT_NEAR(fg[0], 6.0, 1e-12);
    EXPECT_NEAR(fg[1], 10.0, 1e-12);
    EXPECT_NEAR(fg[2], 14.0, 1e-12);

    auto cg = difference::central_gradient(y, 1.0);
    EXPECT_EQ(cg.size(), y.size());
    EXPECT_NEAR(cg[0], 3.0, 1e-12);
    EXPECT_NEAR(cg[1], 4.0, 1e-12);
    EXPECT_NEAR(cg[2], 6.0, 1e-12);
    EXPECT_NEAR(cg[3], 7.0, 1e-12);

    matrix::matrixd mat(3, 3);
    for (size_t i = 0; i < mat.rows(); ++i) {
        for (size_t j = 0; j < mat.cols(); ++j) {
            mat(i, j) = 10.0 * static_cast<double>(i) + static_cast<double>(j);
        }
    }

    auto row_diff = difference::diff(mat, 0);
    EXPECT_EQ(row_diff.rows(), 2);
    EXPECT_EQ(row_diff.cols(), 3);
    for (size_t i = 0; i < row_diff.rows(); ++i) {
        for (size_t j = 0; j < row_diff.cols(); ++j) {
            EXPECT_NEAR(row_diff(i, j), 10.0, 1e-12);
        }
    }

    auto col_diff = difference::diff(mat, 1);
    EXPECT_EQ(col_diff.rows(), 3);
    EXPECT_EQ(col_diff.cols(), 2);
    for (size_t i = 0; i < col_diff.rows(); ++i) {
        for (size_t j = 0; j < col_diff.cols(); ++j) {
            EXPECT_NEAR(col_diff(i, j), 1.0, 1e-12);
        }
    }

    auto lap = difference::laplacian(mat);
    EXPECT_EQ(lap.rows(), 3);
    EXPECT_EQ(lap.cols(), 3);
    EXPECT_NEAR(lap(1, 1), 0.0, 1e-12);

    bool threw = false;
    try {
        (void)difference::diff(std::vector<double>{1.0});
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE(threw);

    auto compare_dir =
        project_root() / "test_result" / "difference" / "matlab_compare";
    std::filesystem::create_directories(compare_dir);

    std::ofstream vec_file(compare_dir / "difference_vector.txt");
    vec_file << std::setprecision(17);
    for (size_t i = 0; i < y.size(); ++i) {
        vec_file << y[i] << " " << (i < d.size() ? d[i] : 0.0) << " " << cg[i]
                 << "\n";
    }

    std::ofstream mat_file(compare_dir / "difference_matrix.txt");
    mat_file << std::setprecision(17);
    for (size_t i = 0; i < mat.rows(); ++i) {
        for (size_t j = 0; j < mat.cols(); ++j) {
            mat_file << mat(i, j) << (j + 1 == mat.cols() ? '\n' : ' ');
        }
    }

    std::ofstream row_file(compare_dir / "difference_row_diff.txt");
    row_file << std::setprecision(17);
    for (size_t i = 0; i < row_diff.rows(); ++i) {
        for (size_t j = 0; j < row_diff.cols(); ++j) {
            row_file << row_diff(i, j)
                     << (j + 1 == row_diff.cols() ? '\n' : ' ');
        }
    }

    std::ofstream col_file(compare_dir / "difference_col_diff.txt");
    col_file << std::setprecision(17);
    for (size_t i = 0; i < col_diff.rows(); ++i) {
        for (size_t j = 0; j < col_diff.cols(); ++j) {
            col_file << col_diff(i, j)
                     << (j + 1 == col_diff.cols() ? '\n' : ' ');
        }
    }

    std::cout << "Total checks: " << g_total << ", failures: " << g_failures
              << "\n";
    return g_failures == 0 ? 0 : 1;
}

int main() { return test_difference(); }
