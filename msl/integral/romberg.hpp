/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: romberg.hpp
** -----
** File Created: Friday, 9th January 2026 14:58:10
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Thursday, 5th March 2026 14:57:54
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_ROMBERG_HPP
#define MSL_ROMBERG_HPP

#include <cmath>
#include <span>
#include <stdexcept>
#include <vector>

#include "matrix/real_matrix_base.hpp"
#include "trapz.hpp"

namespace msl::integral {

// ============================================================================
// Total Integral (Romberg Integration)
// ============================================================================

/**
 * @brief Romberg integration (adaptive, high accuracy)
 *
 * Uses Richardson extrapolation on trapezoidal rule for high-precision total
 * integral computation.
 *
 * @param y Function values (size must be 2^k + 1)
 * @param dx Uniform spacing
 * @param tol Tolerance for convergence (default: 1e-10)
 * @return Integral value
 */
inline double
romberg(std::span<const double> y, double dx, double tol = 1e-10) {
    size_t n = y.size();

    // Check if n = 2^k + 1
    size_t k = 0;
    size_t check = 1;
    while (check < n) {
        check *= 2;
        ++k;
    }
    if (check + 1 != n) {
        throw std::invalid_argument("Romberg: needs 2^k + 1 points");
    }

    // Romberg table
    std::vector<std::vector<double>> R(k + 1);
    for (size_t i = 0; i <= k; ++i) {
        R[i].resize(i + 1);
    }

    // R[0,0] = trapezoidal with all points
    R[0][0] = trapz(y, dx);

    // Fill Romberg table
    size_t stride = 1;
    for (size_t i = 1; i <= k; ++i) {
        stride *= 2;

        // Trapezoidal with stride
        double sum = 0.5 * (y.front() + y.back());
        for (size_t j = stride; j < n; j += stride) {
            sum += y[j];
        }
        R[i][0] = sum * dx * stride;

        // Richardson extrapolation
        double power = 4.0;
        for (size_t j = 1; j <= i; ++j) {
            R[i][j] = (power * R[i][j - 1] - R[i - 1][j - 1]) / (power - 1.0);
            power *= 4.0;
        }

        // Check convergence
        if (i > 0 && std::abs(R[i][i] - R[i - 1][i - 1]) < tol) {
            return R[i][i];
        }
    }

    return R[k][k];
}

/**
 * @brief Romberg integration for matrix (column-wise, uniform spacing)
 *
 * Zero-copy operation: directly fills the provided result span without
 * internal allocation. Integrates each column independently using Romberg's
 * method.
 *
 * @param mat Input matrix (each column is a function, rows must be 2^k + 1)
 * @param result Output buffer for integral values (must have same size as
 * number of columns)
 * @param dx Spacing between rows
 * @param tol Tolerance for convergence (default: 1e-10)
 */
inline void romberg(const matrix::real_matrix_base &mat,
                    std::span<double> result,
                    double dx,
                    double tol = 1e-10) {
    if (mat.rows() < 3) {
        throw std::invalid_argument("Romberg: needs at least 3 rows");
    }

    if (result.size() != mat.cols()) {
        throw std::invalid_argument(
            "Romberg: result span must have same size as number of columns");
    }

    for (size_t j = 0; j < mat.cols(); ++j) {
        std::vector<double> col(mat.rows());
        for (size_t i = 0; i < mat.rows(); ++i) {
            col[i] = mat(i, j);
        }
        result[j] = romberg(col, dx, tol);
    }
}

/**
 * @brief Romberg integration for matrix (column-wise, uniform spacing)
 *
 * Convenience wrapper that allocates and returns a vector. For zero-copy
 * operations, use the void version with span parameter.
 *
 * Integrates each column independently using Romberg's method.
 *
 * @param mat Input matrix (each column is a function, rows must be 2^k + 1)
 * @param dx Spacing between rows
 * @param tol Tolerance for convergence (default: 1e-10)
 * @return Vector of integral values for each column
 */
inline std::vector<double>
romberg(const matrix::real_matrix_base &mat, double dx, double tol = 1e-10) {
    std::vector<double> result(mat.cols());
    romberg(mat, result, dx, tol);
    return result;
}

} // namespace msl::integral

#endif // MSL_ROMBERG_HPP
