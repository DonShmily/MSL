/**
**  MSL - Modern Scientific Library
**
**  Copyright 2025 - 2026, Dong Feiyue, All Rights Reserved.
**
** Project: MSL
** File: simpson.hpp
** -----
** File Created: Friday, 9th January 2026 14:58:10
** Author: Dong Feiyue (FeiyueDong@outlook.com)
** -----
** Last Modified: Thursday, 5th March 2026 14:58:00
** Modified By: Dong Feiyue (FeiyueDong@outlook.com)
*/

#ifndef MSL_SIMPSON_HPP
#define MSL_SIMPSON_HPP

#include <span>
#include <stdexcept>
#include <vector>

#include "matrix/real_matrix_base.hpp"

namespace msl::integral
{

// ============================================================================
// Total Integral (Simpson's Rule Integration)
// ============================================================================

/**
 * @brief Simpson's 1/3 rule (requires odd number of points)
 *
 * More accurate than trapezoidal for smooth functions
 *
 * @param y Function values (size must be odd)
 * @param dx Uniform spacing
 * @return Integral value
 */
inline double simpson(std::span<const double> y, double dx)
{
    if (y.size() < 3)
    {
        throw std::invalid_argument("Simpson: needs at least 3 points");
    }

    size_t n = y.size();

    // If even number of points, use trapezoidal for last segment
    bool use_trap_last = (n % 2 == 0);
    size_t n_simp = use_trap_last ? n - 1 : n;

    double sum = y[0] + y[n_simp - 1];

    // Odd indices (weight = 4)
    for (size_t i = 1; i < n_simp - 1; i += 2)
    {
        sum += 4.0 * y[i];
    }

    // Even indices (weight = 2)
    for (size_t i = 2; i < n_simp - 1; i += 2)
    {
        sum += 2.0 * y[i];
    }

    double result = sum * dx / 3.0;

    // Add trapezoidal correction for last segment if needed
    if (use_trap_last)
    {
        result += 0.5 * (y[n - 2] + y[n - 1]) * dx;
    }

    return result;
}

/**
 * @brief Cumulative Simpson's rule for vector (uniform spacing, with output
 * buffer)
 *
 * Computes cumulative integral: result[i] = ∫[0 to i] y dx using Simpson's
 * rule. For even indices, uses Simpson's rule; for odd indices, uses linear
 * interpolation.
 *
 * Zero-copy operation: directly fills the provided result span without internal
 * allocation.
 *
 * @param y Function values at equally spaced points
 * @param result Output buffer for cumulative integral values (must have same
 * size as y)
 * @param dx Spacing between points
 */
inline void
cumsimpson(std::span<const double> y, std::span<double> result, double dx)
{
    if (y.size() < 3)
    {
        throw std::invalid_argument("Simpson: needs at least 3 points");
    }

    if (result.size() != y.size())
    {
        throw std::invalid_argument(
            "Simpson: result span must have same size as input");
    }

    // First point
    result[0] = 0.0;

    // Use Simpson's rule for pairs of intervals
    for (size_t i = 2; i < y.size(); i += 2)
    {
        double simp = (y[i - 2] + 4.0 * y[i - 1] + y[i]) * dx / 3.0;
        result[i] = result[i - 2] + simp;

        // Linear interpolation for odd index
        if (i > 2)
        {
            result[i - 1] = 0.5 * (result[i - 2] + result[i]);
        }
        else
        {
            result[i - 1] = 0.5 * simp;
        }
    }

    // Handle last point if even number of points
    if (y.size() % 2 == 0)
    {
        size_t n = y.size();
        result[n - 1] = result[n - 2] + 0.5 * (y[n - 2] + y[n - 1]) * dx;
    }
}

/**
 * @brief Cumulative Simpson's rule for vector (uniform spacing)
 *
 * Computes cumulative integral: output[i] = ∫[0 to i] y dx using Simpson's
 * rule. For even indices, uses Simpson's rule; for odd indices, uses linear
 * interpolation.
 *
 * @param y Function values at equally spaced points
 * @param dx Spacing between points
 * @return Cumulative integral values
 */
inline std::vector<double> cumsimpson(std::span<const double> y, double dx)
{
    std::vector<double> result(y.size());
    cumsimpson(y, result, dx);
    return result;
}

/**
 * @brief Simpson's rule for matrix (column-wise, uniform spacing)
 *
 * Integrates each column independently using Simpson's rule.
 *
 * @param mat Input matrix (each column is a function)
 * @param dx Spacing between rows
 * @return Vector of integral values for each column
 */
inline void simpson(const matrix::real_matrix_base &mat,
                    std::span<double> result,
                    double dx)
{
    if (mat.rows() < 3)
    {
        throw std::invalid_argument("Simpson: needs at least 3 rows");
    }
    if (result.size() != mat.cols())
    {
        throw std::invalid_argument(
            "Simpson: result span must have same size as number of columns");
    }

    for (size_t j = 0; j < mat.cols(); ++j)
    {
        std::vector<double> col(mat.rows());
        for (size_t i = 0; i < mat.rows(); ++i)
        {
            col[i] = mat(i, j);
        }
        result[j] = simpson(col, dx);
    }
}

/**
 * @brief Simpson's rule for matrix (column-wise, uniform spacing)
 *
 * Integrates each column independently using Simpson's rule.
 *
 * @param mat Input matrix (each column is a function)
 * @param dx Spacing between rows
 * @return Vector of integral values for each column
 */
inline std::vector<double> simpson(const matrix::real_matrix_base &mat,
                                   double dx)
{
    std::vector<double> result(mat.cols());
    simpson(mat, result, dx);
    return result;
}

} // namespace msl::integral

#endif // MSL_SIMPSON_HPP